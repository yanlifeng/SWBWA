# cgs_cross + pool CPE 热点分析与低风险优化

## 1. 目标与实验边界

本轮工作回到单节点、单主进程的 `cgs_cross + pool` 配置，使用 LWPF3 将
Stage2 中的 CPE 路径细分到 kernel 级别，并在保持 SAM 输出逐字节一致的前提下
进行低风险优化。

实验配置：

- 执行模式：`EXEC_MODE=cgs_cross`
- CPE 分配器：`CPE_ALLOCATOR=pool`
- MPI：关闭
- LWPF3：开启，只采样 CG5，输出 64 个 CPE 的 avg/min/max
- 作业资源：单节点、单主进程、6 个 CG，每个 CG 64 个 CPE
- BWA 参数：`-t 1 -1 -I 170,80,500,1`
- 代表数据：
  - `ERR1203383_*fastq_1`：75 bp，短 read 路径代表
  - `SRR7963242_*fastq_1`：150 bp，mate rescue 重负载代表
  - `SRR2496709_*fastq_1`：用于补充 PE 正确性与 Stage2 验证

LWPF 会给热点路径增加计数开销，因此图表适合比较同一插桩布局下的相对变化，
不应直接当作关闭 `CPE_PROFILE` 后的发布版绝对性能。每个 A/B case 当前只运行一
次，结论定位为方向性结果，而不是统计显著性结论。

## 2. 新增的父子 kernel 层次

```text
WORKER_ALIGNMENT
├── MEM_CHAIN
│   ├── MEM_CHAIN_COLLECT
│   └── MEM_CHAIN_BUILD
├── CHAIN_FILTER
├── CHAIN_EXTENSION
│   └── CHAIN_EXTENSION_DP
├── ALIGNMENT_FINALIZE
└── SAM_FORMAT
    ├── MATE_RESCUE
    │   ├── MATE_REF_FETCH
    │   ├── MATE_KSW_ALIGN
    │   │   ├── KSW_QUERY_INIT_FORWARD
    │   │   ├── KSW_DP_FORWARD
    │   │   ├── KSW_QUERY_INIT_REVERSE
    │   │   └── KSW_DP_REVERSE
    │   └── MATE_DEDUP
    └── PAIRING

SAM_COPY  # 输出复制阶段，位于 WORKER_ALIGNMENT 之外
```

父 kernel 包含子 kernel 的时间，父子计数不能相加。可视化中的连线和缩进表达
这一包含关系。

## 3. 热点结论

### ERR1203383 PE（75 bp）

基线中每个被采样 CPE 的平均周期数：

| 路径 | 平均周期 | 占直接父路径比例 |
| --- | ---: | ---: |
| `WORKER_ALIGNMENT` | 8.716 G | 100% |
| `MEM_CHAIN` | 4.703 G | 54.0% |
| `MEM_CHAIN_COLLECT` | 3.052 G | `MEM_CHAIN` 的 64.9% |
| `MEM_CHAIN_BUILD` | 1.406 G | `MEM_CHAIN` 的 29.9% |
| `CHAIN_EXTENSION` | 0.893 G | 10.2% |
| `CHAIN_EXTENSION_DP` | 0.451 G | `CHAIN_EXTENSION` 的 50.5% |
| `SAM_FORMAT` | 2.644 G | 30.3% |
| `MATE_RESCUE` | 1.657 G | `SAM_FORMAT` 的 62.7% |
| `MATE_KSW_ALIGN` | 1.044 G | `MATE_RESCUE` 的 63.0% |

短 read 的首要热点是 `mem_chain()`，其中 seed interval 收集比后续链构建更重。

### SRR7963242 PE（150 bp）

使用包含 `MATE_DEDUP` 计数器的 `baseline_dedup` 作为匹配基线：

| 路径 | 平均周期 | 占直接父路径比例 |
| --- | ---: | ---: |
| `WORKER_ALIGNMENT` | 113 G | 100% |
| `MEM_CHAIN` | 25 G | 22.1% |
| `CHAIN_EXTENSION` | 22 G | 19.5% |
| `CHAIN_EXTENSION_DP` | 18 G | `CHAIN_EXTENSION` 的 81.8% |
| `SAM_FORMAT` | 60 G | 53.1% |
| `MATE_RESCUE` | 57 G | `SAM_FORMAT` 的 95.0% |
| `MATE_KSW_ALIGN` | 35 G | `MATE_RESCUE` 的 61.4% |
| `KSW_DP_FORWARD` | 26 G | `MATE_KSW_ALIGN` 的 74.3% |
| `KSW_DP_REVERSE` | 7.155 G | `MATE_KSW_ALIGN` 的 20.4% |
| `MATE_DEDUP` | 21 G | `MATE_RESCUE` 的 36.8% |
| `MATE_REF_FETCH` | 0.325 G | `MATE_RESCUE` 的 0.6% |

150 bp PE 的决定性热点是 mate rescue。参考序列获取很轻，主要成本是 KSW 动态规
划和每轮候选结果去重，KSW 内部又以前向 DP 为主。

## 4. 实施的优化

### 4.1 mate rescue 只做一次去重

原实现每处理一个方向后，只要之前任一方向产生过候选结果，就会再次调用
`mem_sort_dedup_patch()`。现在记录是否新增过候选，四个方向处理完成后统一去重
一次。候选集合不变，排序与去重函数不变，仅去掉对逐步中间集合的重复处理。

在 `SRR7963242 PE` 上：

- `MATE_DEDUP`：约 21 G -> 19 G 平均周期，下降约 9.5%
- Stage2 part3：54.640 s -> 53.830 s，下降约 1.5%
- Stage2 总计：55.968 s -> 54.637 s，下降约 2.4%

### 4.2 去掉静态 SAM buffer 的无效初始化

首次创建三个 512 MiB SAM buffer 时，原代码先将总计 1.5 GiB 内存全部填充为
`'A'`。随后 CPE 会覆盖所有实际输出字节，主核也会设置记录终止位置，因此这次
整块初始化不参与任何结果。

去掉该初始化后，Stage2 part1 的一次性准备开销：

- `ERR1203383 PE`：0.614 s -> 0.097 s，下降 84.2%
- `SRR7963242 PE`：0.703 s -> 0.187 s，下降 73.4%

这是主核侧最明确且与数据内容无关的收益。

## 5. Stage2 A/B 结果

| 数据 | 版本 | Stage2 | part1 准备 | part2 格式化 | part3 比对与长度 | part4 SAM 切片 | part5 SAM 生成 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 PE | baseline | 5.088 s | 0.614 s | 0.033 s | 4.250 s | 0.151 s | 0.039 s |
| ERR1203383 PE | optimized | 4.295 s | 0.097 s | 0.028 s | 3.980 s | 0.151 s | 0.039 s |
| SRR7963242 PE | baseline_dedup | 55.968 s | 0.703 s | 0.129 s | 54.640 s | 0.380 s | 0.110 s |
| SRR7963242 PE | optimized | 54.637 s | 0.187 s | 0.129 s | 53.830 s | 0.380 s | 0.110 s |

整体变化：

- `ERR1203383 PE` Stage2 减少 15.6%，约 1.18x
- `SRR7963242 PE` Stage2 减少 2.4%，约 1.02x

ERR 的 part3 也出现约 6.4% 变化，但单次测试无法排除节点状态和计数扰动，不能
把这一部分全部归因于代码修改。SRR 的匹配插桩基线更能反映去重优化本身。

## 6. 正确性

优化后对三个代表数据进行了 PE 验证，并补测两个 SE 数据。所有输出均与项目的
既有标准 SAM MD5 完全一致：

| 数据 | 模式 | MD5 | 结果 |
| --- | --- | --- | --- |
| ERR1203383 | PE | `c2af4bf0b057d5125ce2a9d770f13741` | PASS |
| ERR1203383 | SE | `473eec3972fbd651d8b911b9fb5c6e25` | PASS |
| SRR2496709 | PE | `20f980ce3955d09e7c133ba26ddf2b77` | PASS |
| SRR7963242 | PE | `a799dc7268f389120ec92820e04b5118` | PASS |
| SRR7963242 | SE | `0acfb1f46abd9fed3862c28a33e26da4` | PASS |

## 7. 结果文件

- `analysis/cpe_kernel_hierarchy.svg`：CPE 父子 kernel 周期对比
- `analysis/stage2_breakdown.svg`：Stage2 主核六部分对比
- `analysis/kernel_profile.tsv`：所有 LWPF avg/min/max 原始解析结果
- `analysis/stage2_profile.tsv`：Stage1/2/3 与 Stage2 六部分计时
- `analysis/correctness.tsv`：MD5 检查结果
- `scripts/analyze_cpe_profile.py`：位于项目根目录下的解析和 SVG 生成程序

PNG 是由相同 SVG 转换得到的便览文件，SVG/TSV 是可复现的源结果。

## 8. 复现方法

在神威项目目录中编译：

```bash
bash build_cross.sh \
  EXEC_MODE=cgs_cross \
  CPE_ALLOCATOR=pool \
  CPE_PROFILE=1 \
  CPE_PROFILE_CG=5 \
  USE_MPI=0
```

单节点、单主进程运行示例：

```bash
bsub -I -b -q q_share -n 1 -cgsp 64 \
  -share_size 2000 -mpecg 6 -xmalloc -cross_size 42000 \
  -cache_size 128 -priv_size 16 \
  ./SWBWA mem -v 3 -t 1 -1 \
  -I 170,80,500,1 \
  -o result.sam \
  ../data/GRCh38.d1.vd1.fa \
  ../data/SRR7963242_1.fastq_1 \
  ../data/SRR7963242_2.fastq_1
```

在项目根目录重新生成 TSV 和 SVG：

```bash
python3 scripts/analyze_cpe_profile.py \
  --root correctness_results/cpe_hotspot_opt_20260825
```

如系统安装了 librsvg，可生成 PNG：

```bash
rsvg-convert -o correctness_results/cpe_hotspot_opt_20260825/analysis/cpe_kernel_hierarchy.png \
  correctness_results/cpe_hotspot_opt_20260825/analysis/cpe_kernel_hierarchy.svg
rsvg-convert -o correctness_results/cpe_hotspot_opt_20260825/analysis/stage2_breakdown.png \
  correctness_results/cpe_hotspot_opt_20260825/analysis/stage2_breakdown.svg
```

## 9. 下一步方向

本轮没有冒险改动 KSW 算法或评分语义。若继续优化单进程版本，优先级应为：

1. 研究 `MATE_KSW_ALIGN` 的前向 DP，尤其是 150 bp PE 路径；
2. 分析 `mem_collect_intv()` 的访存和分支行为，面向 75 bp 短 read；
3. 在 `CPE_PROFILE=0` 下做多次交错 A/B，给出发布版的中位数和波动范围；
4. 任何 KSW 或 seed 搜索改动都继续用 PE/SE 五组标准 MD5 做回归。
