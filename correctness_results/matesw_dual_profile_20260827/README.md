# Mate-SW 双 alignment 可行性统计

## 目的

评估 `ksw_align2_matesw()` 是否值得改成一次同时处理两条 alignment、
从而利用 FP16 向量当前空闲的另外 16 条物理 lane。

统计只在 `CPE_PROFILE=1` 时启用。每个 CPE 在本地累计，单个 CPE
任务结束后才写回一次主存，避免在 `mem_matesw()` 热路径上产生共享写竞争。
`CPE_PROFILE=0` 时相关调用为空操作。

## 配置

- `EXEC_MODE=cgs_cross`
- `CPE_ALLOCATOR=pool`
- 当时使用 `FLOAT16_VECTOR=1`，等价于当前的
  `KSW_U8_MODE=float16_16`
- `CPE_PROFILE=1`
- 1 个节点、1 个 MPE 进程、6 个 CG、384 个 CPE
- `-t 1 -1 -v 4 -o /dev/null -I 170,80,500,1`

选取两个 read 长度和工作负载特征不同的 PE 子集：

- `ERR1203383`: 75 bp
- `SRR7963242`: 150 bp

## 结果

| 数据集 | `mem_matesw` 调用 | 0 候选 | 1 候选 | 2--4 候选 | 可配对候选 | 理想加速上限 | Stage 2 | CPE part 3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 PE | 3,107,306 | 2,426,930 | 680,376 | 0 | 0.00% | 1.000x | 3.694 s | 3.380 s |
| SRR7963242 PE | 15,743,913 | 6,371,538 | 9,372,375 | 0 | 0.00% | 1.000x | 44.778 s | 43.890 s |

这里的“候选”是实际进入 `ksw_align2_matesw()` 的方向；不是调用入口处
理论上的四个方向。现有 `skip[]`、参考区间和最小种子长度筛选完成后，两组
数据的每次 `mem_matesw()` 调用都最多只剩一个有效候选。

SRR7963242 本轮 Stage 1 遇到共享文件系统异常波动，读取带宽仅
`6.120 MiB/s`，但 Stage 2 是独立累计的活跃计算时间，不影响上述候选分布
和 CPE part 3 结论。

## 结论

不建议实现“单次 `mem_matesw()` 内双 alignment”。两组代表性数据均没有
任何可配对候选，理论加速上限就是 `1.000x`。要使用空闲的 16 条物理 lane，
必须跨不同 read/pair 批处理多个独立 `mem_matesw()` 调用；这需要改变当前
控制流、临时内存和结果回填方式，复杂度与风险明显更高，不能视为当前内核
的轻量优化。

## 验证

- `CPE_PROFILE=1` 的 `cgs_cross + pool` 构建和两组运行通过。
- `CPE_PROFILE=0` 的完整 `build_cross.sh` 构建通过，见
  `build_profile_off.log`。
- 原始运行输出见 `ERR1203383_PE.run.log` 和
  `SRR7963242_PE.run.log`。
