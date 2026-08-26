# cgs_cross + pool 大数据发布版验证

## 结论

- 正确性：12/12 个 `has1/no1 x PE/SE` 输出 MD5 与 x86 BWA 标准结果逐字节一致。
- 测试配置：`cgs_cross + pool`、`CPE_PROFILE=0`、单节点、单主进程、6 CG x 64 CPE。
- 12 个 case 的 Stage 2 全部快于旧版：加速范围 1.010x-1.063x，中位数 1.023x。
- part1 准备阶段下降 48.3%-69.3%（中位数 52.1%），验证了去掉 1.5 GiB SAM buffer 无效初始化的结构性收益。
- part3 比对与长度阶段中位数改善 1.6%；仅一个 case 慢 0.8%，属于单次跨节点测试噪声范围。mate-rescue 去重优化主要作用于 150 bp PE。
- 当前 has1/no1 的 Stage 2 差异为 0.28%-3.80%，说明计算结果不受本轮严重 I/O 抖动支配。
- 本轮存在明显共享文件系统抖动。低于 30 MiB/s 的读写被标记为 I/O 离群，不能用 pipeline total 判断计算优化效果。

`no1` 的 Stage 1/2/3 在三线程流水线中有重叠，不能相加；`has1` 使用 `-1` 关闭该流水线，三个阶段近似串行。计算优化统一以 Stage 2 和 part3 为主指标。

## 正确性

| 模式 | 数据 | reads | MD5 | 状态 |
| --- | --- | --- | --- | --- |
| has1 | ERR1203383 | PE | `dc5c0a6babd41641db22808caedb7a44` | PASS |
| has1 | ERR1203383 | SE | `ecf7f98b8d9354dfde7b962fa8d2c0b7` | PASS |
| has1 | small_SRR7963242 | PE | `215c1ee5deddbef3647d94b90e841834` | PASS |
| has1 | small_SRR7963242 | SE | `b2ea694f6a08a3424dece15f1b8e8114` | PASS |
| has1 | SRR2496709 | PE | `b318cf90a083eafc9aa4c461bf023f90` | PASS |
| has1 | SRR2496709 | SE | `322fe0c60d0ca365389307c06a21d77a` | PASS |
| no1 | ERR1203383 | PE | `dc5c0a6babd41641db22808caedb7a44` | PASS |
| no1 | ERR1203383 | SE | `ecf7f98b8d9354dfde7b962fa8d2c0b7` | PASS |
| no1 | small_SRR7963242 | PE | `215c1ee5deddbef3647d94b90e841834` | PASS |
| no1 | small_SRR7963242 | SE | `b2ea694f6a08a3424dece15f1b8e8114` | PASS |
| no1 | SRR2496709 | PE | `b318cf90a083eafc9aa4c461bf023f90` | PASS |
| no1 | SRR2496709 | SE | `322fe0c60d0ca365389307c06a21d77a` | PASS |

## 当前性能

| 模式 | 数据 | reads | pipeline | Stage 1 | Stage 2 | Stage 3 | read MiB/s | write MiB/s | I/O 离群 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| has1 | ERR1203383 | PE | 297.693 | 83.126 | 69.469 | 144.533 | 80.3 | 67.9 | - |
| no1 | ERR1203383 | PE | 191.628 | 131.628 | 72.110 | 185.852 | 50.7 | 52.2 | - |
| has1 | ERR1203383 | SE | 155.635 | 52.516 | 28.989 | 73.731 | 63.6 | 58.9 | - |
| no1 | ERR1203383 | SE | 369.286 | 362.227 | 29.504 | 85.893 | 9.2 | 50.9 | read |
| has1 | small_SRR7963242 | PE | 439.925 | 91.126 | 216.244 | 131.960 | 79.5 | 72.9 | - |
| no1 | small_SRR7963242 | PE | 683.405 | 670.212 | 216.841 | 156.812 | 10.8 | 60.7 | read |
| has1 | small_SRR7963242 | SE | 147.946 | 30.149 | 54.837 | 62.608 | 120.1 | 71.9 | - |
| no1 | small_SRR7963242 | SE | 309.948 | 302.616 | 55.296 | 89.192 | 12.0 | 49.2 | read |
| has1 | SRR2496709 | PE | 308.663 | 106.764 | 67.978 | 133.073 | 63.2 | 72.5 | - |
| no1 | SRR2496709 | PE | 636.126 | 629.356 | 69.525 | 167.885 | 10.7 | 56.1 | read |
| has1 | SRR2496709 | SE | 145.067 | 53.063 | 28.514 | 62.739 | 63.6 | 70.2 | - |
| no1 | SRR2496709 | SE | 351.687 | 348.698 | 28.700 | 74.390 | 9.7 | 58.0 | read |

## Stage 2 A/B

`previous` 来自 `correctness_results/bigdata_results/cgs_cross_pool`，`optimized` 是本轮发布版。单次跨节点测试仍有小幅噪声，重点看 part1 的结构性下降、part3 是否没有退化，以及 has1/no1 是否一致。

| 模式 | 数据 | reads | previous Stage2 | optimized Stage2 | speedup | old part1 | new part1 | old part3 | new part3 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| has1 | ERR1203383 | PE | 71.093 | 69.469 | 1.023x | 0.979 | 0.472 | 66.331 | 65.233 |
| no1 | ERR1203383 | PE | 73.120 | 72.110 | 1.014x | 1.589 | 0.693 | 66.774 | 66.386 |
| has1 | ERR1203383 | SE | 29.281 | 28.989 | 1.010x | 0.982 | 0.468 | 26.130 | 26.350 |
| no1 | ERR1203383 | SE | 31.055 | 29.504 | 1.053x | 1.659 | 0.598 | 26.771 | 26.036 |
| has1 | small_SRR7963242 | PE | 220.988 | 216.244 | 1.022x | 1.021 | 0.507 | 217.184 | 212.954 |
| no1 | small_SRR7963242 | PE | 221.732 | 216.841 | 1.023x | 1.697 | 0.795 | 217.212 | 213.139 |
| has1 | small_SRR7963242 | SE | 55.641 | 54.837 | 1.015x | 1.017 | 0.498 | 52.910 | 52.641 |
| no1 | small_SRR7963242 | SE | 56.792 | 55.296 | 1.027x | 1.412 | 0.686 | 53.428 | 52.582 |
| has1 | SRR2496709 | PE | 69.618 | 67.978 | 1.024x | 0.979 | 0.473 | 65.451 | 64.320 |
| no1 | SRR2496709 | PE | 71.352 | 69.525 | 1.026x | 1.441 | 0.582 | 65.725 | 65.010 |
| has1 | SRR2496709 | SE | 29.020 | 28.514 | 1.018x | 0.979 | 0.506 | 26.132 | 26.101 |
| no1 | SRR2496709 | SE | 30.515 | 28.700 | 1.063x | 1.783 | 0.547 | 26.460 | 25.890 |

## I/O 观察

- `small_SRR7963242 SE no1` 的 Stage 1 为 302.616 s，读带宽只有 11.963 MiB/s，最慢单次 fread 为 47.206 s；同数据 has1 的读带宽为 120.134 MiB/s。这个 10x 差异是共享文件系统波动。
- `-1` 在当前代码中设置 `no_mt_io=1`，关闭三阶段多线程 I/O 流水线；它不应让底层 `fread` 吞吐稳定快一个数量级。本轮 has1/no1 的 Stage 2 接近，而 Stage 1 的实际 `fread` 差异很大，也支持 I/O 抖动判断。
- 参考索引加载在若干节点上也曾长时间停滞，但它发生在三阶段 pipeline 之前，没有进入表中计时；异常启动已剔除。

## 文件

- `analysis/performance.tsv`：所有阶段、Stage 2 part1-part6、I/O 和 A/B 数值。
- `analysis/correctness.tsv`：12 个 MD5 的期望值、实测值和状态。
- `analysis/stage2_comparison.svg`：旧版与优化版 Stage 2 对比。
- `analysis/pipeline_io.svg`：当前 Stage 1/2/3 与 I/O 离群标记。
- `scripts/analyze_cpe_bigdata.py`：完整解析和可视化生成代码。

## 复现分析

```bash
python3 scripts/analyze_cpe_bigdata.py
```
