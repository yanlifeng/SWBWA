# MPI 尾部微块大数据测试

## 测试目标

验证动态 MPI 输入调度的“主体大块 + 最后 10% 微块”策略：

- 主体 chunk：`19,660,800` B（18.75 MiB）
- 尾部 chunk：`4,915,200` B（4.6875 MiB，主体的 1/4）
- 尾部范围：输入文件最后 10%
- MPI rank：6
- 输入模式：`dynamic`
- 输出模式：`split`、`single_unordered`
- 运行模式：`no1`
- 正确性模式：`MPI_EXACT_READ_INDEX=1`

本轮没有测试 `has1` 和 `static`。校验阶段会合并并按 read ID 排序 SAM，计算归一化 MD5；排序耗时不计入 SWBWA pipeline。

## 正确性

三个数据集的 PE/SE 在两种输出模式下全部通过，共 **12/12 PASS**。

| 数据 | 模式 | 期望及实际 MD5 |
|---|---|---|
| ERR1203383 | PE | `dc5c0a6babd41641db22808caedb7a44` |
| ERR1203383 | SE | `ecf7f98b8d9354dfde7b962fa8d2c0b7` |
| small_SRR7963242 | PE | `215c1ee5deddbef3647d94b90e841834` |
| small_SRR7963242 | SE | `b2ea694f6a08a3424dece15f1b8e8114` |
| SRR2496709 | PE | `b318cf90a083eafc9aa4c461bf023f90` |
| SRR2496709 | SE | `322fe0c60d0ca365389307c06a21d77a` |

split 与 single_unordered 对同一输入得到相同 MD5。成功后产生的 SAM 已删除，目录中只保留日志、MD5 和分析结果。

## 本轮性能

时间为 6 个 rank 的 `median [min, max]`，单位为秒。流水线各 stage 会重叠，不能相加。`CPE align` 是 stage2 中的 CPE alignment/SAM-length pass。

| 输出 | 数据 | 模式 | pipeline | stage1 | stage2 | CPE align | stage3 | 读带宽中位数 | 写带宽中位数 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| split | ERR1203383 | PE | 148.6 [148.6, 148.6] | 76.5 [72.6, 78.4] | 65.5 [58.4, 78.0] | 62.6 [55.6, 75.1] | 121.5 [117.4, 122.3] | 19.1 MiB/s | 11.2 MiB/s |
| single | ERR1203383 | PE | 176.3 [176.3, 176.3] | 34.1 [31.9, 36.5] | 67.0 [58.6, 73.4] | 64.2 [55.8, 70.5] | 156.4 [153.3, 157.1] | 39.3 MiB/s | 8.5 MiB/s |
| split | ERR1203383 | SE | 67.3 [67.3, 67.3] | 41.5 [40.2, 42.3] | 26.4 [24.9, 28.0] | 24.2 [22.6, 25.8] | 55.4 [55.0, 57.3] | 18.3 MiB/s | 10.7 MiB/s |
| single | ERR1203383 | SE | 83.1 [83.1, 83.1] | 18.9 [17.1, 20.5] | 26.9 [23.9, 30.0] | 24.4 [21.7, 27.4] | 71.2 [69.9, 72.3] | 37.4 MiB/s | 8.0 MiB/s |
| split | small_SRR7963242 | PE | 293.2 [293.2, 293.2] | 62.1 [56.7, 69.1] | 287.6 [286.3, 289.2] | 284.4 [283.5, 286.2] | 118.0 [85.9, 122.0] | 24.1 MiB/s | 11.9 MiB/s |
| single | small_SRR7963242 | PE | 277.4 [277.4, 277.4] | 55.8 [48.7, 65.0] | 271.4 [270.7, 272.3] | 268.2 [267.7, 269.5] | 137.2 [103.1, 190.1] | 27.2 MiB/s | 9.3 MiB/s |
| split | small_SRR7963242 | SE | 75.7 [75.7, 75.7] | 42.9 [39.3, 47.5] | 53.5 [52.9, 54.5] | 51.2 [50.3, 52.0] | 61.6 [60.8, 62.0] | 18.5 MiB/s | 10.2 MiB/s |
| single | small_SRR7963242 | SE | 92.1 [92.1, 92.1] | 17.5 [16.0, 18.9] | 53.3 [51.9, 56.2] | 50.8 [49.5, 53.5] | 76.2 [74.8, 79.2] | 43.8 MiB/s | 8.1 MiB/s |
| split | SRR2496709 | PE | 139.9 [139.9, 139.9] | 82.3 [74.2, 92.0] | 65.1 [59.8, 79.9] | 62.2 [57.0, 77.2] | 131.2 [129.4, 131.6] | 17.9 MiB/s | 10.3 MiB/s |
| single | SRR2496709 | PE | 176.4 [176.4, 176.4] | 36.6 [34.3, 38.5] | 64.2 [59.2, 80.2] | 61.4 [56.4, 77.4] | 164.9 [162.7, 166.5] | 37.7 MiB/s | 8.0 MiB/s |
| split | SRR2496709 | SE | 64.9 [64.9, 64.9] | 40.5 [39.9, 43.1] | 26.6 [24.9, 28.5] | 24.3 [22.8, 26.3] | 57.3 [52.7, 57.5] | 18.3 MiB/s | 10.8 MiB/s |
| single | SRR2496709 | SE | 82.7 [82.7, 82.7] | 18.7 [17.5, 20.2] | 27.1 [24.1, 29.5] | 24.7 [21.9, 27.1] | 71.4 [69.2, 72.8] | 37.9 MiB/s | 8.3 MiB/s |

## 与原动态大块版本对比

对照数据来自 `correctness_results/bigdata_mpi_lwpf_off_20260823` 中相同的 `dynamic + no1` case。

### ERR1203383 PE

这是尾部微块的目标工作负载，收益明确：

| 输出 | pipeline | 最慢 rank 的 stage2 | stage2 max/median |
|---|---:|---:|---:|
| split 原版 | 161.6 s | 99.2 s | 1.665x |
| split 尾部微块 | 148.6 s (-8.1%) | 78.0 s (-21.4%) | 1.191x |
| single 原版 | 188.2 s | 99.1 s | 1.672x |
| single 尾部微块 | 176.3 s (-6.3%) | 73.4 s (-25.9%) | 1.095x |

原版最慢的末尾大块约 42 秒。当前末尾重数据被拆开后，最慢微块约 12.2 秒；仍能看到 chunk 229/230/231 较重，但已经分散到多个 rank。

### 其他数据

- `ERR1203383 SE` 的 stage2 最大值下降 4.3%-10.2%，pipeline 基本持平或小幅改善。
- `SRR2496709 PE` 的 stage2 最大值下降 2.9%-4.9%，pipeline 波动为 +0.8%-1.6%。它的慢块主要在 chunk 100 等文件中部，尾部策略不会专门拆这些块。
- `small_SRR7963242 PE` 本身已经高度均衡，stage2 max/median 从约 1.01x 变为 1.00x，性能基本不变。
- 其余 SE case 的 pipeline 变化约在 -1.0% 到 +3.1%，属于本轮 I/O 和运行波动范围。

尾部微块增加了 ticket 数，但每个 rank 的 scheduler RMA ticket path 最多约 0.48 秒；相对几十到数百秒的 pipeline 很小。它保住了主体大 chunk 对 CPE kernel 的吞吐，同时明显缓解了 ERR1203383 PE 的尾部负载。

## single_unordered 输出诊断

single 相对 split 的 stage2 基本一致，差异主要来自输出：

- single 的 `MPI_Win_flush` 原子 offset 路径，每个 case 六个 rank 合计只有 0.34-0.69 秒，最慢 rank 不超过 0.16 秒。
- 真正耗时的是并行 `pwrite`：single 中位带宽只有 8.0-9.3 MiB/s/rank；split 为 10.2-11.9 MiB/s/rank。
- 除 `small_SRR7963242 PE` 这次受 stage2 运行波动影响外，single pipeline 比 split 慢约 18.7%-27.4%。

因此，这轮结果不支持“single 慢在 RMA offset 原子操作”；主要限制仍是共享文件系统上多个 rank 对同一文件做随机位置写入的吞吐。

## 结论

当前“最后 10% 使用 1/4 chunk”是比全局缩小 `-K` 更合适的折中：

1. 12 个正确性 case 全部通过。
2. ERR1203383 PE 的尾部失衡显著改善，pipeline 缩短 6%-8%。
3. 对已经均衡的数据，性能基本保持，最差观察到约 3% 波动。
4. 调度 RMA 开销仍然很小，没有抵消微块收益。
5. single_unordered 的后续优化重点应放在文件写入布局和共享文件系统行为，而不是继续压缩 offset RMA。

详细数据：

- `md5_results.tsv`：正确性结果
- `baseline_comparison.tsv`：与原动态大块版本的逐 case 对照
- `analysis/mpi_rank_timings.tsv`：rank 级原始指标
- `analysis/mpi_case_summary.tsv`：case 汇总
- `analysis/mpi_stage2_balance.svg`：stage2 rank 均衡图
- `analysis/mpi_outlier_matrix.svg`：stage/I/O 热力图
