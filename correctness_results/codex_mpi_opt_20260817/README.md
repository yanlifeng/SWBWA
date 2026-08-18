# MPI 动态输入与输出优化报告

## 测试配置

- 数据集：`SRR7963242`，分别测试 PE 和 SE
- 计算资源：一个申威节点，六个 MPI rank
- 输入模式：动态调度，`MPI_EXACT_READ_INDEX=1`，不使用 `-1`
- 输出模式：split 或 single unordered
- 性能分析：使用 `-v 4`
- 正确性检查：对每个保留方案计算标准化 SAM 的 MD5

流水线的不同阶段会相互重叠，因此各阶段耗时不能直接相加。

## 测试结果

| 方案 | PE Pipeline | PE Stage 2 范围 | SE Pipeline | SE Stage 2 范围 |
| --- | ---: | ---: | ---: | ---: |
| 原始 split | 74.61 s | 67.42-69.90 s | 20.59 s | 13.38-14.38 s |
| split + progress 线程 | 73.73 s | 65.58-70.63 s | 20.04 s | 13.34-14.21 s |
| 原始 single，使用 MPI-IO | 75.80 s | 66.50-70.46 s | 27.67 s | 17.56-21.84 s |
| single + `pwrite` | 73.67 s | 63.39-67.23 s | 26.90 s | 13.54-14.10 s |
| single + `pwrite` + progress 线程 | 74.17 s | 63.69-67.98 s | 27.13 s | 13.55-14.08 s |
| 最终 single + 保留 exact 边界 | 73.05 s | 62.56-67.81 s | 27.43 s | 13.39-14.18 s |

最终方案中 SE 的 Pipeline 时间主要受共享文件系统的运行波动影响。更稳定、
更有参考价值的改进指标是 Stage 2 的时间范围和 RMA 操作耗时。

## 分析结论

1. `MPI_File_write_at` 会干扰申威 MPI 的目标端 progress。在原始 single
   输出的 PE 测试中，动态调度器的 RMA 时间最高达到 19.67 秒。改为使用 MPI
   RMA 申请互不重叠的文件区间，再通过 POSIX `pwrite` 写入后，RMA 时间下降到
   0.20-0.41 秒。
2. 独立的 progress 线程在主线程等待从核计算期间调用 `MPI_Iprobe`，能够使输入
   和输出 RMA 的执行时间更加稳定。最终方案中，每个 rank 的调度器 RMA 时间为
   0.13-0.26 秒。
3. exact-index 初始化阶段原本已经计算了全部对齐后的 FASTQ chunk 边界，但随后
   丢弃了这些结果。现在保留并广播这些边界后，流水线中的边界查找时间从每个
   rank 的 0.8-3.8 秒下降到约 10-12 微秒。
4. 将输出缓冲区从 64 MiB 缩小到 16 MiB 后，PE 每个 rank 的 flush 次数由约
   6 次增加到 21-22 次，部分 rank 累积的 `pwrite` 时间上升到约 50 秒。因此
   放弃该方案，继续使用默认的 64 MiB 缓冲区。
5. split 输出原本就不存在 MPI-IO 与 RMA 之间的相互影响。progress 线程可以消除
   偶发的调度器等待，但不会显著改变整体吞吐率。目前运行波动的主要来源已经变为
   共享 FASTQ 和 SAM 文件的存储带宽。

## 正确性与编译验证

- PE MD5：`a799dc7268f389120ec92820e04b5118`
- SE MD5：`0acfb1f46abd9fed3862c28a33e26da4`
- 最终 single 输出的 PE 和 SE 正确性检查均通过，生成的 SAM 文件已删除。
- split 输出的 PE 和 SE 正确性检查也已通过。
- 非 MPI、MPI exact-index 和 MPI inexact-index 三种配置均编译通过。
- 编译过程中唯一的链接警告是 `kopen.c` 原有的静态 glibc `getaddrinfo` 警告。

## 后续验证

当前正确性和性能测试使用一个申威节点、六个 MPI rank。single unordered 模式
依赖共享文件系统正确支持多个进程向互不重叠的文件区间执行 `pwrite`。该语义已在
单节点上验证，后续还需要补充多节点环境下的正确性和吞吐率测试。
