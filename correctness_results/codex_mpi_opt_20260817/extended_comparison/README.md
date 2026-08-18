# SWBWA MPI 动态输入优化：正确性与性能对比报告

## 1. 报告结论

本轮对比上一提交 `4fb1ff5` 与当前未提交优化版本，重点检查
`MPI_INPUT_MODE=dynamic` 下的 `split` 和 `single_unordered` 两种输出模式。

主要结果如下：

1. 两个版本、两种输出模式在全部 16 个子集正确性任务上均通过既有 MD5
   基准，没有发现输出回归。
2. 完整 `SRR2496709` PE 数据在四种版本/输出组合下得到相同的归一化
   MD5：`f90eb3fae1e001d710e3bdbce5d1ffbf`。该值证明四组结果彼此一致，
   但尚未与独立外部基准比较，不能单独作为完整数据绝对正确性的证明。
3. 当前版基本消除了旧版 `single_unordered` 特有的 Stage 2 失衡：完整数据
   的各 rank Stage 2 平均值由 `105.273 s` 降至 `65.711 s`，下降
   `37.6%`，并与当前 `split` 的 `66.238 s` 基本一致。
4. 当前版将完整数据 single 输出中的动态调度时间从
   `12.630-32.704 s/rank` 降到 `0.132-0.141 s/rank`；输出 offset 的
   `MPI_Win_flush` 最大累计时间从 `94.435 s` 降至 `0.130 s`。
5. Stage 2 已经不再是 single 输出相对 split 的主要差距。当前版完整数据
   pipeline 仍为 `177.741 s`，而 split 为 `139.451 s`，剩余差距主要来自
   多 rank 向同一文件执行 `pwrite` 时的文件系统吞吐。

## 2. “子集”到底指什么

这里的“子集”不是测试时从完整 FASTQ 临时 `head` 出来的数据，也不是每个
MPI rank 分到的某个 chunk。它指项目正确性脚本长期固定使用的、文件名以
`fastq_1` 结尾的测试输入：

| 测试名 | 实际 R1 输入 | PE 时使用的 R2 | R1 文件大小 | R1 FASTQ records | 动态 chunks |
| --- | --- | --- | ---: | ---: | ---: |
| `SRR2496709` 子集 | `../data/SRR2496709_1.fastq_1` | `../data/SRR2496709_2.fastq_1` | 583,901,384 B（556.852 MiB） | 2,101,235 | 30 |
| `ERR1203383` 子集 | `../data/ERR1203383_1.fastq_1` | `../data/ERR1203383_2.fastq_1` | 223,555,584 B（213.199 MiB） | 1,000,000 | 12 |
| `SRR2496709` 完整数据 | `../data/SRR2496709_1.fastq` | `../data/SRR2496709_2.fastq` | 3,536,060,312 B（3.293 GiB） | 12,607,411 | 180 |

表中的 record 数是 R1 FASTQ record 数：SE 时表示 read 数，PE 时表示 read
pair 数。PE 分块由 R1 计算 offset，并把相同 offset 用于 R2，因此本测试依赖
R1/R2 具有一致的 record 布局。文件大小也只列 R1，避免把一对 PE 文件误写成
单文件大小。

子集运行同时覆盖：

- `SRR2496709` PE：R1 子集 + R2 子集；
- `SRR2496709` SE：只使用 R1 子集；
- `ERR1203383` PE：R1 子集 + R2 子集；
- `ERR1203383` SE：只使用 R1 子集。

完整数据只运行 `SRR2496709` PE，没有运行完整数据 SE。

## 3. 对比版本和代码差异

### 3.1 上一版本

- Git 提交：`4fb1ff5`（`feat: progress MPI during CPE execution`）。
- single 输出使用 `MPI_File_write_at`。
- MPI progress 主要由执行流程中的显式探测承担。
- exact read index 保存累计 record offset，但领取 chunk 时仍可能重新打开文件、
  查找 FASTQ 边界。

### 3.2 当前优化版本

当前版本是测试时工作区内尚未提交的源码，涉及 `bwamem.c`、`fastmap.c`、
`swbwa_mpi.c`、`swbwa_mpi.h` 和 `swbwa_output.c`。主要变化为：

1. 启动独立 MPI progress 线程，约每 1 ms 调用一次 `MPI_Iprobe`，为 Sunway
   上依赖目标进程进入 MPI 的 RMA 请求提供 progress。
2. CPE 工作等待期间，如果 progress 线程已启动，主线程不再重复调用
   `MPI_Iprobe`。
3. `MPI_EXACT_READ_INDEX=1` 时同时缓存每个 chunk 的 FASTQ 字节边界和累计
   record offset；pipeline 中领取 chunk 后直接查表，不再重复扫描边界。
4. single 输出仍用 RMA 原子累加申请互不重叠的全局文件区间，但数据写入由
   `MPI_File_write_at` 改为 POSIX `pwrite`，避免大块文件 I/O 继续依赖 MPI
   progress 路径。
5. 增加 scheduler、chunk、输出和 progress 线程的详细统计。

## 4. 测试配置

### 4.1 编译和运行模式

| 项目 | 配置 |
| --- | --- |
| MPI ranks | 6 |
| MPI 输入模式 | `dynamic` |
| 精确 read 索引 | `MPI_EXACT_READ_INDEX=1` |
| 输出模式 | `split`、`single_unordered` |
| I/O 模式 | `no1`，即命令行不传 `-1`，使用默认三阶段 pipeline |
| SWBWA verbosity | `-v 4` |
| 比对线程 | `-t 1` |
| PE insert size | `-I 170,80,500,1` |
| 参考基因组 | `../data/GRCh38.d1.vd1.fa` |
| 动态 chunk 目标 | 19,660,800 B（18.75 MiB）/chunk |

比对任务使用的调度资源为：

```text
bsub -I -b -q q_share -N 1 -np 6 -cgsp 64 \
     -share_size 12000 -cache_size 128 -priv_size 16
```

完整数据的核心命令为：

```text
./SWBWA mem -v 4 -t 1 \
    -o <output.sam> -I 170,80,500,1 \
    ../data/GRCh38.d1.vd1.fa \
    ../data/SRR2496709_1.fastq \
    ../data/SRR2496709_2.fastq
```

### 4.2 测试矩阵

总计执行 20 个 SWBWA 比对任务：

- 2 个代码版本；
- 2 种输出模式；
- 每个组合包含 4 个子集任务（两个数据集各 PE/SE）和 1 个完整 PE 任务。

也就是 16 个有外部预期 MD5 的子集任务，加 4 个只比较跨版本一致性的完整
数据任务。

## 5. 正确性检查方法

MPI dynamic 会让不同 rank 动态领取 chunk，因此原始 SAM 的记录顺序不固定。
本轮没有直接对原始 SAM 做 MD5，而是先归一化顺序：

1. `split`：先按 rank 编号依次拼接 6 个输出文件，再调用
   `../fucking_sam_sort` 按 read ID 归一化排序。
2. `single_unordered`：直接对共享单文件调用同一个 sorter。
3. 对排序后的 SAM 执行 `md5sum`。
4. 子集结果与脚本内既有基准比较；完整结果只比较四组测试之间是否一致。
5. 校验完成后删除原始、拼接和归一化 SAM，只保留运行日志、排序日志和
   `md5_results.tsv`。

排序和 MD5 发生在 SWBWA 退出之后，不计入下面的 pipeline 性能。

## 6. 正确性结果

### 6.1 子集基准

| 数据集 | 模式 | 预期 MD5 | `4fb1ff5` split | 当前 split | `4fb1ff5` single | 当前 single |
| --- | --- | --- | --- | --- | --- | --- |
| `SRR2496709` | PE | `20f980ce3955d09e7c133ba26ddf2b77` | PASS | PASS | PASS | PASS |
| `SRR2496709` | SE | `2ca59c689707cc0aa6e0ce1fe30ae145` | PASS | PASS | PASS | PASS |
| `ERR1203383` | PE | `c2af4bf0b057d5125ce2a9d770f13741` | PASS | PASS | PASS | PASS |
| `ERR1203383` | SE | `473eec3972fbd651d8b911b9fb5c6e25` | PASS | PASS | PASS | PASS |

16 个子集任务全部 PASS。

### 6.2 完整 SRR2496709 PE

| 版本 | 输出模式 | 归一化 SAM records | MD5 |
| --- | --- | ---: | --- |
| `4fb1ff5` | split | 25,396,635 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| 当前版 | split | 25,396,635 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| `4fb1ff5` | single_unordered | 25,396,635 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| 当前版 | single_unordered | 25,396,635 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |

输入共有 12,607,411 个 read pairs，而 SAM 有 25,396,635 条记录。SAM 条数可以
高于输入 pair 数，因为一个 read 可能生成 supplementary/secondary 等额外记录。

四组完整结果完全一致，但由于没有提供完整数据的独立标准 MD5，这里只记录为
“跨实现一致”，不标记为外部基准 PASS。

## 7. 性能指标怎么读

日志中的时间含义如下：

- `Pipeline total`：从三阶段 pipeline 启动到结束的实际墙钟时间，是判断整次
  pipeline 快慢的主要指标。
- `Stage 1`：各 rank 累计用于申请缓冲区、领取输入范围和读取 FASTQ 的活跃时间。
- `Stage 2`：各 rank 累计用于 CPE 格式化、比对和生成 SAM 的活跃时间。
- `Stage 3`：各 rank 累计用于写 SAM、释放 batch 的活跃时间。
- `SAM output writes`：Stage 3 中调用统一输出接口的累计时间，包含必要的输出
  offset 申请及实际文件写入。
- `Scheduler`：pipeline 内动态领取 chunk、RMA ticket 和边界查询的累计时间。

三个 stage 会流水重叠，因此不能把 Stage 1、2、3 相加当作 pipeline 总时间。
表中 `min-max` 表示 6 个 MPI rank 的最小值和最大值，`avg` 是 6 个 rank 的
算术平均值。

另外，`MPI_EXACT_READ_INDEX=1` 的全文件预扫描发生在 pipeline 之前，本轮日志
没有为它提供独立总计时。因此下面比较的是 pipeline 性能，不是包含 exact-index
预扫描、SAM 排序和 MD5 的端到端作业时间。

## 8. 完整数据性能

### 8.1 Pipeline 和三个阶段

| 版本 | 输出 | Pipeline (s) | Stage 1 (s) | Stage 2 (s) | Stage 2 avg (s) | Stage 3 (s) | SAM writes (s) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `4fb1ff5` | split | 137.846 | 63.816-90.411 | 58.327-82.620 | 66.384 | 124.256-130.108 | 121.972-127.971 |
| 当前版 | split | 139.445-139.454 | 71.783-83.585 | 58.617-83.199 | 66.238 | 131.146-132.493 | 129.127-130.576 |
| `4fb1ff5` | single | 180.086 | 43.175-63.803 | 78.296-158.976 | 105.273 | 147.569-170.983 | 144.094-168.486 |
| 当前版 | single | 177.736-177.746 | 33.238-38.823 | 58.042-82.710 | 65.711 | 165.093-169.765 | 163.310-167.948 |

关键观察：

1. 当前 single 的 Stage 2 平均值比旧版下降 `37.6%`，最大值从
   `158.976 s` 降到 `82.710 s`，single 特有的巨大计算长尾基本消失。
2. 当前 single 与当前 split 的 Stage 2 平均值只差 `0.527 s`，说明两种输出
   模式下的比对计算效率已经接近。
3. 当前 split 与旧 split 的 Stage 2 平均值只差 `0.146 s`（约 `0.2%`），
   没有看到 progress 线程导致的可见计算回退。
4. 当前 split pipeline 比旧 split 慢约 `1.2%`，对应本轮 SAM 写入时间增加约
   3 秒。只有各一次运行，且文件系统写入本身波动明显，不能据此认定代码回归。
5. 当前 single pipeline 只比旧版缩短约 `1.3%`。计算和 RMA 等待已经明显改善，
   但同一共享文件的实际写入时间接近整个 pipeline，掩盖了 Stage 2 收益。

### 8.2 子集 Stage 2

以下均为 6 个 rank 的 `min-max`，单位为秒：

| 版本 | 输出 | SRR2496709 PE | SRR2496709 SE | ERR1203383 PE | ERR1203383 SE |
| --- | --- | ---: | ---: | ---: | ---: |
| `4fb1ff5` | split | 9.392-19.865 | 4.461-7.140 | 3.270-4.683 | 1.950-2.361 |
| 当前版 | split | 9.437-18.675 | 4.351-6.881 | 3.200-4.610 | 1.660-2.351 |
| `4fb1ff5` | single | 12.223-20.652 | 4.432-10.091 | 3.169-4.658 | 1.819-2.266 |
| 当前版 | single | 9.423-18.603 | 4.321-6.921 | 3.190-4.621 | 1.710-2.250 |

子集上也能看到相同趋势：当前 single 已回到与当前 split 接近的范围；数据较小
时，启动、同步和文件系统波动所占比例更高，因此不适合仅凭单个子集计算稳定的
加速比。

## 9. 动态调度与 chunk 负载

### 9.1 Scheduler 时间

| 版本 | 输出 | Scheduler/rank (s) | FASTQ 边界对齐/rank | 领取 chunks/rank |
| --- | --- | ---: | ---: | ---: |
| `4fb1ff5` | split | 3.962-21.065 | 3.782-13.898 s | 29-31 |
| 当前版 | split | 0.139-0.277 | 40-282 us | 30 |
| `4fb1ff5` | single | 12.630-32.704 | 1.790-6.772 s | 27-32 |
| 当前版 | single | 0.132-0.141 | 40-42 us | 30 |

完整数据共 180 个 chunk。当前版六个 rank 都正好领取 30 个，并且主要从各自
的 round-robin queue 领取。旧 single 中每个 rank 为 27-32 个，说明某些 rank
在输出/RMA progress 上停顿时，其他 rank 已经继续领取并“偷走”了更多工作。

缓存 exact chunk 边界后，pipeline 内边界查询从秒级降到几十微秒。注意这不是
免费消失：建立 exact index 时仍需在 pipeline 之前由 rank 0 扫描完整 R1，随后
广播边界和累计 record 信息；这里只是避免 pipeline 对同一边界重复做文件扫描。

### 9.2 Stage 2 的数据相关长尾

| 版本 | 输出 | 各 rank 最慢 chunk 的 Stage 2 (s) |
| --- | --- | ---: |
| `4fb1ff5` | split | 3.360-9.661 |
| 当前版 | split | 3.350-9.761 |
| `4fb1ff5` | single | 7.200-17.881 |
| 当前版 | single | 3.310-9.581 |

当前版消除了旧 single 中额外放大的 chunk 长尾，但约 `9.7 s` 的最慢 chunk 在
旧 split 中已经存在。这部分更像 read 本身的比对复杂度差异，不是 single 输出
新增的问题。按字节等分能让输入量接近，但不能保证每个 chunk 的比对成本相同。

## 10. RMA progress 效果

Sunway MPI 的 RMA 完成依赖目标 rank 定期进入 MPI。旧 single 中，rank 可能长时间
执行 CPE 计算或文件 I/O，导致其他 rank 的 `MPI_Win_flush` 等待目标 rank 提供
progress。当前版由独立线程持续 `MPI_Iprobe`，把 progress 与主 pipeline 解耦。

| 指标 | `4fb1ff5` single | 当前 single |
| --- | ---: | ---: |
| Scheduler/rank | 12.630-32.704 s | 0.132-0.141 s |
| 输出 `MPI_Win_flush(rank 0)` 最大累计时间 | 94.435 s | 0.130 s |
| Stage 2 avg/rank | 105.273 s | 65.711 s |
| Stage 2 max rank | 158.976 s | 82.710 s |

当前完整数据中 progress 线程统计为：

| 输出 | `MPI_Iprobe` 次数/rank | MPI 调用累计/rank | 平均单次 | 最大单次 |
| --- | ---: | ---: | ---: | ---: |
| split | 13,895-13,906 | 0.235-0.321 s | 16.9-23.1 us | 13.013 ms |
| single | 17,732-17,743 | 0.223-0.369 s | 12.6-20.8 us | 3.745 ms |

single 中没有一次 `MPI_Iprobe` 达到 10 ms；split 中只有一个 rank 记录到 2 次
不低于 10 ms 的调用。累计 MPI 调用时间不到 `0.4 s/rank`。线程调度本身的成本
没有单独计时，但当前 split 的 Stage 2 与旧版持平，说明本轮没有观察到明显的
计算侧开销。

## 11. 输出路径与剩余瓶颈

### 11.1 当前版实际写入

| 输出 | 写入方式 | 调用/rank | 系统调用累计/rank | 有效带宽/rank |
| --- | --- | ---: | ---: | ---: |
| split | 每 rank 独立文件，`write` | 21 | 125.309-126.996 s | 10.16-10.26 MiB/s |
| single | RMA 申请 offset，共享文件 `pwrite` | 21 | 164.106-165.043 s | 7.78-7.90 MiB/s |

single 的每次 flush 大约对应一个 64 MiB 输出缓冲区。RMA offset 申请已经很快，
真正占时间的是 6 个 rank 并发写同一个文件：每个 rank 的 `pwrite` 累计约
164-165 秒，带宽低于 split 的独立文件写入。

旧版 single 的 `MPI_File_write_at` 系统调用累计为 `69.750-138.514 s/rank`，
但它之外还存在大量 RMA/progress 等待，因此不能只拿该系统调用时间与当前
`pwrite` 直接比较。统一输出接口记录的完整 SAM write 时间才更接近 pipeline
实际承受的写出成本。

### 11.2 当前能下的结论

- progress 线程和边界缓存解决的是动态调度及 RMA progress 引起的计算失衡；
- `pwrite` 避免了大块文件 I/O 走 MPI 文件接口，但没有提升本轮共享单文件的
  实际文件系统吞吐；
- 若目标是当前机器上的总运行时间，`split` 仍明显优于 `single_unordered`；
- 若必须生成单文件，下一步应针对共享文件并发写入策略和文件系统行为优化，
  而不是继续压缩已经低于 0.15 秒的 scheduler/RMA offset 开销。

## 12. 局限与复现材料

1. 每个版本/输出组合只运行一次完整数据，I/O 数字不具备多次重复实验的置信区间。
2. exact-index 全文件预扫描未单独计时，也不在 pipeline 时间内。
3. 完整数据 MD5 没有独立外部基准，只能证明四组输出彼此一致。
4. 本轮只测试 1 个节点、6 ranks，不代表跨节点扩展性。
5. SAM 排序和 MD5 不计入性能结果。

目录结构：

```text
extended_comparison/
  previous_4fb1ff5/
    split/{subsets,full}/
    single_unordered/{subsets,full}/
  current_optimized/
    split/{subsets,full}/
    single_unordered/{subsets,full}/
  README.md
```

每个 `subsets/` 和 `full/` 目录保留：

- `*.run.log`：SWBWA 运行日志和各 rank 计时；
- `*.sort.log`：归一化排序记录；
- `md5_results.tsv`：最终 MD5；
- `build.log`：对应子集批次的构建日志（存在时）。

所有远端和本地 SAM 中间文件均已删除。
