# MPI 动态输入粒度对照实验

## 目的

此前针对 ERR1203383 PE 的实现只拆分最后一个大 chunk，性能很好，但策略明显
依赖该数据的慢区间恰好位于文件末尾。这里把它改成更通用的两级调度：

- 文件主体继续使用 `-K` 控制的大 chunk；
- 从文件最后约 10% 开始，切换为 `big_chunk / 4` 的 micro-chunk；
- 10% 的起点向前取整到大 chunk 边界，因此细粒度区域不会小于文件的 10%；
- 所有名义边界仍由原有 FASTQ record 对齐逻辑修正。

另测一组关闭尾部细分、直接令 `-K=4 MiB` 的配置，用来判断能否把全文件都
改成小 chunk。

## 测试条件

- 数据：完整 ERR1203383 PE
- 输入：`../data/bwa_test_big_data/ERR1203383_{1,2}.fastq`
- 参考：`../data/GRCh38.d1.vd1.fa`
- MPI：单节点 6 ranks，每 rank 使用 1 个 CG
- 模式：dynamic input + split output
- 参数：`-t 1 -1 -I 170,80,500,1 -v 4`
- 精确索引：`MPI_EXACT_READ_INDEX=1`
- CPE profiling：关闭；CPE 主路径采用 stage2 的 part 3 累计时间

两组配置为：

1. `tail_10_percent`：默认 `-K=19,660,800` bytes，最后约 10% 使用
   `4,915,200`-byte micro-chunk。
2. `all_4_mib`：`SWBWA_MPI_TAIL_PERCENT=0`，`-K=4,194,304` bytes，
   全文件使用同一种小 chunk。

## 结果

| 策略 | 逻辑 chunks | 每 rank stage2 范围 | stage2 均值 | stage2 最大值 | 最大/均值 | CPE part 3 均值 | CPE part 3 最大值 | 单 rank 最大 RMA |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 尾部 10% micro-chunk | 232 | 57.850-71.196 s | 65.703 s | 71.196 s | 1.084 | 63.256 s | 68.615 s | 0.229 s |
| 全文件 4 MiB chunk | 835 | 71.979-82.453 s | 76.775 s | 82.453 s | 1.074 | 72.081 s | 77.835 s | 0.153 s |

全文件 4 MiB 的负载比值只从 1.084 小幅改善到 1.074，但 stage2 均值增加
16.9%，最大值增加 15.8%。主要损失来自 CPE 批次变小：part 3 均值增加
14.0%；part 1 的准备时间由 1.001 秒增加到 2.234 秒，CPE FASTQ 格式化和
SAM 生成的累计固定开销也大约翻倍。

RMA 不是这组退化的原因。尽管有效任务数由 232 增至 835，全文件 4 MiB 配置
的单 rank 最大 RMA 累计时间仍只有 0.153 秒，甚至低于尾部 10% 配置的
0.229 秒。小 chunk 的主要问题是 CPE kernel 吞吐和每批固定开销，而不是 MPI
ticket 本身。

## 与上一轮实验对照

| 策略 | 逻辑 chunks | stage2 均值 | stage2 最大值 | 最大/均值 |
| --- | ---: | ---: | ---: | ---: |
| 原始大 chunk，不细分 | 178 | 64.888 s | 98.629 s | 1.520 |
| 仅最后一个大 chunk 拆 6 份 | 183 | 64.991 s | 67.842 s | 1.044 |
| 最后约 10% 使用 1/4 chunk | 232 | 65.703 s | 71.196 s | 1.084 |
| 全文件 4 MiB chunk | 835 | 76.775 s | 82.453 s | 1.074 |

只拆最后一块仍是 ERR1203383 上的数据定制最优点。新的 10% 两级策略比它慢
约 4.9%，但不再假定慢区间只存在于最后一个大 chunk；相对完全不细分，最大
stage2 仍降低 27.8%，同时均值只增加约 1.3%。因此它更适合作为默认通用策略。

## 正确性

两组配置都输出 30,558,104 条 SAM records。对尾部 10% 配置的 6 个 split
文件做稳定 QNAME 排序后得到：

- 实际 MD5：`dc5c0a6babd41641db22808caedb7a44`
- 基准 MD5：`dc5c0a6babd41641db22808caedb7a44`

结果一致。远端测试生成的 SAM 已删除，仅保留日志和 MD5 摘要。

## 文件

- `summary.tsv`：核心指标汇总
- `tail10/run.log`：尾部 10% micro-chunk 完整运行日志
- `tail10/md5.log`、`tail10/sorted.md5`：排序正确性验证
- `chunk4m/run.log`：全文件 4 MiB chunk 完整运行日志
