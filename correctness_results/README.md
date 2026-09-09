# Correctness Results

本目录保存 SWBWA 的正确性、性能和调度实验结果。常见目录结构为：

```text
<experiment>/<configuration>/<has1|no1>/<dataset>.*
```

- `*.run.log`：程序和作业运行日志，包含 pipeline/stage 计时。
- `*.sort.log`、`*.merge.log`：MPI split 输出合并或排序过程日志。
- `*.md5`、`*.concat.md5`、`md5_results.tsv`：标准化 SAM 的校验值及对照结果。
- `.status/`、`run.ok`、`md5.ok`：增量测试脚本的完成状态。
- `*.tsv`、`*.csv`：调度、分块、profiling 和消融实验的结构化数据。
- `*.svg`、`*.png`：由 `scripts/` 中分析程序生成的图表。
- `legacy_logs/`：早期按执行模式保存的历史日志归档。

`*.sam` 通常是较大的中间输出；正确性脚本在校验成功后会清理它们。目录名中的
日期或实验名称用于区分代码版本和测试目的，`has1`/`no1` 表示是否使用 `-1`。
