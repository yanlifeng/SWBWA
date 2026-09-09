# Correctness Results

本目录保存 SWBWA 的正确性、性能和调度实验结果。文件通常按以下层次组织：

```text
<experiment>/<configuration>/<has1|no1>/<dataset>.*
```

- `*.run.log`：一次运行的程序和作业日志，包含 pipeline/stage 计时。
- `*.sort.log`、`*.merge.log`：MPI split 输出合并或排序过程日志。
- `*.md5`、`*.concat.md5`、`md5_results.tsv`：标准化 SAM 的校验值及对照结果。
- `.status/`：增量正确性脚本的完成标记；`run.ok` 表示运行完成，`md5.ok` 表示校验通过。
- `*.tsv`、`*.csv`：调度、分块、profiling 和消融实验的结构化数据。
- `*.svg`、`*.png`：由 `scripts/` 中分析程序生成的图表。
- `legacy_logs/`：早期 `logs/{single,cgs,cross_cgs}` 的历史日志归档。

这里不保存 SAM 大文件；正确性脚本在校验成功后会清理中间 SAM。目录名中的日期表示实验日期，后缀如 `has1`/`no1` 表示是否使用 `-1` 运行参数。
