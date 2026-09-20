# 2026-09-14 起 A/B 部署实验归档

A：每节点一个 `cgs_cross + pool` 进程；B：每节点六个 `single + system` 进程。
实验以各日志内的构建参数和 `CMD` 为准；原始日志中的旧账号/绝对路径作为溯源保留。

- `build_*`、`print-config_*`：构建与配置记录。
- `prep_data.log`：tiny/big4 数据生成；big4 是 small_SRR7963242 连续复制四份，非四个独立样本。
- `smoke_*`、`correct_*`、`md5_*`、`sort_*`：冒烟与正确性检查过程。
- `startup*`：多轮启动成本实验；原始 epoch、bjobs、stdout 保留，不能混作同一轮。
- `perf/`：discard 的 A/B 重复运行和汇总。缺失 rank 的计时行只能作为观察样本，
  不代表全 rank 最大值；完整哈希检查以 `scripts/check_discard_hash.py` 为准。

不包含 SAM、二进制或活动锁。此目录数据早于 40 KiB 预算及 2026-09-20 审查修订，
不能用它证明最新默认参数的性能或正确性。
