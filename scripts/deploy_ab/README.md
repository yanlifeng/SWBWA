# A/B 部署实验工具

A 为每节点一个 cgs_cross + pool 进程，B 为每节点六个 single + system 进程。
历史四节点实验归档于 `correctness_results/deploy_ab_20260914` 和
`correctness_results/deploy_tailoff_20260916`，并非当前未测代码的回归证明。

- `run_discard_ab.sh`：ABBAAB 顺序运行，失败或缺少完整哈希则停止。预先准备
  `$WORK_ROOT/runs/{A,B}_discard` 二进制及匹配的 cross 重定位文件。
- `parse_perf_ab.py`：汇总上述日志，保留已观察样本的 min/max/median。
- `parse_tailoff.py`：汇总 small/big4、A/B、各次重复的尾部关闭实验。
- `verify_big4.py DATA_DIRECTORY`：检查大小和各段首尾，**不是全文件字节校验**。

运行须显式提供 `WORK_ROOT`、`DATA_ROOT`，默认 `NODES=1`。历史实验使用 `NODES=4`，
重跑前必须另外取得资源授权。驱动和 `../startup_driver3.sh` 共用 `$WORK_ROOT/driver.lock`，
禁止多个控制器同时运行；断连/状态未知后先查 `bjobs`。不会自动覆盖已有汇总。
每次运行前后查看 `gfsquota`，不要提交超出配额的作业。

性能解析允许缺失的计时行，仅给出观察值，不能把观察到的最大值当成全体 rank 的最大值。
指纹是否有效使用公共完整性检查器判断；缺 rank、重复 rank、关闭哈希或截短哈希均不能证明
正确性。过去日志未记录 `hash_prefix_bytes` 的，只能结合对应构建日志确定其语义。
