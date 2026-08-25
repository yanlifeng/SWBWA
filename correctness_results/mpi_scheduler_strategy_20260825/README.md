# MPI dynamic 调度策略对比

## 测试目的

本轮测试只观察 MPI dynamic 输入调度对计算负载的影响，不测试真实文件写入性能。
固定配置如下：

- 6 个 MPI rank，单节点，`q_share`
- `MPI_INPUT_MODE=dynamic`
- `MPI_EXACT_READ_INDEX=1`
- `OUTPUT_MODE=discard`
- `-1 -v 4`
- 每种正式策略、每个数据运行 3 次
- 测试数据为负载较差的 `ERR1203383 PE` 和 `SRR2496709 PE`

discard 模式不创建 SAM 文件。默认对每条 SAM 记录计算 64 位 hash，各 rank 的
sum/XOR 在运行结束后离线合并，用于验证不同调度策略没有丢失或修改输出。

## 策略

| 策略 | queue 选择 | 主体块 | 后 10% 块 | 最后细粒度区域 |
|---|---|---:|---:|---:|
| baseline | 成功后继续当前 queue | 19,660,800 B | 4,915,200 B | 无 |
| queue rotate | 成功后切换到下一个 queue | 19,660,800 B | 4,915,200 B | 无 |
| fine 1 wave | queue rotate | 19,660,800 B | 4,915,200 B | 约 1 个 rank wave，1,228,800 B/块 |
| fine 2 waves | queue rotate | 19,660,800 B | 4,915,200 B | 约 2 个 rank wave，1,228,800 B/块 |

对应的逻辑 chunk 数：

| 数据 | baseline / rotate | fine 1 wave | fine 2 waves |
|---|---:|---:|---:|
| ERR1203383 PE | 232 | 253 | 271 |
| SRR2496709 PE | 237 | 256 | 274 |

## 结果

下表中的“中位最大”是三次运行中，每次最慢 rank 的 stage2 时间再取中位数；
“最差最大”是三次运行里最慢的一次。它们比单独看 `max/min` 更接近整个计算阶段的
makespan。

| 策略 | 数据 | 中位最大 stage2 (s) | 最差最大 stage2 (s) | 中位 max/min | 中位最慢 chunk (s) | 中位最大 RMA (s) |
|---|---|---:|---:|---:|---:|---:|
| baseline | ERR1203383 | 69.778 | 70.212 | 1.163 | 12.172 | 0.200 |
| queue rotate | ERR1203383 | 69.905 | 69.967 | 1.241 | 12.186 | 0.995 |
| fine 1 wave | ERR1203383 | 68.276 | 69.589 | 1.048 | 5.042 | 1.176 |
| fine 2 waves | ERR1203383 | **66.471** | **66.903** | 1.060 | **5.014** | 1.214 |
| baseline | SRR2496709 | 70.364 | 71.112 | 1.147 | 9.414 | 0.320 |
| queue rotate | SRR2496709 | **65.283** | 67.745 | 1.121 | 9.442 | 1.044 |
| fine 1 wave | SRR2496709 | 65.975 | **66.759** | **1.087** | 9.392 | 1.100 |
| fine 2 waves | SRR2496709 | 67.358 | 67.889 | 1.100 | 9.412 | 1.234 |

所有 hash-enabled 运行的 calls、bytes、hash sum 和 hash XOR 都与 baseline 一致：

| 数据 | Calls | Bytes | Hash sum | Hash XOR |
|---|---:|---:|---|---|
| ERR1203383 PE | 30,549,360 | 8,029,871,401 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| SRR2496709 PE | 25,214,822 | 8,119,600,200 | `56efd396c8f85a1e` | `9edd66de273acc48` |

## 分析

1. queue rotate 不是无效修改。它没有切开 ERR 文件末尾的重块，因此 ERR 的最大
   stage2 基本不变；但它让其他 rank 更早消费不同 queue，对位于文件中段的
   SRR2496709 重块有效，中位最大 stage2 从 70.364 s 降至 65.283 s。

2. queue rotate 增加了远程 ticket 操作。最大 rank 的 scheduler RMA 从约
   0.2-0.3 s 增至约 1.0-1.2 s，但仍只占 stage2 的约 2%，不是主要瓶颈。

3. ERR1203383 的原始瓶颈是文件末尾约 12.2 s 的 5 MB chunk。将最终区域拆为
   1.25 MB 后，最慢 chunk 降至约 5.0 s。two-wave 三次最大 stage2 都稳定在
   66.4-66.9 s，明显比 one-wave 的 65.5-69.6 s 更稳定。

4. SRR2496709 的最慢 chunk 位于文件中段，约 9.4 s。细化文件末尾无法切开它，
   只会增加少量小块开销，因此 queue-only 或 fine 1 wave 在该数据上更快。

5. 从未知数据集上的稳健性考虑，fine 2 waves 更适合作为默认值。四种策略在两组
   数据、三次重复中的最差最大 stage2 分别为 71.112、69.967、69.589 和
   67.889 s，fine 2 waves 的最坏情况最好。`SWBWA_MPI_FINE_TAIL_WAVES=1`
   可以保留为偏吞吐、较少小块的实验配置，设置为 0 可退回 queue-only。

## discard hash 隔离测试

设置 `SWBWA_DISCARD_HASH=0` 后各运行一次。此时仍逐条调用输出接口并统计
calls/bytes，但不扫描 SAM 内容：

| 数据 | Stage2 范围 (s) | max/min | SAM output_write 范围 (s) |
|---|---:|---:|---:|
| ERR1203383 PE | 64.643-68.147 | 1.054 | 3.919-4.646 |
| SRR2496709 PE | 62.400-69.932 | 1.121 | 3.335-4.870 |

开启 hash 时，output_write 通常约为 6.5-8.8 s/rank。hash 内容扫描增加约
3-4 s/rank，并会改变下一次领取 chunk 的时刻，但关闭 hash 后 ERR 仍然均衡、
SRR 中段重块仍然存在，因此 hash 不是计算长尾的根因。由于无 hash 测试运行在
不同节点，绝对 stage2 时间不用于比较策略速度。

## 结论

当前建议保留以下默认组合：

```text
成功领取后轮转到下一个 queue
主体 20 MB
后 10% 使用 5 MB
最后 2 个 rank wave 使用 1.25 MB
```

这套组合同时改善了中段负载扩散和末尾不可再分重块，且额外 RMA 开销较小。
后续若继续优化 SRR2496709 一类“中段单块极重”的数据，需要让任意位置的 20 MB
chunk 可再分或可协作处理；继续只调尾部参数不会解决该问题。
