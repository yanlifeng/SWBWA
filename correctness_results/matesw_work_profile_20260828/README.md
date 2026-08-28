# Same-PE mate-SW 长度与 DP 工作量 profiling

## 目的

上一轮计数表明，150 bp PE 数据中同一 read pair 的两个 mate-rescue 方向
可覆盖 57.61% 的 KSW 候选。这里继续测量这些候选对的实际执行形态，回答：

1. 两个 alignment 的 `qlen/tlen` 和 KSW kernel 维度是否兼容。
2. forward/reverse DP 与 Lazy-F 的实际 SIMD 迭代数是否接近。
3. 假设用 32 个物理 lane 同时运行两个独立 16-lane 状态，理论工作量能减少多少。

## 配置

- `EXEC_MODE=cgs_cross`, `CPE_ALLOCATOR=pool`
- `CPE_PROFILE=1`, sampled CG 5
- `KSW_U8_MODE=int32_16`
- `-t 1 -1 -v 4 -I 170,80,500,1`
- 输出写入 `/dev/null`
- 75 bp：`ERR1203383` PE 子集
- 150 bp：`SRR7963242` PE 子集

## 结果

| 数据 | 有序候选对 | forward qlen 相同 | forward tlen 相同 | forward 维度相同 | reverse 维度相同 | 两遍维度均相同 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 PE (75 bp) | 89,443 | 100.00% | 99.91% | 99.91% | 45.28% | 45.27% |
| SRR7963242 PE (150 bp) | 2,699,775 | 100.00% | 99.52% | 99.52% | 66.71% | 66.66% |

| 数据 | 总 DP work 比 <=1.10 | forward Lazy-F | reverse Lazy-F | forward 成对 work 上界 | reverse 成对 work 上界 |
| --- | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 PE (75 bp) | 67.20% | 27.75% | 71.96% | 1.938x | 1.754x |
| SRR7963242 PE (150 bp) | 69.57% | 37.71% | 74.77% | 1.907x | 1.806x |

| 数据 | 全部 mate-SW DP 投影 | DP 占 worker cycle | stage2 理想上界 |
| --- | ---: | ---: | ---: |
| ERR1203383 PE (75 bp) | 1.169x | 11.88% | 1.017x |
| SRR7963242 PE (150 bp) | 1.391x | 31.12% | 1.096x |

`stage2` 实测分别为 3.865 s 和 46.686 s；开启细粒度统计会增加少量开销，
这些时间不用于和非 profiling 版本做性能比较。

## 解释

1. **forward 长度不是实现障碍。** 两组数据的 `qlen` 全部相同，至少
   99.5% 的 paired candidates 具有完全相同的 forward
   `qlen/tlen/slen/score_size`。
2. **reverse 不能直接假定同维度。** 第一次 DP 的 `qe` 会决定 reverse
   query 的截断长度；两条 alignment 的 reverse 维度完全相同率只有 45.28%
   和 66.71%。双状态 reverse kernel 必须支持独立边界和 lane mask。
3. **控制流并非完全同步。** 150 bp 中只有 69.68% 的 forward work、
   61.14% 的 reverse work 相差不超过 10%；reverse 的 Lazy-F 占比约 75%，
   是主要分歧来源。
4. **分歧没有吃掉主要收益。** 按真实 vector-step 数取两条 alignment 的
   `max(work0, work1)`，150 bp 的 paired forward/reverse 仍分别接近 1.91x、
   1.81x；折算所有未配对候选后，mate-SW DP 自身约为 1.391x。
5. **整个 stage2 的空间温和。** 150 bp 的 forward+reverse DP 约占
   `WORKER ALIGNMENT` cycle 的 31.1%，按该工作量模型折算，stage2 理想上界
   约 1.096x。75 bp 只有约 1.017x，不值得优先实现。
6. **建议先做 150 bp forward 双状态 kernel。** 它覆盖主要 cycle，状态和
   维度几乎完全兼容，且比同时处理 reverse/Lazy-F 的实现风险小；测到真实
   收益后，再决定是否扩展 reverse。

这里的 work 是实际执行的主循环与 Lazy-F **vector-step 数**。主循环和
Lazy-F 每一步的指令成本并不完全相等，也未计入 query profile 构造、双状态
打包和额外分支，因此 `1.391x/1.096x` 是筛选实现方向的上界，不是承诺值。

## 复现解析

```bash
python3 scripts/analyze_matesw_work_profile.py \
  correctness_results/matesw_work_profile_20260828/ERR1203383_PE.log \
  correctness_results/matesw_work_profile_20260828/SRR7963242_PE.log
```
