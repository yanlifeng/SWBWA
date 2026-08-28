# Mate-SW 跨调用双 alignment 可行性统计

## 目的

现有 `mem_matesw()` 每次最多产生一个实际 KSW alignment，无法直接利用
32 个物理 lane 同时执行两个 16-lane DP。本轮把统计提升到 `mem_sam_pe()`：

- `same-PE`：配对同一 PE read 两个 mate-rescue 方向中的候选。
- `adjacent-PE`：配对同一 CPE 连续处理的两次 `mem_sam_pe()` 的候选。

这里只统计候选数量，是实现前的快速上界评估；没有检查两个 alignment 的
`qlen/tlen` 是否接近，也没有计入 SIMD 分支分歧和打包开销。

## 配置

- `EXEC_MODE=cgs_cross`
- `CPE_ALLOCATOR=pool`
- `CPE_PROFILE=1`, sampled CG 5
- `KSW_U8_MODE=int32_16`
- `-t 1 -1 -v 4 -I 170,80,500,1`
- 输出写入 `/dev/null`
- 数据：`ERR1203383` PE 75 bp 子集、`SRR7963242` PE 150 bp 子集

## 结果

| 数据 | PE 调用数 | KSW 候选数 | 两方向均活跃 | same-PE 候选覆盖 | same-PE 理想加速上界 | adjacent-PE 候选覆盖 | adjacent-PE 理想加速上界 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 PE (75 bp) | 1,000,000 | 680,376 | 1.79% | 26.29% | 1.151x | 9.87% | 1.052x |
| SRR7963242 PE (150 bp) | 2,500,000 | 9,372,375 | 15.71% | 57.61% | 1.405x | 6.71% | 1.035x |

理想加速上界按 `candidates / (candidates - pairs)` 计算，假定被配对的两个
alignment 工作量完全一致。因此它只能用于筛选方向，不能当作预计实测加速。

## 结论

1. 单个 `mem_matesw()` 内依旧没有双 alignment 机会。
2. 150 bp 数据中，同一 PE 的两个方向覆盖 57.61% 候选，值得继续做
   `qlen/tlen` 兼容性和工作量分布统计；理想上界约 1.405x。
3. 75 bp 的 same-PE 上界只有 1.151x，收益空间较小。
4. 简单固定配对相邻 PE 仅覆盖 6.71%--9.87%，不值得优先实现。
5. 下一步若继续，应只细化 same-PE 路径：统计候选对的长度比和 DP work
   比，再决定是否承担双状态 KSW 的实现复杂度。

## 复现解析

```bash
python3 scripts/analyze_matesw_cross_call.py \
  correctness_results/matesw_cross_call_profile_20260828/ERR1203383_PE.log \
  correctness_results/matesw_cross_call_profile_20260828/SRR7963242_PE.log
```
