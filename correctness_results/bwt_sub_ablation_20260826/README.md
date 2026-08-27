# BWT 子项消融实验

## 实验目的

将原先由 `SWBWA_ENABLE_BWT_CHAIN_OPT` 一次性控制的四项修改拆开，判断
各自对 CPE alignment 主路径的真实贡献：

1. BWT SA 四路批处理；
2. contig/RID offset 常驻 LDM；
3. 24 KiB chain 临时分配 arena；
4. 裁掉前向 BWT 搜索中的一组预取。

所有变体均保留已经通过正确性验证的 FP16 KSW 优化。`ksw_only` 是本次
基线，只关闭上述四个 BWT/chain 子项。

## 实验设计

- 平台配置：`cgs_cross + pool`，单节点、单进程、6 CG、384 CPE；
- 数据：`ERR1203383 PE` 子集，75 bp reads；
- 原因：该数据在此前 profiling 中的 BWT/chain 热点占比最高，适合快速
  放大子项差异；
- 性能输出：`/dev/null`，排除实际 SAM 文件写入；
- 每个变体先跑 2 次；若 `part3` 相对极差超过 1%，才补第 3 次；
- 本轮全部变体的相对极差均小于 0.3%，没有触发补跑；
- 每个变体额外生成一次完整 SAM，核对行数、字节数和 MD5。

## 性能结果

以下为两次运行的中位数；变化均相对 `ksw_only`。

| 变体 | Stage 2 (s) | Stage 2 变化 | Part 3 (s) | Part 3 变化 | Part 3 极差 |
| --- | ---: | ---: | ---: | ---: | ---: |
| KSW only | 3.9005 | 基线 | 3.530 | 基线 | 0.000% |
| + SA batch | 3.8040 | -2.47% | 3.484 | -1.30% | 0.230% |
| + RID offsets in LDM | 4.0300 | +3.32% | 3.675 | +4.11% | 0.272% |
| + chain arena in LDM | 3.8465 | -1.38% | 3.475 | -1.56% | 0.288% |
| + prefetch trim | 3.9000 | -0.01% | 3.530 | 0.00% | 0.000% |
| + all BWT/chain changes | 3.7915 | -2.80% | 3.475 | -1.56% | 0.288% |

### BWT/chain 热点

单位为 sampled CPE 的平均十亿 cycles。

| 变体 | Worker | Chain build | SA | RID | Insert | Finalize |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KSW only | 7.3375 | 1.5675 | 0.9755 | 0.1060 | 0.2610 | 0.1455 |
| + SA batch | 7.2450 | 1.6085 | 1.0105 | 0.0990 | 0.2475 | 0.1345 |
| + RID offsets in LDM | 7.6705 | 1.7800 | 0.9395 | 0.1060 | 0.4300 | 0.2070 |
| + chain arena in LDM | 7.2265 | 1.4610 | 0.9575 | 0.1050 | 0.1195 | 0.2020 |
| + prefetch trim | 7.3435 | 1.5675 | 0.9740 | 0.1060 | 0.2625 | 0.1450 |
| + all | 7.2245 | 1.5440 | 1.0180 | 0.0970 | 0.1240 | 0.1635 |

## 正确性

六个变体全部得到相同结果：

- 行数：`2,000,259`
- 字节数：`523,048,834`
- MD5：`c2af4bf0b057d5125ce2a9d770f13741`
- 状态：全部 `PASS`

通过后生成的 SAM 已删除，只保留日志和统计表。

## 结论

1. **chain arena 值得保留。** 它将 chain insert cycles 从 `0.2610G`
   降到 `0.1195G`，使 part3 稳定下降约 `1.56%`。
2. **SA batch 有小幅稳定收益。** part3 下降约 `1.30%`，完整输出正确。
   不过 profile 中 SA 区域自身 cycles 没有下降，收益来自更广的 worker
   执行路径，因此暂不应继续扩大 batch 大小。
3. **RID LDM 单独开启有明显负收益。** RID 区域自身没有变快，反而伴随
   insert/finalize 和 worker 总 cycles 上升；当前实现不应作为默认优化。
4. **prefetch trim 可删除。** 两轮 wall time 和热点数据均与基线等价，
   没有可测收益。
5. **四项收益不线性叠加。** 后续交互实验表明，在 chain arena 上增加
   SA batch 会使 part3 变慢约 `1.01%`。最终仅保留 chain arena；SA
   batch、RID LDM 和 prefetch trim 均从生产代码删除。详见
   `../bwt_sa_arena_interaction_20260826/README.md`。

原始逐轮数据见 `analysis/runs.tsv`，汇总见 `analysis/summary.tsv`，完整
正确性结果见 `correctness/results.tsv`。
