# SA batch 与 chain arena 交互实验

## 目的

上一轮子项消融显示 SA batch 和 chain arena 单独开启时均有小幅收益，但
`all` 组合没有表现出线性叠加。本轮只比较：

- `chain_arena`：FP16 KSW + 24 KiB chain arena；
- `sa_chain_arena`：在上述配置上额外开启四路 SA batch。

## 方法

- `cgs_cross + pool`，单节点、单进程、6 CG、384 CPE；
- `ERR1203383 PE`，75 bp reads；
- 每项运行两次，part3 相对极差超过 1% 才补第三次；
- 两项波动均低于 0.3%，未补跑；
- 每项额外生成完整 SAM 并核对 MD5。

## 结果

| 变体 | Stage 2 中位数 | 相对变化 | Part 3 中位数 | 相对变化 | Worker cycles |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain arena | 3.8415 s | 基线 | 3.480 s | 基线 | 7.2330G |
| SA batch + chain arena | 3.8350 s | -0.17% | 3.515 s | +1.01% | 7.3095G |

SA batch 使 chain build cycles 从 `1.4615G` 增至 `1.5730G`，worker cycles
也增加约 `1.06%`。Stage 2 的 `0.17%` 差异小于可解释波动，并且方向与
part3 和 CPE cycles 相反，不能视为收益。

两项完整输出均为：

- `2,000,259` 行；
- `523,048,834` 字节；
- MD5 `c2af4bf0b057d5125ce2a9d770f13741`；
- 状态 `PASS`。

## 结论

最终 BWT/chain 优化仅保留 **chain arena**。SA batch 没有增量收益，连同
上一轮已判定无效的 RID LDM 和 prefetch trim 一并从生产代码删除。

原始逐轮数据见 `analysis/runs.tsv`，汇总见 `analysis/summary.tsv`。
