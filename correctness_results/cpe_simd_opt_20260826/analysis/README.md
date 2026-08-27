# CPE SIMD、BWT 与 chain 优化实验

> **后续结论：** 本文记录的是第一轮组合实验。当时列出的 SA batch、RID
> offsets LDM 和 BWT prefetch 修改，已在后续独立消融与交互实验中判定为
> 无稳定增量收益并从生产代码删除。最终只保留 FP16 16-lane backend、
> 精确 KSW 分配和 24 KiB chain arena；详见
> `../../bwt_sub_ablation_20260826/README.md` 与
> `../../bwt_sa_arena_interaction_20260826/README.md`。

## 实验范围

- 执行模式：`cgs_cross + pool`，单节点、单进程、6 个 CG、384 个 CPE。
- 输入：`ERR1203383` PE 子集（75 bp）与 `SRR7963242` PE 子集（150 bp）。
- 正确性：保留方案均以完整无 header SAM 的行数、字节数和 MD5 校验；性能结果来自 `-1` 模式。
- 从核统计：LWPF3 的 cycle、instruction、D-cache 和 instruction-buffer stall；子 region 的计时包含在父 region 内，不能相加到父 region 上。

结果可由仓库根目录执行以下命令重新生成：

```bash
python3 scripts/analyze_cpe_simd_opt.py
```

## 最终正确性

最终 `CPE_PROFILE=0` 版本生成完整、无 header 的 SAM，并与既有基准逐字节一致：

| 数据 | SAM 行数 | SAM 字节数 | MD5 |
| --- | ---: | ---: | --- |
| ERR1203383 PE（75 bp） | 2,000,259 | 523,048,834 | `c2af4bf0b057d5125ce2a9d770f13741` |
| SRR7963242 PE（150 bp） | 5,033,412 | 2,107,655,733 | `a799dc7268f389120ec92820e04b5118` |

## 第一轮组合中保留、后续继续消融的修改

1. **KSW u8 使用 FP16 算术，保持 16 条逻辑 lane。**
   神威整数 SIMD 是 16 x int32，缺少原生 int8/int16 饱和算术；FP16 则有 32 条物理 lane，并且 0 到 255 的整数可以精确表示。实现只使用低 16 条 FP16 lane，高 16 条保持为 0，因此 stripe 宽度、lazy-F 次数、最大值扫描顺序和 tie 语义都与原实现一致。
2. **BWT SA 以 4 个查询为一批做 inverse-Psi。**
   每轮先预取 4 个相互独立的 Occ block，再按原 seed 顺序消费 SA 结果。batch 4、8、16 均通过正确性校验，batch 4 的 Stage2 最低且寄存器压力更小。
3. **contig offset 表复制到 LDM。**
   当 `n_seqs + 1` 个 offset 不超过 8 KiB 时，从核用本地二分完成 `bns_intv2rid`；否则回退原函数。
4. **每个 CPE 使用 24 KiB 的单-read chain seed arena。**
   chain seed 初始分配和增长优先走 LDM bump allocator，容量不足时回退现有堆分配。arena 每条 read 重置，堆 fallback 仍按原生命周期释放。
5. **去掉 forward BWT 热循环中不会被下一步立即使用的 x0 预取。**
   保留 x1 的即时预取；距离 8 的超前预取实测变慢，未保留。

显式 LDM 峰值预算约为：SMEM 临时向量 16 KiB + chain arena 24 KiB + contig offsets 至多 8 KiB + KSW 临时结构通常约 6 KiB，明显低于 128 KiB 手动 LDM 预算；另 128 KiB 由 `-cache_size` 用作自动 cache。KSW 分配公式原先残留了 int32 模拟时期的 4 倍系数，典型 150 bp query 会申请约 23 KiB；改为按 `sizeof(__m128i)` 精确计算后约为 6 KiB。

## KSW 结论

| 数据 | 实现 | Stage2 | CPE part3 | 正确性 |
| --- | ---: | ---: | ---: | --- |
| ERR 75 bp | int32 baseline | 3.759 s | 3.400 s | 通过 |
| ERR 75 bp | FP16 32 logical lanes | 3.715 s | 3.360 s | **失败** |
| ERR 75 bp | FP16 16 logical lanes | 3.925 s | 3.570 s | 通过 |
| SRR 150 bp | FP16 32 logical lanes | 43.201 s | 42.221 s | **失败** |
| SRR 150 bp | FP16 16 logical lanes | 45.281 s | 44.321 s | 通过 |

32-lane 版本虽然更快，但改变了 stripe 宽度和 KSW 的启发式遍历顺序。ERR 输出中有 259 条 SAM 发生 MAPQ、MQ、XS 或 tie 选择变化，因此已拒绝。16-logical-lane 版本将 mate KSW cycle 在 ERR 上从约 932M 降至 816M，在 SRR 上从约 32G 降至 29G，同时保持逐字节输出一致。

150 bp 的 mate rescue 满足 `l_ms * a < 250`，实际走 `ksw_u8`；本轮数据并不进入 `ksw_i16`。当前 CPE 的 `ksw_i16` SSE 兼容函数仍是历史遗留的 `TODO` 退出路径，因此本轮结论只覆盖已验证的 u8 路径。要利用高 16 条 FP16 lane 而不改变语义，需要把两个相互独立的 alignment 打包进同一向量，并分别维护 target 长度、lazy-F、最大值和提前退出状态，改动面和风险都明显更大，本轮未实现。

## 75 bp BWT 与 chain 结论

最终 profile 中，`MEM_CHAIN_BUILD` 占 worker cycle 的 22.5%，其中：

- `CHAIN_BUILD_SA`：约 66.0% 的 build cycle，是首要热点；
- `CHAIN_BUILD_INSERT`：约 8.4%；
- `CHAIN_BUILD_FINALIZE`：约 9.3%；
- `CHAIN_BUILD_RID`：约 7.1%。

SA batch 测试的 Stage2 分别为 batch 4：3.822 s、batch 8：3.843 s、batch 16：3.843 s。contig offset 本地化使 RID cycle 约下降 34%。组合 FP16x16、SA batch 4、RID LDM 和 chain arena 后，ERR Stage2 为 3.783 s；修正 KSW 分配后复测为 3.786 s，两者可视为持平。相对详细插桩基线 4.070 s 提升约 7.0%。SRR 150 bp 最终复测为 45.355 s，和修正前 45.407 s 持平。

距离 8 的 BWT lookahead prefetch 将 ERR Stage2 从约 3.785 s 增加到 3.869 s，说明额外请求造成的 cache 污染大于隐藏的延迟，因此已回退。

## Mate dedup 结论

候选数组超过 16 个元素的调用比例：

- ERR 75 bp：约 68.8%；
- SRR 150 bp：约 84.4%。

尤其在 150 bp 数据中，`65+` 单独占 50.24%。只覆盖很小数组的专用排序收益面有限，而当前两个排序还参与同分元素的先后语义；本轮不实现专用 sorter。

## 文件

- `run_summary.tsv`：各实验 Stage2/part3 和正确性状态。
- `stage2_comparison.svg`：实验阶梯对比，红色是因输出变化而拒绝的方案。
- `hotspot_hierarchy.tsv` / `.svg`：75 bp 与 150 bp 的父子 kernel 热点。
- `mate_dedup_distribution.tsv` / `.svg`：dedup 候选数组规模分布。
- 上级目录各子文件夹：对应实验的原始轻量日志。
