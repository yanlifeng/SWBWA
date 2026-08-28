# 150 bp same-PE forward 16+16 KSW

## 目的

在 paired-end mate rescue 中，将同一 PE 两个方向、相同 ordinal 的 forward
KSW 候选合并到一次 512-bit SIMD 计算中。低 16 个 FP16 lane 保存方向 0，
高 16 个 FP16 lane 保存方向 1；reverse 起点搜索仍逐候选调用原
`int32_16` 实现。

该路径由 `MATESW_DUAL_FORWARD=1` 控制。只有满足以下条件才启用：

- query 长度均为 150 bp；
- 两个候选的 query 长度和 target 长度分别相等；
- forward 使用 u8 分数；
- 主 KSW backend 为 `int32_16`。

其他长度、不同 DP 维度、i16 分数以及其他 KSW backend 均退回原来的两次
单状态调用。

## 语义处理

- 两个状态各自维护 `gmax`、`te`、次优峰值和提前停止状态。
- Lazy-F 中一个状态先收敛后冻结其 `H1`，另一个状态可继续传播。
- 每个 16-lane 半区独立移位，显式清除 lane 0 和 lane 16 的跨状态边界。
- `qe` 按原 `int32_16` kernel 的物理槽位扫描和 tie 顺序选取，包括 padding
  槽位，以保持现有 SAM 的逐字节兼容性。
- FP16 对本路径使用的 0--255 整数分数可精确表示。

## 正确性

配置：`cgs_cross + pool`、`-t 1 -1 -I 170,80,500,1`。

| 数据 | 长度 | 结果 | 行数 | 字节数 | MD5 |
| --- | ---: | --- | ---: | ---: | --- |
| SRR7963242 PE 子集 | 150 bp | baseline = dual | 5,033,412 | 2,107,655,733 | `a799dc7268f389120ec92820e04b5118` |
| ERR1203383 PE 子集 | 75 bp fallback | 通过既有基准 | 2,000,259 | 523,048,834 | `c2af4bf0b057d5125ce2a9d770f13741` |
| SRR7963242 前 1000 对 | 150 bp | baseline = dual | 2,012 | 854,033 | `bf487e58f3ae0a031074a304fea0a4cb` |

完整 SAM 只在神威端校验，没有拷回本机。
前 1000 对数据另以不带 `-1` 的流水线运行一次，MD5 仍为
`bf487e58f3ae0a031074a304fea0a4cb`。

## 性能

数据为 SRR7963242 150 bp PE 子集。为隔离共享文件系统写入波动，正式复测将
SAM 写到 `/dev/null`；stage2 仍包含完整 alignment 和 SAM record generation。

| 版本 | stage2 run 1 | stage2 run 2 | stage2 均值 | CPE part3 均值 |
| --- | ---: | ---: | ---: | ---: |
| baseline | 46.295 s | 46.292 s | 46.294 s | 45.005 s |
| FP16 16+16 forward | 44.288 s | 44.236 s | 44.262 s | 43.299 s |

- stage2 平均缩短约 **4.39%**；
- CPE alignment/SAM-length pass 平均缩短约 **3.79%**；
- 写普通 SAM 的独立对照为 `46.325 -> 44.245 s`，趋势一致。

逐次原始值见 `summary.tsv`，运行日志见 `benchmark/`、`baseline/`、
`dual_fp16.run.log` 和 `ERR1203383_PE.run.log`。

## 构建

```bash
bash build_cross.sh 8 \
  EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
  HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 \
  CPE_PROFILE=0 KSW_U8_MODE=int32_16 \
  MATESW_DUAL_FORWARD=1 USE_MPI=0
```

使用 `MATESW_DUAL_FORWARD=0` 可构建同源基线。`CPE_PROFILE=0/1` 两种配置均
已通过神威目标编译；profiling 构建仅出现 lwpf3 和既有静态链接警告。
