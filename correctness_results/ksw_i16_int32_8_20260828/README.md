# KSW i16 int32_8 后端

## 目的

神威兼容层此前没有实现八 lane 有符号 i16 KSW 操作。本实验比较：

- `KSW_I16_MODE=scalar_8`：逐 lane 标量参考实现；
- `KSW_I16_MODE=int32_8`：使用低 8 个 32-bit SIMD lane 表示 8 个逻辑
  i16 lane，高 8 个 lane 始终清零。

测试通过 `-A 2` 强制 150 bp read 进入 i16 路径
（`l_seq * a >= 250`）。两种版本均使用 `cgs_cross + pool`、单 MPE 线程、
从核格式化，以及标准的两遍 `build_cross.sh` 构建流程。

## 正确性

完整 `SRR7963242_1.fastq_1` SE 输入在两个后端上产生完全一致的输出：

| 后端 | SAM 记录数 | SAM 字节数 | MD5 |
| --- | ---: | ---: | --- |
| `scalar_8` | 2,514,042 | 979,861,986 | `5e9580e6df968a4c82ad3875fe2a1bda` |
| `int32_8` | 2,514,042 | 979,861,986 | `5e9580e6df968a4c82ad3875fe2a1bda` |

10,000 read 冒烟测试同样逐字节一致：10,048 条记录、3,948,750 字节，
MD5 为 `7fbbaf7d1b534f39b6eb54ff78a2cd7c`。

## 性能

| 后端 | stage2 (s) | CPE part3 (s) |
| --- | ---: | ---: |
| `scalar_8` | 13.280 | 12.650 |
| `int32_8`, run 1 | 12.850 | 12.070 |
| `int32_8`, run 2 | 12.835 | 12.077 |

相对 `scalar_8`，SIMD 两次运行的均值使 stage2 缩短 3.29%，CPE part3
缩短 4.56%。因此默认使用 `int32_8`，同时保留 `scalar_8` 作为正确性参考。

目录内的 PE 和 helper 日志是定位后端问题时留下的探索性诊断，不参与上述
正确性和性能结论。

## 构建

```bash
bash build_cross.sh 8 \
  EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
  HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 \
  CPE_PROFILE=0 KSW_U8_MODE=int32_16 \
  KSW_I16_MODE=int32_8 MATESW_DUAL_FORWARD=1 USE_MPI=0
```
