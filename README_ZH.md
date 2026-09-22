# SWBWA 中文说明

SWBWA 是面向新一代神威平台优化的高精度、高性能短序列比对软件，基于 BWA-MEM 实现。

## 项目特点

- 面向神威异构架构重新设计并行框架，在大共享模式下使用软件预取和访存优化。
- 在双路 48 核 x86 服务器上，相比未优化的单线程版本最高可获得 **330 倍**加速；相比 [bwamem](https://github.com/lh3/bwa) 可获得 **1.2–1.4 倍**加速，结果基本一致。

## 目录结构

- `src/host/`：主核 C/C++ 源代码，包括比对流程、MPI 输入输出和主核辅助模块。
- `src/slave/`：从核 CPE 源代码和从核专用头文件。
- `include/`：公共头文件、配置头文件和生成的 CPE 布局头文件。
- `tools/`：交叉段构建所需的地址和 TLS 信息提取脚本。
- `tests/`：不参与默认构建的 MPI/RMA 和运行时诊断程序。
- `scripts/`：正确性检查、性能测试和结果分析脚本。
- `correctness_results/`：运行日志、正确性结果和测试说明。

## 构建

SWBWA 仅支持新一代神威平台。

### 依赖

- `sw9gcc`（7.1.0 或更高版本）
- `zlib`

### 编译

```bash
git clone https://github.com/RabbitBio/SWBWA.git
cd SWBWA
make -j4
```

### 构建配置

默认构建使用单 CG 执行、CPE FASTQ 格式化、系统 CPE 分配器、动态 MPI 输入和无序单文件 MPI 输出：

```bash
make print-config
```

FASTQ 格式化始终在 CPE 上执行，没有单独的格式化模式开关。

支持的构建变量如下：

| 变量 | 可选值 | 默认值 |
| --- | --- | --- |
| `EXEC_MODE` | `single`、`cgs`、`cgs_cross` | `single` |
| `CPE_ALLOCATOR` | `system`、`pool` | `system` |
| `HOST_MALLOC_WRAPPER` | `0`、`1` | `1` |
| `HOST_MALLOC_STATS` | `0`、`1` | `0` |
| `CPE_KERNEL_OPT` | `0`、`1` | 非 MPI `cgs_cross + pool` 为 `1`，其他为 `0` |
| `CPE_LDM_MODE` | `0` 全关、`1` 分级 malloc 池、`2` 手工优化 | `2`；模式1仅支持非 MPI `cgs_cross + pool` |
| `CPE_DISCARD_DIGEST` | `0`、`1` | 非 MPI `cgs_cross + pool` 且 `OUTPUT_MODE=discard DISCARD_HASH_BYTES=0` 时为 `1`，其他为 `0` |
| `USE_MPI` | `0`、`1` | `1` |
| `MPI_INPUT_MODE` | `static`、`dynamic` | MPI 构建时为 `dynamic` |
| `OUTPUT_MODE` | `split`、`single_unordered`、`discard` | MPI 为 `single_unordered`，非 MPI 为 `split` |
| `MPI_EXACT_READ_INDEX` | `0`、`1` | `1` |

`MPI_INPUT_MODE` 和 `MPI_EXACT_READ_INDEX` 仅在 `USE_MPI=1` 时生效。
非 MPI 也支持 `OUTPUT_MODE=discard`；非 MPI 的 `split` 保持普通单文件输出，
`single_unordered` 则必须启用 MPI。

`CPE_LDM_MODE` 只控制显式私有对象放置，不改变运行时 cache、栈、SIMD 或算法。
模式0关闭新池和历史手工 LDM 申请；模式1关闭旧手工放置，使用每 CPE 32 KiB
分级池接管 malloc/calloc/realloc/free，不适合或放不下的请求回退原交叉段 pool；
模式2关闭新池，保留已验证的手工工作区/阶段复用优化。模式1的 SAM 和公共缓存
query 仍走 heap，不是无条件将所有 malloc 放入 LDM。池元数据计入统一40 KiB预算。
旧的 `CPE_LDM_ALLOC`、`CPE_MANUAL_LDM`、`CPE_LDM_BYTES`、`LDM_SCRATCH_BUDGET`
不再是 Makefile 参数，使用时会明确报错。详细约束见 [LDM 说明](docs/LDM_ALLOCATOR.md)。

```bash
# 将 CPE_LDM_MODE 换成 0/1/2；cgs_cross 必须完整两遍编译。
bash build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 CPE_LDM_MODE=2
```

`MPI_EXACT_READ_INDEX=1` 用于正确性检查。程序会在比对前由 rank 0 扫描完整 FASTQ，建立精确的记录前缀索引。

动态 MPI 输入在 FASTQ 前 90% 使用配置的大块；最后 10% 使用四分之一大小的尾部块，最后两轮 rank 再使用四分之一大小的细尾块。可通过 `SWBWA_MPI_TAIL_PERCENT=0` 关闭尾部细分，或通过 `SWBWA_MPI_FINE_TAIL_WAVES=0` 保留中等尾部而关闭细尾区域。

编译时可用 `MPI_TAIL_PERCENT=0..100` 设置尾部比例默认值，运行环境变量仍可覆盖。
`OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`（默认）对每次提交的完整 SAM blob 计算指纹；
正整数 `DISCARD_HASH_BYTES=N` 只计算前 N 字节，**仅用于测量哈希成本，不能做全输出正确性验证**。
日志包含 `hash_prefix_bytes`，`scripts/check_discard_hash.py` 会拒绝截短哈希。
`EXTRA_CPPFLAGS` 可传递受 `#ifndef` 保护的实验参数；修改参数后必须完整重编译，
cross 构建仍需 `build.sh` 的两遍流程。

当前 CPE 分配器对已包装的 LDM scratch 使用 40 KiB 总预算，分配失败回退 heap，
KSW profile 和 mate-dedup 的单次 LDM 上限各为 16 KiB。这不包含完整栈/静态数据用量。
关闭 `CPE_KERNEL_OPT` 时，context + SMEM 会使独立的 24 KiB chain arena 通常回退 heap。
开启后，chain arena 复用收集阶段结束后不再使用的 16 KiB SMEM scratch；
实际容量按 backing allocation 计算，不把 16 KiB 当成 24 KiB 使用。
`-v 4` 的 pool 诊断输出每批次 LDM 峰值、拒绝次数和结束时未归还字节数。

`OUTPUT_MODE=discard` 是 profiling 模式，不生成 SAM 文件。默认会对提交到输出接口的非空 SAM 数据计算与顺序无关的 sum/XOR 指纹；只有在测量哈希开销时才建议设置 `SWBWA_DISCARD_HASH=0`。

`CPE_DISCARD_DIGEST` 仅在同时满足
`EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`
时默认设为 `1`，其他情况默认 `0`。比较纯核心 stage2 性能或验证最终复制路径时，
显式设置 `CPE_DISCARD_DIGEST=0`：主核在 stage3 对复制后的最终 SAM 缓冲区计算哈希。
符合条件时的默认值 `CPE_DISCARD_DIGEST=1` 面向整体 discard 吞吐：
哈希移至 stage2 part5 的 CPE 复制 worker，最终 SAM 复制由逐 read 摘要元数据替代。
报告时分别列出 part3、part5 和 stage3，不将输出工作迁移当作比对算法加速。

启用路径要求 `EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0`
以及 `OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`。它检查**生成的完整 SAM blob 字节**，仅排除末尾 NUL，
保留逐 read/mate 的 blob 边界及 `calls/bytes/sum/xor` 语义；
**不执行也不验证被跳过的最终 SAM 复制**，日志明确标记
`hash_scope=generated_sam final_sam_copy=0`。正常 SAM 输出模式不变。
`SWBWA_DISCARD_HASH=0` 同时关闭 CPE 哈希，保留 calls/bytes，但不构成 FULL 正确性证据。
构建时传 Make 变量 `CPE_DISCARD_DIGEST`，不要在 `EXTRA_CPPFLAGS` 中重复定义；
修改开关后仍须完整执行两遍构建。

### 单进程 CPE 核心优化

`CPE_KERNEL_OPT=1` 组合启用下列等价实现，默认只针对非 MPI `cgs_cross + pool`：

- SMEM/chain 复用同一段 LDM scratch，不改变 seed 的生命周期和算法顺序。
- `ksw_extend2()` 的 EH 和 QP 合计不超过 4096 字节时优先放入 LDM；
  预算不足或请求过大时保留原来的可复用 heap scratch。
- 五字符 alphabet 的 QP 构建共享 query load，保持原来的行布局。
- 整数 SIMD 合并 gap 下限截断，利用 predicate/XOR 实现选择。
- 主核复用 SAM 指针/长度 scratch，消除主核重复写入终止符。

不改变 lane 数量、打分、lazy-F 终止、tie 比较、任务划分，也不提高 40 KiB LDM 总预算。
用 `CPE_KERNEL_OPT=0` 构建核心对照组；输出 hash 加速和按批次计时独立于此开关。

纯核心示例：单进程六 CG 的 stage2 测量与 FULL 最终复制校验，显式设置 `CPE_DISCARD_DIGEST=0`：

```bash
./build.sh 6 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 \
    CPE_KERNEL_OPT=1 OUTPUT_MODE=discard DISCARD_HASH_BYTES=0 CPE_DISCARD_DIGEST=0 CPE_PROFILE=0
SWBWA_DISCARD_HASH=1 bsub -I -b -q q_share -n 1 -cgsp 64 -mpecg 6 \
    -share_size 2000 -xmalloc -cross_size 42000 -cache_size 128 -priv_size 16 \
    ./SWBWA mem -v 4 -t 1 -1 -I 170,80,500,1 -o ignored.sam \
    ref.fa read1.fq read2.fq
```

`-1` 将三个 stage 串行执行。比较 stage2 和其内部 part3，不把输入波动或 hash 加速
算作 CPE 比对加速。正确性检查必须同时匹配 `calls/bytes/sum/xor`，
并确认 `enabled=1 hash_prefix_bytes=0`。指纹是回归检查，不是无碰撞证明或 SAM MD5。

采用默认 discard 吞吐配置时，用下列命令重新构建，再使用上面的运行命令。
此配置下 `CPE_KERNEL_OPT` 和 `CPE_DISCARD_DIGEST` 均默认 `1`；FULL 校验覆盖生成的 SAM
（`hash_scope=generated_sam`），不验证被跳过的最终复制。

```bash
./build.sh 6 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 \
    OUTPUT_MODE=discard DISCARD_HASH_BYTES=0 CPE_PROFILE=0
```

本机可运行 `bash tests/run_host_checks.sh` 和 `bash tests/run_cpe_kernel_checks.sh`。
后者使用编译器向量适配层和 sanitizer 比较实际 CPE 函数，不替代神威端完整输出验证。

示例：

```bash
make clean
make -j4 EXEC_MODE=cgs CPE_ALLOCATOR=system
```

切换构建模式前请先执行 `make clean`，因为 Make 不会根据所有编译宏的变化自动重新编译旧对象。

交叉段执行需要重新生成 CPE 布局和重定位数据：

```bash
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=system
```

`build.sh` 的第一个参数是并行编译数，后面直接传入 Make 变量。它会自动完成占位布局构建、ELF 段地址提取、布局头文件生成、第二遍构建以及交叉段重定位信息生成。

常用配置示例：

```bash
./build.sh 8 EXEC_MODE=single CPE_ALLOCATOR=system USE_MPI=0
./build.sh 8 EXEC_MODE=cgs CPE_ALLOCATOR=pool USE_MPI=0
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=1 \
    MPI_INPUT_MODE=dynamic OUTPUT_MODE=single_unordered
```

## 基本用法

### 建立参考序列索引

SWBWA 兼容由 [bwamem](https://github.com/lh3/bwa) 生成的索引文件。

```bash
./SWBWA index ref.fa
```

### 比对 reads

单端测序：

```bash
./SWBWA mem ref.fa reads.fq -o aln.sam
```

双端测序：

```bash
./SWBWA mem ref.fa read1.fq read2.fq -o aln.sam
```

MPI 作业中每个 rank 独立读取自己的 FASTQ 区间，不由 rank 0 把整个输入读入内存：

```bash
bsub -I -b -q q_share -N 1 -np 6 -cgsp 64 \
  -share_size 12000 -cache_size 128 -priv_size 16 \
  ./SWBWA mem -t 1 -1 -K 5000000 -o out.sam \
  ref.fa read1.fq read2.fq
```

`OUTPUT_MODE=split` 为每个 rank 生成独立输出；`single_unordered` 使用 MPI RMA 原子申请单文件区间，但不保证记录顺序；`discard` 不写 SAM，适合测量 stage2 性能。

统一运行入口只负责提交作业，不负责重新编译：

```bash
./run.sh single -- ./SWBWA mem -t 1 -o out.sam ref.fa reads.fq
./run.sh cgs_cross -- ./SWBWA mem -t 1 -o out.sam ref.fa reads.fq
./run.sh mpi --nodes 1 --ranks 6 -- \
  ./SWBWA mem -t 1 -K 5000000 -o out.sam ref.fa reads.fq
```

测试脚本和结果目录约定见 [`scripts/README.md`](scripts/README.md) 与 [`correctness_results/README.md`](correctness_results/README.md)。

远程登录地址、账号、凭据和私人运维说明请在仓库之外保存，不要提交到 Git。实验脚本所需的远程连接参数应通过本地环境变量显式传入。

## 获取帮助

```bash
./SWBWA mem
```

该命令会显示完整的参数说明。
