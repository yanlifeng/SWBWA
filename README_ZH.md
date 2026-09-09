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
| `USE_MPI` | `0`、`1` | `1` |
| `MPI_INPUT_MODE` | `static`、`dynamic` | MPI 构建时为 `dynamic` |
| `OUTPUT_MODE` | `split`、`single_unordered`、`discard` | `single_unordered` |
| `MPI_EXACT_READ_INDEX` | `0`、`1` | `1` |

`MPI_INPUT_MODE`、`OUTPUT_MODE` 和 `MPI_EXACT_READ_INDEX` 仅在 `USE_MPI=1` 时生效。

`MPI_EXACT_READ_INDEX=1` 用于正确性检查。程序会在比对前由 rank 0 扫描完整 FASTQ，建立精确的记录前缀索引。

动态 MPI 输入在 FASTQ 前 90% 使用配置的大块；最后 10% 使用四分之一大小的尾部块，最后两轮 rank 再使用四分之一大小的细尾块。可通过 `SWBWA_MPI_TAIL_PERCENT=0` 关闭尾部细分，或通过 `SWBWA_MPI_FINE_TAIL_WAVES=0` 保留中等尾部而关闭细尾区域。

`OUTPUT_MODE=discard` 是 profiling 模式，不生成 SAM 文件。默认会对提交到输出接口的非空 SAM 数据计算与顺序无关的 sum/XOR 指纹；只有在测量哈希开销时才建议设置 `SWBWA_DISCARD_HASH=0`。

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

## 获取帮助

```bash
./SWBWA mem
```

该命令会显示完整的参数说明。
