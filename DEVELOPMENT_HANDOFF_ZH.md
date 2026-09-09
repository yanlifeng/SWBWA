# SWBWA 开发与研究交接文档

本文档用于在新的机器、账号或开发者接手 SWBWA 时，快速恢复项目背景、构建方式、实验方法和优化脉络。内容以当前仓库 `master` 分支的代码为准，不能替代源代码本身；当文档与代码不一致时，以代码和 `Makefile` 为准。

## 1. 项目定位

SWBWA 是基于 BWA-MEM 的 Sunway 异构平台实现。程序由主核（MPE）和从核（CPE）共同完成：

- 主核负责命令行、索引加载、FASTQ 输入、流水线组织、MPI 输入输出和最终资源回收。
- CPE 负责 FASTQ 格式化、SMEM/BWT seeding、chain build、chain extension、mate rescue 和 SAM 记录生成。
- MPI 版本把输入 FASTQ 划分到多个 rank，并在 rank 层面做静态或动态调度。
- CPE 层面仍然在每个 rank 内部以 read 或 PE pair 为工作单位动态领取任务。
- `cgs_cross` 版本把部分 CPE 代码和私有运行时数据复制到交叉段，以缓解默认大共享模式下单一 CG 主存和代码访问拥堵的问题。

当前工程的主要研究问题可以概括为：

1. 如何把 BWA-MEM 的不规则 seeding、chaining 和 SW 工作映射到 Sunway CPE。
2. 如何利用 LDM、预取、局部内存池和交叉段，降低主存访问与运行时开销。
3. 如何在 MPI rank 之间平衡 FASTQ 计算负载，同时维持足够大的批次以摊薄 RMA、格式化和启动开销。
4. 如何在不改变 SAM 语义的前提下，探索 KSW SIMD 和 PE mate rescue 的批处理优化。

## 2. 当前仓库和 Git 状态

当前推荐使用 `master` 分支：

    0a08920 docs: clarify repository result layout

当前重要的组织提交和性能提交包括：

| 提交 | 含义 |
| --- | --- |
| `4258ae4` | 整理结果目录，统一 `build.sh`/`run.sh` 入口 |
| `029631b` | 整理主核源码和独立测试程序 |
| `2ac0f20` | 将主核、从核、公共头文件放入规范目录 |
| `0a08920` | 进一步说明结果目录布局 |
| `7c3ebe6` | CPE alignment 热路径的 LDM 优化 |
| `f2134d9` | FP16 KSW 实验后端和 chain arena |
| `d6e1582` | 显式选择 CPE KSW u8 后端 |
| `953b2e3` | PE mate rescue 批处理和 i16 SIMD KSW |
| `7d536cc` | 复用标量 `ksw_extend2()` 的临时 scratch buffer |
| `61221c2` | MPI 动态 FASTQ 调度 |
| `3dc7e31` | MPI 输入分块和输出框架 |
| `4fb1ff5` | CPE 执行期间让主核进入 MPI progress |
| `30d8722` | MPI 动态调度尾部细分策略 |
| `3129ab3` | MPI dynamic I/O 和 discard profiling 调整 |
| `fc658b8` | 大数据正确性测试与 CPE profiling |

`perf/cpe-kernel-profile-opt` 是性能实验分支，包含部分与当前 `master` 相近但历史基线不同的提交。日常维护、构建和新实验以 `master` 为准；需要比较历史优化时再切换该分支。

当前仓库已经完成以下整理：

    src/host/       主核 C/C++ 源码
    src/slave/      从核 CPE 源码和从核专用头文件
    include/        公共头文件、配置和生成的 CPE 布局头文件
    tests/          不参与默认构建的 MPI/RMA/进度诊断程序
    tools/          cgs_cross 两阶段构建所需的地址/TLS 工具
    scripts/        正确性、性能、profiling 和结果分析脚本
    correctness_results/
                    运行日志、MD5、TSV、图表和实验说明

早期对话中曾有 `AGENTS.md` 和若干 `.codex_*` 临时文件；当前整理后的可达历史和工作树不再依赖这些文件。不要把旧工作区的临时备份重新拷贝进仓库。

## 3. 五分钟恢复环境

### 3.1 克隆和查看配置

    git clone <repository-url> SWBWA
    cd SWBWA
    git checkout master
    make print-config

需要注意：`include/swbwa_config.h` 提供的是源码层默认值，而当前 `Makefile` 还会通过 `?=` 给出构建入口默认值。因此直接执行 `make` 时，实际配置以 `make print-config` 的输出为准。当前 Makefile 的默认组合是：

    EXEC_MODE=single
    CPE_ALLOCATOR=system
    HOST_MALLOC_WRAPPER=1
    HOST_MALLOC_STATS=0
    CPE_PROFILE=0
    KSW_U8_MODE=int32_16
    KSW_I16_MODE=int32_8
    MATESW_DUAL_FORWARD=1
    USE_MPI=1
    MPI_INPUT_MODE=dynamic
    OUTPUT_MODE=single_unordered
    MPI_EXACT_READ_INDEX=1

这意味着“默认构建”是单 CG 执行模式的 MPI-enabled 可执行文件；`USE_MPI=0` 才是明确的非 MPI 构建。MPI 的 `MPI_INPUT_MODE` 和 `OUTPUT_MODE` 在非 MPI 构建中没有意义，也不会被加入源码编译宏。

### 3.2 常用构建命令

切换构建模式时先清理，因为 Make 不会自动根据所有编译宏变化重新编译旧对象：

    make clean

非 MPI、单 CG、系统 CPE allocator：

    make -j8 EXEC_MODE=single CPE_ALLOCATOR=system USE_MPI=0

六 CG、不使用交叉段、系统 allocator：

    make -j8 EXEC_MODE=cgs CPE_ALLOCATOR=system USE_MPI=0

六 CG、交叉段、CPE pool：

    ./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0

MPI dynamic + split：

    make clean
    make -j8 EXEC_MODE=single CPE_ALLOCATOR=system USE_MPI=1 \
        MPI_INPUT_MODE=dynamic OUTPUT_MODE=split MPI_EXACT_READ_INDEX=1

MPI dynamic + 单文件无序输出：

    make clean
    make -j8 EXEC_MODE=single CPE_ALLOCATOR=system USE_MPI=1 \
        MPI_INPUT_MODE=dynamic OUTPUT_MODE=single_unordered \
        MPI_EXACT_READ_INDEX=1

MPI enabled 时 Makefile 会把 `mpicxx -show -mhybrid` 展开的链接命令中的
`single_static` 替换为 `multi_static`，因为当前程序申请
`MPI_THREAD_MULTIPLE`，并且 CPE 执行期间有独立的主核 progress thread 进入 MPI。
不要修改系统 MPI 安装目录；这是本项目链接阶段的选择。

### 3.3 `cgs_cross` 的两阶段构建

`cgs_cross` 不是只加一个编译宏。它依赖最终 ELF 中的 CPE text/data 地址，并需要根据这些地址生成 `include/swbwa_cpe_layout.h`，再做第二遍构建；随后还要生成重定位信息。

`build.sh` 的入口形式是：

    ./build.sh [parallel_jobs] [Make variable assignments ...]

例如：

    ./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0

当前脚本会先写占位布局头、构建、调用 `swreadelf` 提取布局、重写布局头、再次构建，并调用 `tools/xlink.py` 与 `tools/get_tls.py` 生成交叉段所需信息。 `include/swbwa_cpe_layout.h` 是生成文件，不要手工填写地址。修改影响 CPE text/data、TLS 或重定位的源码后必须重新跑该流程。

运行交叉段二进制时，作业资源通常还需要对应的 `-xmalloc` 和 `-cross_size` 参数。 `run.sh cgs_cross` 和正确性脚本已经提供了这套资源模板，但直接手写 `bsub` 时不要漏掉。

## 4. 程序架构和数据流

### 4.1 主核流水线

主入口链路：

    src/host/main.c
      -> main_mem()                         [src/host/fastmap.c]
      -> 索引/输入/MPI/输出初始化
      -> kt_pipeline_queue(..., process, ..., 3)
           step 0: 读取 FASTQ 批次
           step 1: CPE 格式化、比对、SAM 长度/记录生成
           step 2: 输出 SAM、释放批次
      -> swbwa_output_close()
      -> MPI scheduler/progress/profile/timing 收尾

`src/host/fastmap.c` 中的 `process()` 是理解整个主核执行顺序的第一入口。当前三阶段的含义是：

- **stage1**：申请输入 buffer，按范围 `fseeko`/`fread`，必要时把 FASTQ 末端对齐到完整 record。
- **stage2**：调用 `mem_process_seqs_merge2()`，完成 CPE FASTQ 格式化、CPE alignment、SAM 长度计算和 SAM 记录生成。
- **stage3**：把每条 read 的 SAM 字符串交给统一输出接口，并释放本批次的 `seqs`/元数据。

`-1` 不是“是否使用 CPE”的开关。在 `fastmap.c` 中它被解析为 `no_mt_io=1`，从而选择 `kt_pipeline_single(1, ...)`；不带 `-1` 则选择 `kt_pipeline_queue(3, ...)`。因此 correctness 结果中的 `has1/no1` 实际是在比较两种主核输入/流水线路径，必须把它当作独立实验变量，不能简单宣称其中一个永远更快。

### 4.2 CPE alignment 链路

核心函数在 `src/slave/bwamem.c`：

    mem_align1_core_impl()
      -> mem_chain()
           -> mem_collect_intv()
           -> bwt_smem1() / BWT Occ / SA
           -> chain build
      -> mem_chain_flt()
      -> mem_flt_chained_seeds()
           -> mem_seed_sw()
                -> ksw_align2()
      -> mem_chain2aln()
           -> ksw_extend2()                 [主 chain extension]
      -> mem_sort_dedup_patch()

PE 额外经过 `src/slave/bwamem_pair.c`：

    mem_sam_pe()
      -> mate rescue
           -> ksw_align2_matesw()
                -> ksw_u8() 或 ksw_i16()
      -> mem_pair()
      -> mem_reg2sam()

一个 read 或 pair 在 chaining 后可能有多个候选 `mem_alnreg_t`。候选会经历过滤、extension 和去重；PE mate rescue 并不是“重复输出已经完成的 alignment”，而是在一个 mate 已经有可靠位置、另一个 mate 位置不完整或不满足配对约束时，根据 insert-size/orientation 约束搜索可能的 mate 区域，再用短 SW 判断候选。这个过程因此可能触发额外的 KSW 工作。

### 4.3 两层动态调度

当前系统有两个完全不同的调度层：

1. **MPI rank 层**：以 FASTQ byte range/chunk 为单位领取任务。dynamic 模式通过 MPI RMA ticket window 领取 chunk。
2. **CPE 层**：一个 CPE kernel 收到一个较大的输入批次后，再以 read（SE）或 read pair（PE）为单位动态领取任务。

因此 CPE 层可以把一个 chunk 内部的 read 级不规则性摊平，但 MPI 层最后仍可能被某个慢 chunk 的尾部工作限制。 `-K` 太小会增加 MPI/RMA/格式化和 kernel 启动次数；`-K` 太大则会减少 rank 级可调度单元，放大最后慢 chunk 的影响。

## 5. 配置开关完整说明

### 5.1 Makefile 对外接口

| Make 变量 | 可选值 | 当前默认 | 作用与注意事项 |
| --- | --- | --- | --- |
| `EXEC_MODE` | `single`, `cgs`, `cgs_cross` | `single` | 分别是 1 CG/64 CPE、6 CG/384 CPE、6 CG + 交叉段。 |
| `CPE_ALLOCATOR` | `system`, `pool` | `system` | CPE 临时分配走系统 allocator 或 segment-tree pool。 |
| `HOST_MALLOC_WRAPPER` | `0`, `1` | `1` | 主核 malloc wrapper。主要为了可选统计和 calloc 兼容。 |
| `HOST_MALLOC_STATS` | `0`, `1` | `0` | 主核 malloc/free/calloc/realloc 统计；需要 wrapper。 |
| `CPE_PROFILE` | `0`, `1` | `0` | 打开 LWPF3 CPE profiling。会有额外 profiling 开销。 |
| `CPE_PROFILE_CG` | `0..5` | single 为 0，否则通常为 5 | 只采样一个 CG，降低 profiling 成本。 |
| `KSW_U8_MODE` | `int32_16`, `float16_16`, `float16_32` | `int32_16` | mate SW 的 u8 路径后端。默认是结果最稳的 int32-16。 |
| `KSW_I16_MODE` | `scalar_8`, `int32_8` | `int32_8` | i16 路径的标量或 8-lane int32 SIMD。 |
| `MATESW_DUAL_FORWARD` | `0`, `1` | `1` | 150 bp、same-PE、forward 条件满足时尝试两个 16-lane KSW。 |
| `USE_MPI` | `0`, `1` | `1` | 是否编译 MPI 输入/输出和 RMA 路径。 |
| `MPI_INPUT_MODE` | `static`, `dynamic` | dynamic | 仅 MPI 构建有效。 |
| `OUTPUT_MODE` | `split`, `single_unordered`, `discard` | single_unordered | 仅 MPI 构建有效。 |
| `MPI_EXACT_READ_INDEX` | `0`, `1` | `1` | 是否预先建立精确 record index，正确性测试建议为 1。 |

源码层的统一宏都在 `include/swbwa_config.h`。不要在不同 `.c` 文件中重新定义模式含义；新增开关前先判断它是否确实需要独立做消融。如果只是内部固定调参，保留为配置头中的内部常量，不要扩散到 Makefile。

### 5.2 CPE 内部常量

当前配置头中还有以下内部常量：

- `SWBWA_ENABLE_DYNAMIC_SCHEDULING=1`：CPE worker 内部动态领取 read/pair。
- `SWBWA_ENABLE_CPE_PREFETCH=1`：BWT/Occ 相关路径使用软件预取。
- `SWBWA_READS_PER_DYNAMIC_TASK=1`：CPE 动态任务粒度为一个 read 或 pair。
- `SWBWA_MAX_TASKS_PER_CPE=(50<<10)`：保存第二遍 SAM 生成所需的 task 列表上限。
- `SWBWA_DEFAULT_FASTQ_BYTES_PER_CG=64*1024*300`：约 19.2 MiB/CG 的默认 FASTQ byte 目标。
- `SWBWA_CPE_POOL_BYTES_PER_CPE=24 MiB`：CPE pool 每个 CPE 的默认空间。
- `SWBWA_CPE_FORMAT_BUFFER_BYTES=512 MiB`：CPE FASTQ 格式化 buffer 上限。
- `SWBWA_OUTPUT_BUFFER_BYTES=64 MiB`：MPI 输出批量 buffer 大小，不是输入 chunk 大小。
- `SWBWA_PIPELINE_QUEUE_CAPACITY=4`、`SWBWA_PIPELINE_BUFFER_COUNT=5`：主核流水线队列容量。

动态 MPI 的 chunk 总字节数大致是：

    chunk_bytes = SWBWA_DEFAULT_FASTQ_BYTES_PER_CG * SWBWA_CG_COUNT

所以同一个 `-K` 在 single CG 和 6-CG 构建中不是同一个总 chunk 大小。程序帮助中的 `-K` 说明是“每个 CG、每个输入文件的目标原始 FASTQ bytes”。这也是早期实验中“大共享模式 CPE 变多但 chunk 没同步变大”导致理解混乱的根源之一。

## 6. FASTQ 输入和 `n_processed`

### 6.1 MPI 是否由 rank 0 全量读取

正常 MPI 运行不是 rank 0 把整个 FASTQ 读入内存。每个 rank 独立打开自己的输入文件并从自己的 byte range 读取：

    rank 0: open/read [start0, end0)
    rank 1: open/read [start1, end1)
    ...

边界通过 `find_fastq_boundary()`/相关 range 函数调整到完整 FASTQ record，不允许把 FASTQ 四行 record 从中间截断。PE 假设两个输入文件大小相同，并对 read1/read2 使用相同逻辑 offset；实际读取时会检查两个 buffer 得到的字节数一致。

### 6.2 static 和 dynamic

- **static**：根据 rank 数和文件大小计算连续 byte range，再向完整 FASTQ record 边界对齐。简单、RMA 少，但计算负载取决于 byte range 内的 read 难度。
- **dynamic**：把输入划分成较多 chunk，各 rank 通过 RMA ticket window 动态领取 chunk。主体使用大 chunk，默认最后 10% 进入四分之一大小的 medium tail，尾部最后若干 wave 再使用更小的 fine chunk。

dynamic 调度的关键代码在 `src/host/swbwa_mpi.c`：

- `swbwa_fastq_chunk_bytes()`：读取文件大小，估算 chunk bytes。
- `find_fastq_boundary()`：把估算 byte offset 调整到 FASTQ record 边界。
- `swbwa_mpi_fastq_scheduler_open()`：初始化调度器和 RMA ticket window。
- `swbwa_mpi_fastq_scheduler_next()`：领取 chunk、计算范围、记录本 rank 的 chunk 统计。
- `scheduler_chunk_range()`：根据 chunk id 计算主体/tail/fine-tail 的实际范围。

每次领取动态 chunk 时，当前算法大致做：

    ticket = Fetch_and_op(+1, target_rank_ticket)
    index  = target_rank + ticket * mpi_size
    range  = scheduler_chunk_range(index)

`next_queue = (queue + 1) % mpi_size` 只是让下一次尝试轮换目标 rank，避免一个 rank 连续成为所有请求的热点；它不改变 chunk 编号映射和正确性语义。

### 6.3 `MPI_EXACT_READ_INDEX` 和 `n_processed`

`n_processed` 是每个批次在全局输入中的逻辑 read 序号起点。它不是一个简单的“已处理多少条”的局部计数，部分 BWA-MEM 语义会使用这个全局 id：

- SE 中 `mem_mark_primary_se()` 以 `n_processed + i` 生成稳定 id/hash。
- PE 中 `mem_sam_pe()` 使用 pair 级的 `(n_processed >> 1) + i`。
- 某些 primary/tie/候选排序相关路径会间接受到该 id 影响。

`MPI_EXACT_READ_INDEX=1` 时，dynamic scheduler 在正式比对前建立完整 chunk/record index，获得精确的 chunk 首 record 前缀；这是正确性测试推荐的模式，代价是一次额外的 FASTQ 扫描。 `=0` 时，快速路径使用 byte offset 作为稳定但非精确的 seed/index，避免完整扫描，适合性能实验，但不能直接假设 SAM 与标准输出严格一致。

这也是为什么“输入 chunk 的字节数一样”并不能保证多个版本产生字节级相同 SAM：SAM 顺序、primary 选择、tie 语义和 `n_processed` 都可能参与其中。正确性验证应固定 `MPI_EXACT_READ_INDEX=1`，并使用统一的排序/标准化过程。

## 7. 输出模式和排序语义

统一输出接口位于：

    include/swbwa_output.h
    src/host/swbwa_output.c

上层只调用：

    swbwa_output_open()
    swbwa_output_write()
    swbwa_output_flush()
    swbwa_output_close()

### 7.1 `split`

每个 rank 输出独立文件：

    output.rank000000.sam
    output.rank000001.sam
    ...

不写 SAM header，使用批量缓冲。它避免了所有 rank 争抢单文件 offset，适合观察 stage2 和验证 MPI 计算负载。rank 文件按输入范围拼接并不一定恢复原始 SAM 顺序；correctness 脚本会按配置决定是否需要拼接或排序。

### 7.2 `single_unordered`

所有 rank 打开同一个文件。rank 0 建立并初始化一个全局 64 位 offset，其他 rank 在 flush 时通过：

    MPI_Fetch_and_op(... MPI_SUM ...)
    MPI_Win_flush(...)
    pwrite(... reserved_offset ...)

申请互不重叠的文件区间，然后各 rank 用 POSIX `pwrite` 并行写入。RMA 只负责“申请写入区间”，大块数据写入仍是并行的；但 Sunway 当前 MPICH/ch3:swch 的 RMA progress 需要目标 rank/主线程进入 MPI 的现象，会使 offset 申请和 flush 出现明显等待。因此这个模式的 stage3 和总时间可能大幅波动，不能把它的 I/O 时间混入 CPE kernel 性能结论。

该模式**不保证 SAM 记录顺序**，也不是当前的“单文件有序输出”。当前仓库没有完成的 ordered single-file 实现。

### 7.3 `discard`

不落盘 SAM，仅对输出 blob 做可选的 64 位 sum/XOR 指纹：

    OUTPUT_MODE=discard
    SWBWA_DISCARD_HASH=1

它适合隔离 stage2，尤其适合 MPI dynamic 的 chunk/RMA/CPE profiling。sum/XOR 满足结合律，可以由后处理脚本聚合，但它不是 MD5，也不能证明记录顺序和每一条 SAM 都完全相同。最终正确性仍然要跑完整 SAM 的规范化 MD5。

## 8. 内存和代码放置

### 8.1 主核 malloc wrapper

文件：

    include/malloc_wrap.h
    src/host/malloc_wrap.c

主核 wrapper 目前尽量简单地转发到系统 `malloc/free/realloc/strdup`。特殊点是 `calloc`：由于运行时加载器可以替换 malloc，但不一定替换 calloc，wrapper 用 `malloc + memset` 手工实现 calloc。 `HOST_MALLOC_STATS=1` 时收集主核分配统计，并在程序结束输出；默认关闭。

早期版本曾在业务代码中大量使用 `_sw_xmalloc/_sw_xfree` 和运行时动态判断 cross segment，后来统一回到普通 malloc/free，让加载器或 allocator 负责实际放置。不要重新在业务层引入成对但难以审计的特殊 free 路径。

### 8.2 CPE system allocator 和 pool allocator

`CPE_ALLOCATOR=system` 时，从核使用平台系统 allocator。

`CPE_ALLOCATOR=pool` 时，主要代码在：

    src/slave/malloc_wrap.c
    include/malloc_wrap.h

pool 按 4 B 到 256 KiB 的 size class 建立 segment tree，每个 CPE 有独立 pool slice 和独立 tree 元数据。大对象不适合进入 size class 时退回系统 malloc。pool 的特点是：

- 小对象分配/释放不需要频繁访问系统 allocator。
- 每个 CPE 的地址范围可快速判断 free 是否属于 pool。
- 释放只是更新 segment tree 的占用状态，不做通用堆的合并。
- pool 容量和 size class 由固定配置决定，超出容量或大对象路径必须关注日志和 fallback。

pool 与主核 wrapper 是两件事，不要把 `CPE_ALLOCATOR=pool` 理解为“所有主核 malloc 也进入 pool”。

### 8.3 `cgs_cross` 交叉段

交叉段相关接口和数据结构位于：

    include/swbwa_cpe.h
    include/swbwa_runtime.h
    src/host/bwamem.c
    src/slave/slave.c
    include/swbwa_cpe_layout.h

主核运行时会：

1. 分配交叉段目标 buffer。
2. 把 CPE text/data 拷贝到目标地址。
3. 根据 GOT、TLS、GP、CSR 等信息做必要的重定位/地址替换。
4. 使用 `swbwa_cpe_spawn()`/`swbwa_cpe_join()` 执行交叉段入口。

因此 cross 版本的错误常常表现为从核的 `UNALIGN`、`ACV1`、`OVPC1` 或跳转到入口附近，而不一定是算法错误。排查顺序应是：确认第二遍构建、确认 `swbwa_cpe_layout.h` 与可执行文件一致、确认 `-xmalloc/-cross_size` 资源参数、再看业务指针和 buffer 生命周期。

## 9. CPE 性能优化历史和当前结论

下面按“目标、实现、验证状态、论文价值”记录已经讨论和保留的优化。并不是所有后端都应作为默认生产路径。

### 9.1 LDM 和 BWT/Occ 热路径

**代码位置**：`src/slave/bwt.c`、`src/slave/bwamem.c`，相关配置在
`include/swbwa_config.h`，profiling 入口在 `src/slave/swbwa_cpe_profile.h`。

主要思路：

- 把小型、重复使用的 BWT/Occ 辅助数据和计数表尽量放到 LDM。
- 对 BWT/SA 访问增加软件预取，隐藏主存延迟。
- 让短 read 的 SMEM 临时向量优先使用局部空间，较长输入再走 heap。
- 避免把数 GB 的 BWT/SA 整体复制进 LDM；LDM 只保存元数据和可复用的热点小块。

profiling 中需要重点看 `MEM_CHAIN`、`MEM_COLLECT_*`、`CHAIN_BUILD_*`、GLD/GST、D-cache miss 和 memory wait。历史消融显示，SA batch、RID offset cache 等改动没有稳定收益；软件预取和局部小数据优化更值得保留，但应以多个数据集和多次运行的 stage2/CPE cycle 为依据。

**评价**：属于比较稳妥的 Sunway locality 优化，论文中可以作为“BWT-driven irregular memory access 的局部性处理”，但不要声称所有 BWT 数据都在 LDM。

### 9.2 每个 CPE 复用 worker context

**代码位置**：`src/slave/bwamem.c` 中的 `worker12_context_t` 及其初始化/销毁，
`src/slave/slave.c` 中的 `slave_worker12_s_pre_fast*()`。

生命周期不是“每个 read 一个 context”，而是：

    CPE kernel entry
      -> 一个 CPE 建立一个 worker context
      -> 该 context 处理本批次领取到的多个 read/pair
      -> CPE kernel exit 时销毁

context 主要保存：

- `bwt_t` 的小型元数据副本；
- SMEM 临时向量；
- chain seed arena；
- 与本 CPE 相关的复用状态。

原始路径中，部分 per-read alignment 会重复准备小型 `bwt_t`/辅助对象。当前优化把这类准备提升到 worker-batch 生命周期，减少每条 read 的小对象分配、初始化和元数据搬运。数 GB 的 BWT/SA 数组仍然在主存中，并没有被整体复制。

**评价**：这是批次级 runtime amortization，不是算法变化；适合和 CPE 动态调度、批处理一起描述。

### 9.3 24 KiB chain seed bump arena

**代码位置**：`src/slave/bwamem.c` 的 chain arena 辅助函数和 `mem_chain()`。

这是一个小型、临时、按 context 归属的 bump allocator：

1. 每条 read 开始时把 arena offset 复位为 0。
2. chain seed 数组优先从 arena 顺序分配。
3. 扩容时如果 arena 够用继续 bump；不够时退回原 heap 路径。
4. free 时区分 arena 地址和真正 heap 地址；arena 中的对象不逐块 free。
5. 下一条 read 开始时整体复位，因此不产生逐块回收开销。

它适合 chain seed 这种“同一条 read 内创建、read 完成后一起失效”的短生命周期对象。它不是 `CPE_ALLOCATOR=pool` 那种全局/每 CPE segment-tree pool，两者的生命周期和解决的问题不同。

**评价**：历史消融中比 SA batch/RID cache 更稳定，建议保留；论文中可归为“lifetime-aware temporary allocation”。

### 9.4 KSW query profile 和 DP 状态的 LDM 优先路径

**代码位置**：`src/slave/ksw.c` 的 `ksw_qinit_impl()`、`ksw_qdestroy()` 以及 query profile 使用路径。

对于尺寸不大的 query profile，如果估算内存不超过
`SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES`，且调用方允许，当前实现优先放入 LDM；否则退回 heap。销毁时用 `in_ldm` 等状态避免对 LDM 地址调用普通 free。

放入 LDM 的是 KSW query profile 和相关短期 DP 状态，不是整个 reference 或任意长度的 DP 矩阵。这样做的目的，是让 KSW 的内层循环减少对主存的重复读取。

### 9.5 KSW u8/i16 后端

**代码位置**：`src/slave/ksw.c`、`src/slave/scalar_sse.h`，选择宏在
`include/swbwa_config.h`。

当前后端：

    u8:  int32_16    默认、最稳的精确路径
         float16_16  实验路径
         float16_32  实验路径，可能更快但不保证 SAM 等价

    i16: scalar_8    标量基线
         int32_8     默认 i16 SIMD 路径

这里的“16”或“8”是逻辑 lane/向量组织方式；Sunway SIMD 没有直接等价于 x86 某些 int8/int16 packed primitive 的全部能力，所以当前 SIMD 实现常用 int32 元素保存窄整数逻辑值。这会浪费单条向量的位宽，但换来了已有编译器/指令模型下的可控正确性。

KSW 的路径选择仍由 score/overflow 条件决定：正常 75 bp/150 bp 数据在默认参数下通常主要走 u8；read 更长或 `-A` 等分数设置提高时可能进入 i16。调试 i16 时可用较高 `-A` 或实际长 read，但不要把人为触发 i16 的实验时间直接外推到普通 150 bp 数据。

### 9.6 150 bp same-PE forward 双状态 KSW

**代码位置**：`src/slave/ksw.c` 的 `ksw_align2_matesw_dual_forward()` 及相关配置。

只有同时满足以下条件才走快速路径：

- `SWBWA_ENABLE_MATESW_DUAL_FORWARD=1`；
- u8 后端为 `int32_16`；
- mate rescue 的两个 query 都是 150 bp；
- 两个 query/tlen 形状满足实现要求；
- 是 forward 方向的 same-PE 场景。

它使用两个 16-lane 状态组织一次批处理，然后对两个结果分别执行必要的 start/tie 处理，以保持原有语义。reverse 方向和不满足条件的输入回退到原始单次 KSW。

这个优化不是把任意两个 alignment 粗暴塞进 32 lane，也不是通用的 32-lane 等价实现；它是 workload-aware 的 guarded fast path。FP16 32-lane 实验不能直接作为严格 correctness 默认路径，因为浮点舍入、饱和和 tie 选择可能改变 SAM。

**评价**：这是最有“工作量/创新性”表现的 KSW 优化之一，但论文中必须报告触发率、fallback 比例和严格 MD5 结果，不能只展示理论 lane 数。

### 9.7 `ksw_extend2()` scratch 复用

**代码位置**：`src/slave/ksw.c` 的 `ksw_extend2()` 及其线程局部 scratch 管理。

当前为每个执行线程复用两类临时缓冲：`eh_t` 和 int8 query buffer；容量不足时才 `realloc` 扩容，正常调用复用已分配空间。它减少 repeated allocation、free 和初始化调用开销，不改变 extension 的 DP 递推语义。

**评价**：工程上稳妥，通常适合作为默认保留；性能增益应以真实 stage2 和 CPE cycle 说明，不要只看单个函数的微基准。

## 10. MPI 设计、RMA 和性能结论

### 10.1 MPI 初始化和线程级别

`src/host/swbwa_mpi.c` 当前调用：

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

如果运行时提供级别低于 `MPI_THREAD_MULTIPLE`，程序直接报错退出。原因是 dynamic input scheduler、single output RMA 和 progress thread 可能在不同主核线程中进入 MPI。当前 MPI-enabled 链接使用 `multi_static`，不是默认的 `single_static`。

`tests/test_mpi_env.c` 用于查看 MPI 标准版本、线程级别和窗口属性；`tests/test_mpi_rma.c` 用于隔离测试 RMA progress、双窗口和 Fetch-and-op 行为。它们不参与 SWBWA 默认构建。

### 10.2 RMA ticket window

每个 rank 暴露一个 64 位 ticket 变量，所有 rank 对目标 rank 的 ticket 做 `MPI_Fetch_and_op(... MPI_SUM ...)`。 `MPI_Win_lock/unlock` 或 lock-all/flush 的作用是建立 MPI-3 RMA 的访问/完成语义，不是普通 C 原子变量。

重要事实：

- `Fetch_and_op` 是原子的“读取旧值并累加”；多个 rank 不会拿到同一个 ticket。
- `MPI_Win_flush(target, win)` 确认对指定 target 的 RMA 操作已经完成并取得结果。
- 当前 Sunway MPICH `ch3:swch` 的测试表明，目标 rank 长时间不进入 MPI 时，RMA completion 可能延迟；这更像目标侧软件 progress 路径，而不能仅凭 MPI API 判定底层是否有硬件原子。
- x86 Open MPI 测试中 RMA latency 为微秒级；Sunway 测试中可到毫秒甚至更高，且会受到另一个 RMA window 和目标进程状态影响。

### 10.3 MPI progress thread

`swbwa_cpe_run_with_mpi_progress()` 在 CPE kernel 运行期间启动/配合 progress thread。线程使用 `MPI_Iprobe` 进入 MPI，并采用 absolute `clock_nanosleep` 控制间隔，避免相对 sleep 在 EINTR 下累计漂移。当前实现要求 `MPI_THREAD_MULTIPLE`。

这只能帮助 MPI runtime 获得进入机会，不能消除所有 RMA/文件系统瓶颈；`single_unordered` 的 output offset window 和 dynamic input window 仍可能竞争。调试时看：

- CPE completion-observing probe 时间；
- probe 次数和超过 1 ms/10 ms/100 ms/1 s 的次数；
- scheduler RMA lock/fetch/flush 时间；
- output flush 的 reservation 时间和实际 `pwrite` 时间。

### 10.4 负载不均衡的正确解释

曾经看到某些数据的 MPI rank stage2 差异很大，常见原因不是“所有 read 平均难度不同”这么简单，而是：

1. MPI chunk 仍然是不可再分的较大工作单元。
2. BWA-MEM 的 seed 数、重复区域、候选 chain 数和 mate rescue 数有长尾。
3. dynamic 只能在 chunk 领取阶段平衡，不能把一个已领取的慢 chunk 迁移给其他 rank。
4. 如果 stage1 读取慢，rank 领取新 chunk 的速度下降，它拿到的 chunk 数会变少，看起来像 stage2 不均衡。
5. 如果 output 是 single_unordered，RMA reservation 或共享文件系统写入慢，会改变 pipeline 背压，也会间接影响每个 rank 真正参与 stage2 的批次序列。

因此分析 MPI 时必须把 stage2、stage1、stage3 分开。 `correctness_results/timing_summary.md` 中 MPI 单元格用 `min-max` 表示 6 个 rank 的阶段范围，不能把一个范围误读成单个 rank 的平均时间。

历史上，`dynamic + split` 往往比 `dynamic + single_unordered` 更适合观察计算负载；single output 的慢点经常出现在 output/RMA，而不一定在 CPE kernel。 `ERR1203383` 之类的短 read 数据还可能具有不同的 seed/chain/mate 工作分布，不能只按输入文件字节数判断难度。

## 11. Profiling 和正确性验证

### 11.1 主核 stage timer

运行时 `-v 4` 会打开详细诊断；主核 timing report 由 `src/host/fastmap.c` 输出，包含：

    Pipeline
      total
      stage 1 - allocate and read raw FASTQ blocks
      stage 2 - align reads and generate SAM records
      stage 3 - write SAM and release batch data

    Stage 2 worker details
      part 1 - prepare CPE task and reusable buffers
      part 2 - CPE FASTQ formatting and input release
      part 3 - CPE alignment and SAM length pass
      part 4 - assign slices in the shared SAM buffer
      part 5 - CPE SAM record generation
      part 6 - release temporary worker data

stage1 还有 allocation、`fseeko`、`fread`、FASTQ boundary alignment 的调用次数、累计时间和最大单次时间。不要用 stage3 的慢写文件时间推断 CPE alignment 慢。

### 11.2 CPE LWPF profiling

编译时打开：

    ./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
        USE_MPI=0 CPE_PROFILE=1 CPE_PROFILE_CG=5

运行时再使用：

    ./run.sh cgs_cross -- ./SWBWA mem -v 4 -t 1 -1 \
        -o profiling.sam ref.fa reads.fq

`CPE_PROFILE` 和 `-v 4` 是两件事：前者编译 LWPF 插桩，后者控制运行时日志。profiling 只采样一个 CG，当前默认各个非 single 模式通常采样 CG5；需要横向比较时固定采样 CG 和节点条件。

profiling region 在 `include/swbwa_cpe_profile.h`，主要包括：

- `MEM_CHAIN`、`MEM_COLLECT_*`、`CHAIN_BUILD_*`；
- `CHAIN_FILTER`、`CHAIN_EXTENSION`、`CHAIN_EXTENSION_DP`；
- `MATE_RESCUE`、`MATE_REF_FETCH`、`MATE_KSW_ALIGN`；
- `KSW_QUERY_INIT_*`、`KSW_DP_*`；
- `MATE_DEDUP`、`DEDUP_SORT_*`；
- `SAM_FORMAT`、`SAM_COPY`。

关闭 `CPE_PROFILE` 时，profile helper 是空实现，避免把 profiling 成本带进普通构建。完整分析脚本在 `scripts/analyze_*.py`，profiling 结果和图表通常保存到 `correctness_results/` 的相应实验目录。

### 11.3 正确性矩阵

统一入口：

    bash scripts/correctness.sh matrix --io-mode both

该脚本支持 `run`、`verify`、`matrix`、`bigdata`、`bigdata-mpi`，使用 `.status` 下的 checkpoint 增量续跑。成功验证的 SAM 会删除，只保留轻量日志、MD5 和 TSV；失败 SAM 会保留以便诊断。

典型配置矩阵：

    single + system
    cgs_cross + pool
    MPI exact + static  + split
    MPI exact + static  + single_unordered
    MPI exact + dynamic + split
    MPI exact + dynamic + single_unordered

`has1/no1` 需要分别验证。MPI split、single_unordered 生成的记录顺序可能不同，必须根据脚本规则拼接/排序后再比较。 `correctness_results/timing_summary.md` 记录的一轮完整矩阵共有 72 个正确性检查通过。

当前小数据的 header-free normalized SAM MD5：

| 数据集 | PE | SE |
| --- | --- | --- |
| `SRR7963242` | `a799dc7268f389120ec92820e04b5118` | `0acfb1f46abd9fed3862c28a33e26da4` |
| `SRR2496709` | `20f980ce3955d09e7c133ba26ddf2b77` | `2ca59c689707cc0aa6e0ce1fe30ae145` |
| `ERR1203383` | `c2af4bf0b057d5125ce2a9d770f13741` | `473eec3972fbd651d8b911b9fb5c6e25` |

大数据标准值在 `scripts/bigdata_expected_md5.tsv`。这些值是经过统一 header 处理和必要排序后的参考值；不要把 MPI `discard` 的 sum/XOR 指纹称为 MD5。

## 12. 论文写作建议

### 12.1 最值得展开的创新点

下面这些方向既有实现工作量，也能和实验指标对应起来：

#### A. 面向 Sunway 异构存储层次的 BWA-MEM 映射

把 BWA-MEM 的 irregular seeding/chaining/SW 过程拆成主核流水线、CPE batch、LDM 热数据、主存 BWT/SA 和交叉段代码/运行时。重点不是简单地“使用 CPE”，而是说明不同数据的生命周期和访问规律分别放在哪里：

- 大型 reference index 保持主存共享；
- 小型 query profile、SMEM vector、chain seed 临时对象进入 LDM 或 CPE-local area；
- CPE text/private runtime 通过 cross segment 复制，降低单 CG 主存/代码热点；
- 主核利用流水线把输入、CPE 计算和输出重叠。

可报告：stage2、CPE region cycle、GLD/GST、cache miss、memory wait、single/cgs/cgs_cross 对比。

#### B. 两级动态调度

MPI rank 层使用动态 FASTQ chunk，CPE 层使用 read/pair 级动态任务；尾部采用 10% medium chunk + fine-tail wave，降低 MPI 最后不可分割慢 chunk 的影响。

这里要诚实区分两种粒度：MPI 不能迁移已领取的 chunk，CPE 才能在 chunk 内做 read 级平衡。论文可把它表述为“hierarchical scheduling”，并用以下数据支撑：每 rank 领取的 chunk id、chunk stage2 时间、CPE part3 累计时间、最大/最小 stage2、RMA 时间。

#### C. 保序的 CPE FASTQ 格式化和 SAM slice 生成

格式化阶段先统计各 CPE 的 record 数，再做 prefix sum，给每个 CPE 分配稳定的全局 read index 和 SAM buffer slice。这样 CPE 可以动态处理 read，同时最终 SAM slice 的空间分配不需要每条 read 通过主核原子追加。

这是一个很适合写成工程贡献的点：既解决了 CPE 并行格式化，又避免“谁先完成谁先写”带来的顺序不可验证问题。PE/SE 共用同一套框架，SE 通过第二 buffer 为空判断 pairing，不从不可靠的 pointer-embedded option 结构取状态。

#### D. 生命周期感知的 CPE memory hierarchy

把 worker context、24 KiB chain bump arena、KSW LDM query profile、CPE segment-tree pool 和主核 calloc wrapper 放在统一的 allocator/locality 设计中讨论。每个组件解决的问题不同：

- context：减少 batch 内重复元数据准备；
- arena：减少同一 read 生命周期内的短对象 free/metadata 操作；
- LDM profile：降低 DP inner loop 对主存的重复读取；
- pool：减少小对象系统 allocator 访问；
- host calloc：适配加载器只替换 malloc 的限制。

论文中应提供内存开销、fallback 率和 correctness 结果，避免只说“全部放进本地内存”。

#### E. workload-aware PE mate rescue KSW

150 bp same-PE forward 双状态路径是较有特色的优化：它只在形状严格匹配时触发，不满足条件就回退；因此既能利用常见测序 read 形状，又不会把任意 alignment 混装导致 tie 语义变化。建议报告：触发候选比例、forward/reverse/fallback 比例、mate rescue cycle、整体 stage2、严格 MD5。

### 12.2 应谨慎表述的方向

- `float16_16` 和 `float16_32` 是实验后端；不能把“某些输入更快”写成严格等价加速。
- `float16_32` 目前不保证 SAM-equivalent，适合消融，不适合默认 correctness 结论。
- `discard` hash 只是低成本指纹，不是完整正确性证明。
- MPI `single_unordered` 是并行写单文件，但不是有序输出。
- `MPI_EXACT_READ_INDEX=1` 会额外扫描输入，不能把这部分开销混入纯 alignment speedup；论文应分别报告 preprocessing 和 steady-state pipeline。
- stage1/stage3 的共享文件系统抖动不能当作 CPE kernel 变化。CPE 优化应主要看 has1、discard/split、stage2 和 CPE profile。

### 12.3 推荐的实验设计

每个关键配置至少重复 3 次；如果节点噪声明显，报告中位数、范围和异常节点说明。推荐固定：

1. 同一 reference/index、同一输入文件和同一 `-I`。
2. 同一 `-K`，并注明它是每 CG byte target。
3. 同一 `-t 1`、相同 `bsub` 资源；独占/共享节点不能混在一个表里。
4. correctness 使用 `MPI_EXACT_READ_INDEX=1`；性能 profiling 明确是否包含完整预扫描。
5. CPE 研究重点看 has1 + stage2；MPI I/O 研究再单独报告 stage1/stage3。
6. SIMD 消融至少包含默认 `int32_16`、候选后端、运行时间、CPE cycle 和 MD5。
7. MPI 消融至少包含 static/dynamic × split/single_unordered；split 是更干净的计算负载基线，single output 用来研究并行 I/O/RMA。

## 13. 已知限制和后续工作

当前维护者接手后应优先记住这些边界：

1. **没有 ordered single-file MPI output**：现有 `single_unordered` 只保证写区间不重叠，不保证全局 record 顺序。要实现有序单文件，最好从“每个逻辑 chunk 的 record byte length + prefix/extent + ordered commit”重新设计，不要在现有 RMA offset 上简单加全局锁。
2. **MPI dynamic 仍有 chunk tail**：尾部细分缓解但不能消除一个已领取慢 chunk 的长尾。更细粒度会增加 RMA 和 CPE setup 开销，`-K` 不能无限减小。
3. **Sunway RMA progress 是系统级瓶颈**：multi_static 解决线程级别可用性，但不保证 RMA 像 x86 Open MPI 一样微秒级或完全异步。
4. **CPE pool 有固定容量**：如果出现新的大对象/长 read，必须检查 pool exhaustion、fallback 和地址合法性。
5. **CPE profiling 目前按一个 CG 采样**：跨 CG 负载结论不能只根据单个采样 CG 推断。
6. **主核 pipeline queue 有固定数组容量**：`src/host/kthread.cpp` 中的 `DataType` queue 使用固定的 `1<<20` 条目数组；它不是按字节大小限制，但如果长期运行的 in-flight item 数超过设计上限，仍需重新评估计数器和环形索引的边界。
7. **构建依赖 Sunway 工具链**：`swgcc/mpicc/mpicxx/swreadelf/swaddr2line/swobjdump` 等工具在 x86 上不能直接复现最终 CPE 二进制；x86 只能用于独立 MPI/RMA demo 或算法对照。

## 14. 日常维护检查清单

### 改动源码前

- 先执行 `git status --short --branch`，不要覆盖用户未提交修改。
- 确认本次改动属于 host pipeline、CPE kernel、MPI、allocator 还是脚本；避免跨层混改。
- 如果改了编译宏，先在 `include/swbwa_config.h` 和 `Makefile` 找现有定义，避免添加同义 tag。

### 改动 CPE 代码后

- 先用 `EXEC_MODE=single` 或 `cgs` 做快速编译/小数据正确性。
- `cgs_cross` 必须重新跑 `build.sh` 两阶段流程。
- 用 `-v 4` 保留必要日志；长时间测试不要打开逐 read debug print。
- 先跑 `single + system` 基线，再跑 `cgs_cross + pool`，最后才跑 MPI。

### 做性能实验时

- stage2 与 stage1/stage3 分开记录。
- MPI 计算负载优先用 `OUTPUT_MODE=discard` 或 split；single_unordered 只作为 I/O/RMA 实验。
- 记录 `-K`、CG 数、rank 数、节点资源、has1/no1、exact index 和 profile 开关。
- 发现异常时先看 stage1 read bandwidth、stage3 write/RMA，再看 chunk stage2 和 CPE part3。
- 每次实验保留轻量 log/TSV，成功 correctness 后删除 SAM，避免结果目录失控。

### 提交前

    git diff --check
    git status --short
    git diff --stat

确认没有把节点生成的 `.sam`、临时 `.codex_*`、编译中间文件和个人路径提交进去。提交信息要说明：修改层次、开关变化、正确性验证、性能验证和已知限制。

## 15. 推荐阅读顺序

新接手者可以按下面顺序恢复上下文：

1. `README.md`：构建、运行和目录概览。
2. `include/swbwa_config.h`：所有当前配置宏和默认值。
3. `Makefile`：Make 变量如何映射成编译宏、MPI 如何链接 `multi_static`。
4. `src/host/fastmap.c`：命令行、输入、pipeline 和 stage timer。
5. `src/host/bwamem.c`：主核 worker/CPE 调用、`n_processed`、cross runtime。
6. `src/slave/slave.c`：CPE worker 入口、动态任务和格式化流程。
7. `src/slave/bwamem.c`：CPE BWA-MEM 主流程、context、arena、chain/extension。
8. `src/slave/bwamem_pair.c`：PE pairing、mate rescue、dedup 和 KSW 调用。
9. `src/slave/ksw.c`、`src/slave/scalar_sse.h`：KSW 后端和 SIMD 分支。
10. `src/host/swbwa_mpi.c`、`src/host/swbwa_output.c`：MPI 输入调度、RMA progress 和输出。
11. 本文第 9 节：当前 CPE 优化的统一说明；早期独立优化笔记已合并到本交接文档。
12. `correctness_results/timing_summary.md` 和 `scripts/README.md`：历史实验格式、MD5 和测试脚本约定。

这个项目最容易混淆的三件事是：`-K` 的 CG 级字节语义、MPI chunk 与 CPE read 任务的两级粒度、以及 single output 的 I/O/RMA 时间和 stage2 alignment 时间。后续分析只要始终把这三者拆开，很多“看起来像算法退化”的现象就能比较快定位到正确层次。
