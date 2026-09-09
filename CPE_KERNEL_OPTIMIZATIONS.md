# SWBWA CPE Kernel 算法与优化说明

本文档面向需要继续阅读、验证和优化 SWBWA 从核代码的开发者。内容以当前
`perf/cpe-kernel-profile-opt` 分支为准，重点解释：

1. BWA-MEM 在本项目中的实际执行流程；
2. 主核与 CPE 之间如何组织任务；
3. 已经保留的 CPE kernel 优化及其原理；
4. 各项优化对应的代码位置、开关、正确性边界和实测收益；
5. 已经验证无效或不满足严格正确性的方案，避免重复试错。

文中的“严格正确”指归一化 SAM 输出与基准结果一致，而不只是比对率或位置统计
接近。

## 1. 建议阅读顺序

如果准备逐项检查代码，建议按下面的顺序阅读：

| 顺序 | 主题 | 主要文件和函数 |
| --- | --- | --- |
| 1 | 主核三阶段流水线 | `bwamem.c::mem_process_seqs_merge2()` |
| 2 | CPE 动态任务调度 | `slave/slave.c::slave_worker12_s_pre_fast*()` |
| 3 | 单条 read 的 BWA-MEM 主流程 | `slave/bwamem.c::mem_align1_core_impl()` |
| 4 | SMEM、chain、chain extension | `mem_collect_intv()`、`mem_chain()`、`mem_chain2aln()` |
| 5 | PE 配对与 mate rescue | `slave/bwamem_pair.c::mem_sam_pe()` |
| 6 | KSW 后端 | `slave/ksw.c`、`slave/scalar_sse.h` |
| 7 | profiling 插桩 | `slave/swbwa_cpe_profile.[ch]` |
| 8 | 编译配置 | `swbwa_config.h`、`Makefile` |

## 2. BWA-MEM 算法回顾

### 2.1 总体思路

BWA-MEM 并不是直接对整条 read 和整条参考基因组做动态规划。它先利用 FM-index
寻找精确匹配种子，再把位置相容的种子组织成 chain，最后只在有限参考窗口内做
带状动态规划。这样把大多数工作从昂贵的 DP 转换为 FM-index 查询和局部扩展。

从一条 read 到 SAM 记录，大致经过：

```text
FASTQ read
  -> nt4 编码
  -> SMEM/seed 搜索
  -> seed 映射到参考坐标
  -> seed chaining
  -> chain 过滤
  -> chain 局部 DP 扩展
  -> 候选去重、主次比对判定
  -> PE 配对和 mate rescue（PE 才有）
  -> MAPQ、CIGAR、SAM 字段生成
```

### 2.2 单条 read 的入口

从核中的单 read 主入口是
[`slave/bwamem.c`](slave/bwamem.c) 的
`mem_align1_core_impl()`。当前代码的主要步骤是：

1. 将 ASCII 碱基转换成 BWA 使用的 nt4 编码；
2. 调用 `mem_chain()` 搜索种子并建立 chain；
3. 调用 `mem_chain_flt()` 过滤低质量或被包含的 chain；
4. 调用 `mem_flt_chained_seeds()` 进一步精简 chain 内种子；
5. 对保留的 chain 调用 `mem_chain2aln()` 做局部扩展；
6. 调用 `mem_sort_dedup_patch()` 排序、去冗余并修补候选。

这一段对应 BWA-MEM 的“找种子、连 chain、扩展成 alignment”主体，也是 CPE
计算的主要来源。

### 2.3 SMEM 和补充种子

`mem_chain()` 首先调用 `mem_collect_intv()` 收集 FM-index 区间。该函数包含三类
种子搜索：

1. `bwt_smem1()`：寻找 SMEM（super-maximal exact match）；
2. 对较长 SMEM 再做一次分裂搜索，补充被长 SMEM 遮盖的候选；
3. `bwt_seed_strategy1()`：以类似 LAST 的策略补充种子。

这些搜索最终得到 `bwtintv_t` 区间。区间仍是 FM-index 中的位置范围，不是最终
参考坐标。

### 2.4 从 FM-index 区间到 chain

`mem_chain()` 对种子调用 `bwt_sa()`，把 sampled suffix-array 位置恢复为参考坐标，
然后通过 `bns_intv2rid()` 确认所属 contig。位置和方向相容的种子被插入同一条
chain；当前实现使用 BWA 原有的 kbtree 结构寻找附近 chain。

chain 的核心含义是：“这些 seed 很可能来自同一条真实 alignment”。它减少了
后续 DP 的候选数量，也给 DP 提供参考窗口和对角线约束。

`mem_chain_flt()` 会根据 chain 分数、重叠关系和阈值去掉明显较差的 chain；
`mem_flt_chained_seeds()` 则进一步去掉 chain 内不值得单独扩展的种子。

### 2.5 chain 到局部 alignment

`mem_chain2aln()` 为每条保留 chain：

1. 根据 seed 和带宽估计参考区间；
2. 取出对应参考序列；
3. 对 chain 内种子排序；
4. 跳过已经被现有 alignment 覆盖的 seed；
5. 分别向 seed 左、右调用 `ksw_extend2()` 做带状扩展；
6. 生成 `mem_alnreg_t` 候选。

这里的 `ksw_extend2()` 是 chain extension DP。它与后文 mate rescue 使用的
`ksw_align2_matesw()` 不是同一个 kernel。当前 SIMD 后端优化主要集中在 mate
rescue KSW，不能把其收益直接理解成所有 DP 都得到同样加速。

### 2.6 PE 配对和 mate rescue

PE 数据的两个 read 各自完成上述单端 alignment 后，进入
[`slave/bwamem_pair.c`](slave/bwamem_pair.c) 的 `mem_sam_pe()`：

1. 根据一端的高质量 alignment 推测另一端应出现的参考窗口；
2. 在窗口中执行 mate rescue，即 `mem_matesw` 路径；
3. 调用 `mem_mark_primary_se()` 标记主、次 alignment；
4. 调用 `mem_pair()` 根据插入片段分布寻找最佳 read pair；
5. 计算配对分数和 MAPQ；
6. 调用 `mem_reg2aln()`、`mem_reg2sam()` 生成最终记录。

mate rescue 会对候选参考窗口做正向和反向局部 Smith-Waterman。其热点最终落到
[`slave/ksw.c`](slave/ksw.c) 的 `ksw_u8()` 或 `ksw_i16()`。

### 2.7 u8 和 i16 路径如何选择

`swbwa_matesw_run_one()` 按 BWA 原逻辑设置 `KSW_XBYTE`。简化后，条件是：

```text
read_length * match_score < 250  -> u8 路径
否则                            -> i16 路径
```

默认 `-A 1` 时，75 bp 和 150 bp 通常进入 u8；250 bp 及以上进入 i16。测试时可用
`-A 2` 强制 150 bp 数据进入 i16，用于验证 i16 SIMD 后端。

### 2.8 三类容易混淆的 DP

代码中至少有三类动态规划用途：

| DP 用途 | 主要函数 | 当前优化覆盖 |
| --- | --- | --- |
| chain 左右扩展 | `ksw_extend2()` | 保持原实现，局部 O3 等尝试未保留 |
| mate rescue 局部 SW | `ksw_align2_matesw()`、`ksw_u8()`、`ksw_i16()` | 本轮 SIMD/LDM/双状态优化的重点 |
| 最终 CIGAR 生成 | `mem_reg2aln()`、`bwa_gen_cigar2()` | 保持 BWA 语义 |

因此，看到 `KSW_U8_MODE` 或双状态 KSW 开关时，应理解为主要改变 mate rescue，
而不是替换整个 BWA-MEM 的 alignment 算法。

## 3. 主核与 CPE 的执行结构

### 3.1 主核三阶段流水线

主核入口位于 [`bwamem.c`](bwamem.c) 的 `mem_process_seqs_merge2()`。一批数据
大致分成：

1. stage 1：读取和格式化 FASTQ；
2. stage 2：CPE alignment、计算 SAM 长度和生成 SAM；
3. stage 3：写出 SAM 并释放批次数据。

在启用 CPE format 时，stage 2 又可以细分为：

1. 准备 CPE 参数和复用缓冲区；
2. CPE FASTQ 格式化；
3. CPE alignment 和 SAM 长度预计算；
4. 主核根据长度为每条记录分配精确 SAM slice；
5. CPE 将 SAM 写入已分配 slice；
6. 清理临时数据。

这里采用“两遍 SAM”设计：第一遍只确定每条记录的长度，主核做 prefix/slice
分配；第二遍由原 CPE 把文本写入自己的 slice。这样避免 CPE 对共享输出指针做
高频原子追加，同时可以维持批次内确定的记录布局。

### 3.2 CPE 动态调度

调度代码位于 [`slave/slave.c`](slave/slave.c)：

- `slave_worker12_s_pre_fast()`：普通 CGS 路径；
- `slave_worker12_s_pre_fast_cross()`：代码复制到交叉段的 CGS 路径；
- `slave_acquire_task()`：通过原子 fetch-and-add 领取任务；
- `slave_worker12_s_fast*()`：按第一遍记录的任务列表完成 SAM 写入。

当前 `SWBWA_READS_PER_DYNAMIC_TASK=1`。SE 时一个任务是一条 read，PE 时一个任务
是一对 read。细粒度调度能较好吸收 read 难度差异，尤其避免少量重复序列或复杂
read 把整个 CPE 拖成尾部任务。

每个 CPE 在第一遍保存自己领取过的 task ID。第二遍 SAM copy 不重新抢任务，而是
复用该列表，确保访问第一遍建立的数据和对应输出 slice。

## 4. Profiling 基础设施

### 4.1 为什么先做 profiling

主核 stage 2 计时只能看到整个 CPE kernel 用时，无法区分 BWT、chain、KSW、
pairing 或 SAM format。为此项目引入 LWPF3 手动插桩，但把开关集中在 Makefile，
关闭时调用会编译为空操作，不在热点源码中铺大量条件编译。

代码位置：

- [`slave/swbwa_cpe_profile.h`](slave/swbwa_cpe_profile.h)：从核 region 定义和空实现；
- [`swbwa_cpe_profile.c`](swbwa_cpe_profile.c)：主核侧 LWPF 配置、汇总和输出；
- [`Makefile`](Makefile)：`CPE_PROFILE`、`CPE_PROFILE_CG`；
- [`swbwa_config.h`](swbwa_config.h)：编译期配置检查。

### 4.2 region 父子关系

主要 region 可以按下面理解：

```text
WORKER_ALIGNMENT
  MEM_CHAIN
    SMEM_COLLECT
    CHAIN_BUILD
  CHAIN_FILTER
  CHAIN_EXTENSION
    CHAIN_EXTENSION_DP
  ALIGNMENT_FINALIZE
  SAM_FORMAT
    MATE_RESCUE
      MATE_REF_FETCH
      MATE_KSW_ALIGN
        KSW_QUERY_INIT_FORWARD
        KSW_DP_FORWARD
        KSW_QUERY_INIT_REVERSE
        KSW_DP_REVERSE
      MATE_DEDUP
        MATE_DEDUP_SORT_END
        MATE_DEDUP_REDUNDANCY
        MATE_DEDUP_SORT_SCORE
    PAIRING

SAM_COPY
```

`SAM_COPY` 在 alignment 主 region 之外，因为它属于第二遍 SAM 写入。

### 4.3 采集事件

当前 profiling 采集：

- cycle；
- retired instruction；
- 全局 load/store；
- D-cache access/miss；
- memory barrier/wait；
- instruction buffer empty 等平台事件。

profiling 默认关闭，只采一个指定 CG 的 64 个 CPE，并输出平均值、最小值和最大值。
这既控制日志量，也能观察 CPE 之间是否存在严重负载不均。

### 4.4 编译示例

生产模式：

```bash
bash build_cross.sh 8 \
  EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
  HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 \
  CPE_PROFILE=0 \
  KSW_U8_MODE=int32_16 KSW_I16_MODE=int32_8 \
  MATESW_DUAL_FORWARD=1 USE_MPI=0
```

profiling 模式：

```bash
bash build_cross.sh 8 \
  EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool \
  HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 \
  CPE_PROFILE=1 CPE_PROFILE_CG=5 \
  KSW_U8_MODE=int32_16 KSW_I16_MODE=int32_8 \
  MATESW_DUAL_FORWARD=1 USE_MPI=0
```

## 5. 已保留的 CPE kernel 优化

### 5.1 去掉超大静态 SAM 缓冲区清零

**位置**：[`bwamem.c`](bwamem.c)，`mem_process_seqs_merge2()` 的静态 SAM 缓冲区
初始化。

旧代码在初始化时清零三块约 512 MiB 的 SAM 缓冲区，总写流量约 1.5 GiB。后续
代码会按真实 SAM 长度完整覆盖使用区间，因此预先清零没有语义价值。

优化后只完成分配和必要元数据初始化，不再触碰整个缓冲区。这不是 CPE 算法本身
的优化，但会显著降低 stage 2 part 1，尤其是第一次批次。

在早期大数据实验中，ERR 数据的相关 part 1 从约 0.614 s 降到 0.097 s。

### 5.2 每个 CPE 复用 worker context

**位置**：

- [`slave/bwamem.c`](slave/bwamem.c)：`worker12_context_t`、context 初始化/销毁；
- [`slave/slave.c`](slave/slave.c)：`slave_worker12_s_pre_fast*()`。

每次外层 CPE kernel 启动时，每个 CPE 建立一个 context，并在它处理该次
`para->work_item_count` 批次的所有 read 期间复用。这里要区分两种粒度：

- 外层 batch：主核一次提交给 `worker12_s_pre_fast*()` 的工作批次，大小通常受
  `-K`、FASTQ 缓冲区和格式化结果影响；
- 内层 dynamic task：`slave_acquire_task()` 每次领取的一条 SE read 或一对 PE reads。

context 是“每个 CPE、每个外层 batch”一个，不是“每个 dynamic task”一个。context
主要保存：

- `bwt_t` 的小型元数据副本；
- SMEM 临时向量；
- chain seed arena；
- 与本 CPE 相关的复用状态。

这里只复制一次 `bwt_t` 元数据，数 GB 的 BWT 数组仍在主存。这里最重要的收益是
避免每条 read 重复分配、复制和释放 `bwt_t`。需要特别注意：优化前的 `worker12`
路径已经通过 `w->aux[tid]` 复用了部分 SMEM 辅助对象，因此不能把所有 SMEM
初始化开销都归因于这一次 context 改造；SMEM 临时 vector 放入 LDM 是后面的独立
优化。

生命周期是：

```text
CPE kernel entry
  -> worker context init
  -> 动态领取多个 read/pair task
  -> worker context destroy
CPE kernel exit
```

原来的 `worker12` 热路径中，每条 read 都进入 `mem_align1_core()`，由它执行：

```text
ldm_malloc(sizeof(bwt_t))
复制 bwt_t、L2[5] 和 cnt_table[256]
执行 mem_align1_core_impl()
ldm_free(bwt_t)
```

当前 `worker12_pre_fast()` 直接调用 `mem_align1_core_impl()`，把同一个 CPE context
中的 `&context->bwt` 和 `context->smem_aux` 传给每条 read。于是上述 `bwt_t` 的
分配、复制和释放从“每条 read 一次”变成“每个外层 batch 一次”。

当前的 `mem_align1_core()` 包装函数仍然保留，并且仍然包含上述 per-read `bwt_t`
生命周期；它服务于 `worker1_pre_fast()` 等非 `worker12` 路径。也就是说，context
优化主要覆盖当前批量 alignment 的 `worker12` 热路径，并没有把整个项目所有
alignment 入口都改成 context 方式。

### 5.3 SMEM 临时向量放入 LDM

**位置**：[`slave/bwamem.c`](slave/bwamem.c) 的 `smem_aux_t`、`smem_aux_init()` 和
`mem_collect_intv()`。

SMEM 搜索需要两个临时 `bwtintv_t` vector。常见 75/100/150 bp read 的临时元素
上界与 query 长度同量级，因此 context 为每个 vector 预留 256 项，总计约 16 KiB。

处理长度小于 256 的 read 时直接使用 LDM 固定数组，避免动态扩容和主存小块访问；
更长 read 保留原 heap 路径，避免把支持长度写死为 256 bp。

这项优化的主要价值是降低短 read 高频临时分配和缓存抖动。单独消融时收益存在
波动，因此应和 context、chain arena 一起理解，而不是当作独立的大幅加速项。

### 5.4 24 KiB chain seed bump arena

**位置**：[`slave/bwamem.c`](slave/bwamem.c) 的 chain arena 辅助函数和
`mem_chain()`。

chain build 会创建并扩展许多短生命周期 seed 数组。系统 allocator 或 CPE pool
都需要维护元数据，并可能造成主存访问。当前每个 CPE context 预留 24 KiB arena：

1. 每条 read 开始时将 arena offset 复位为 0；
2. chain seed 初始存储和扩容优先从 arena 顺序分配；
3. arena 不足时退回原 heap 路径；
4. free 时区分 arena 地址和真正 heap 地址。

这是 bump allocator，不单独回收每一块；一条 read 完成后整体复位。它适合 chain
seed 这种生命周期一致的小对象。它和 worker context 是配套关系，但 arena 本身
是单独引入的优化，不能和 `bwt_t` per-read 复制消除混为一谈。

BWT 子项消融中，SA batch、RID offset 缓存等方案没有稳定收益，最终只保留了该
arena，因为它对 allocator 开销的改善最稳定，且不改变算法语义。

### 5.5 KSW query profile 和 DP 状态优先放入 LDM

**位置**：[`slave/ksw.c`](slave/ksw.c) 的 `ksw_qinit_impl()`、`ksw_qdestroy()`。

KSW 会为 query 构造 striped profile，并分配 `H0/H1/E/Hmax` 等 DP 状态。旧路径
存在较粗的容量估算和频繁 CPE heap 分配。当前实现：

1. 精确计算 query profile 与 DP 状态总字节数；
2. 在一个连续、按 SIMD 对齐的区域内布置所有数组；
3. mate rescue 临时 query 在容量允许时优先放入 LDM；
4. 超过设定容量时使用 CPE heap；
5. `ksw_qdestroy()` 根据来源正确释放。

典型 150 bp mate rescue 的实际需求约为数 KiB，明显小于旧估算。连续布局也减少
多个 allocator 调用并改善地址局部性。

注意：当前 LDM 选择主要按计算出的容量判断，并依赖已验证的 128 KiB 手动 LDM
预算。它不是任意压力下都能自动从一次失败的 `ldm_malloc()` 无缝重试 heap，因此
修改 context 常驻空间时必须重新核算峰值。

### 5.6 mate rescue 反向互补序列复用

**位置**：[`slave/bwamem_pair.c`](slave/bwamem_pair.c) 的 mate rescue task 准备、
运行和 finish 路径。

同一条待 rescue read 可能对应多个候选 anchor。原流程可能为每个候选重复构造
reverse-complement query。当前 task 在第一次需要反向 query 时延迟生成 `task->rev`，
后续候选复用，task 结束后统一释放。

优化不改变候选顺序、DP 参数或 tie 规则，只消除完全相同的数据变换。

### 5.7 mate rescue 候选统一去重

**位置**：[`slave/bwamem_pair.c`](slave/bwamem_pair.c) 的
`swbwa_matesw_finish()`。

旧组织方式会在处理不同方向/候选后多次调用去重。当前先收集本次 mate rescue
产生的完整候选集，最后只执行一次 `mem_sort_dedup_patch()`。

最终去重函数、比较规则和候选集合没有改变，因此这是等价的调用合并。profiling
中 mate dedup 是明显热点，这项改动减少了重复排序和全记录搬运。

### 5.8 mate rescue 专用紧凑 key 排序

**位置**：[`slave/bwamem.c`](slave/bwamem.c) 的 `mem_sort_dedup_patch()` 及其紧凑
排序辅助函数。

`mem_alnreg_t` 较大，直接排序会反复移动完整结构体。mate rescue 的去重调用满足
`bns == NULL && pac == NULL && query == NULL`，不执行普通路径的 patch，因此可以
安全选择专用排序：

1. 构造约 24 字节的 key，包含坐标、分数、query 起点和原下标；
2. 先按参考终点相关 key 排序；
3. 按结果一次性重排完整 `mem_alnreg_t`；
4. 执行原有 redundancy 判定；
5. 再按分数/坐标 key 排序并重排。

key 和 scratch 在总容量不超过约 64 KiB 时放入 LDM，否则使用 heap。普通
`mem_sort_dedup_patch()` 调用仍走原实现，避免改变 chain alignment 的 patch 和 tie
语义。

这是一轮优化中收益最大的单项之一。150 bp 代表数据的优化阶梯中，加入紧凑 key
后 stage 2 从约 52.6 s 降到 47.4 s，scratch 进一步 LDM 化后约为 46.4 s。

### 5.9 只对 mate-rescue KSW 热函数使用 O3

**位置**：[`slave/ksw.c`](slave/ksw.c) 的 `SWBWA_MATESW_HOT`，作用于 `ksw_u8()`
和 `ksw_i16()`。

整个 `slave/ksw.c` 使用 O3 曾导致 SE 路径回退，说明编译器对 `ksw_extend2()` 等
函数的调度或代码尺寸不适合统一提高优化级别。当前只给 profiling 明确识别出的
mate-rescue SW 热函数添加函数级 O3，其余代码保持工程默认 O2。

这种做法限制了代码尺寸和寄存器压力的影响范围。代表数据中仅热函数 O3 有稳定
收益，而整个文件 O3 已撤销。

### 5.10 u8 KSW 的三种后端

**位置**：

- [`swbwa_config.h`](swbwa_config.h)：模式常量和默认值；
- [`slave/scalar_sse.h`](slave/scalar_sse.h)：Sunway SIMD 兼容操作；
- [`slave/ksw.c`](slave/ksw.c)：query profile 和 `ksw_u8()`。

Sunway CPE SIMD 没有直接等价于 x86 SSE 的 16 路 unsigned int8 饱和运算。项目
保留三种实现用于比较：

#### int32_16（默认）

使用 16 个 32-bit lane 模拟原算法的 16 个 unsigned byte lane。辅助操作显式完成
饱和加减、比较、max、移位和 tie 扫描。

优点是 stripe 宽度、Lazy-F 传播顺序和 tie 行为可以与原 16-lane 语义对应，当前
严格正确性最好。缺点是一个 SIMD 寄存器只承载 16 个有效低精度元素，硬件利用率
不理想。

#### float16_16（实验）

使用硬件 FP16 lane 表示整数分数，但只使用低 16 个逻辑 lane，高半区屏蔽。它
保留 16-lane stripe 和 tie 语义，正确性测试通过，但 FP16 转换、屏蔽和模拟饱和
操作使大数据实测慢于 `int32_16`。

#### float16_32（实验，不保证严格正确）

使用全部 32 个 FP16 lane。理论吞吐更高，但 stripe 宽度从 16 变为 32，会改变
Lazy-F 的传播路径、循环停止时刻以及相同分数下的 tie 选择。部分数据输出一致，
部分数据出现 MD5 差异，因此不能作为严格正确的默认后端。

编译选择：

```bash
make KSW_U8_MODE=int32_16
make KSW_U8_MODE=float16_16
make KSW_U8_MODE=float16_32
```

大数据对比结果如下：

| 后端 | ERR stage 2 | small SRR7963242 stage 2 | SRR2496709 stage 2 | 严格正确性 |
| --- | ---: | ---: | ---: | --- |
| int32_16 | 62.674 s | 182.222 s | 60.550 s | 全部通过 |
| float16_16 | 71.563 s | 189.064 s | 63.759 s | 全部通过 |
| float16_32 | 62.061 s | 170.304 s | 59.365 s | 两组出现差异 |

因此当前默认明确回到 `int32_16`。FP16 后端保留用于后续架构研究，不应在正式
正确性实验中误开。

### 5.11 150 bp same-PE forward 双状态 16+16 KSW

**位置**：

- [`slave/bwamem_pair.c`](slave/bwamem_pair.c)：same-PE 候选匹配和双任务组织；
- [`slave/ksw.c`](slave/ksw.c)：`ksw_align2_matesw_dual_forward()` 及双状态 DP；
- [`swbwa_config.h`](swbwa_config.h)：`SWBWA_ENABLE_MATESW_DUAL_FORWARD`；
- [`Makefile`](Makefile)：`MATESW_DUAL_FORWARD`。

#### 动机

FP16x32 直接把单个 alignment 的 stripe 扩成 32 lane 会改变算法语义。更安全的
利用方式是：低 16 lane 计算一个 alignment，高 16 lane 同时计算另一个完全独立
的 alignment。每个 alignment 内仍保持原 16-lane stripe。

profiling 显示，150 bp PE 数据中，同一对 read 的两个方向能配成双任务的候选约
占 57.61%，且 forward DP 的 `qlen/tlen` 兼容率很高，因此有实际合并空间。

#### 实现

`mem_sam_pe()` 按两个 mate-rescue 方向的同序候选尝试配对。满足条件时：

1. 低 16 lane 保存方向 0 的 H/E/F/profile；
2. 高 16 lane 保存方向 1 的独立状态；
3. 两半分别维护 `gmax`、`te`、peak、active 和 Lazy-F 状态；
4. lane shift 使用掩码，禁止状态跨越 16-lane 边界；
5. 一个方向提前结束后可冻结，另一个继续；
6. tie 扫描顺序仍按各自 16-lane 原路径执行。

当前只合并 forward DP。reverse 起点搜索继续执行两个原始串行调用，这样限制实现
复杂度和正确性风险。

#### 启用条件和回退

运行时同时满足以下条件才走双状态路径：

- 编译期开启 `MATESW_DUAL_FORWARD=1`；
- u8 后端为严格的 `int32_16`；
- 两个任务都设置 `KSW_XBYTE`；
- 两条 query 都是 150 bp；
- 两任务 `qlen` 相同且 `tlen` 相同。

任何条件不满足都回退为两个原 KSW 调用。因此 75 bp、i16、尺寸不兼容以及其他
后端不会被强行合并。

#### 效果

150 bp 代表数据的完整正确性检查通过。多次运行平均：

- stage 2：46.294 s -> 44.262 s，约 4.39% 改善；
- CPE part 3：45.005 s -> 43.299 s，约 3.79% 改善。

这低于理论上界，因为只有 forward、只有兼容候选能合并，BWT、chain、reverse、
SAM 等工作完全不变。

### 5.12 i16 的 int32_8 SIMD 后端

**位置**：

- [`slave/scalar_sse.h`](slave/scalar_sse.h)：8-lane i16 兼容辅助操作；
- [`slave/ksw.c`](slave/ksw.c)：`ksw_i16()`；
- [`swbwa_config.h`](swbwa_config.h)：i16 模式；
- [`Makefile`](Makefile)：`KSW_I16_MODE`。

原 i16 KSW 的逻辑 stripe 是 8 个 signed 16-bit lane。Sunway 没有完全对应的
指令，因此 `int32_8` 使用低 8 个 32-bit lane 表示原 8 个 i16 lane，高 8 lane
清零。辅助函数显式实现：

- 有符号 i16 加法并 clamp 到 `[-32768, 32767]`；
- 减法并按 KSW 语义 clamp；
- max 和比较；
- 8-lane shift；
- max reduction 和 tie 扫描。

与 u8 的 16->32 lane 尝试不同，这里逻辑 stripe 仍是原来的 8，因此可以保持原
DP 传播和 tie 语义。

可选模式：

```bash
make KSW_I16_MODE=scalar_8
make KSW_I16_MODE=int32_8
```

用 `-A 2` 强制 150 bp 数据进入 i16 后，完整 SE 正确性检查通过。代表测试中：

- stage 2：13.280 s -> 约 12.843 s，约 3.29% 改善；
- CPE part 3：12.650 s -> 约 12.074 s，约 4.56% 改善。

当前默认是 `int32_8`。

## 6. LDM 预算和使用约束

当前常用配置把 256 KiB CPE 本地空间中的 128 KiB 交给编译器 cache，剩余部分供
手动 LDM 和栈使用。主要手动对象约为：

| 对象 | 生命周期 | 典型/上限 |
| --- | --- | ---: |
| SMEM 两个固定 vector | 整个 worker context | 约 16 KiB |
| chain seed arena | 整个 worker context | 24 KiB |
| context 其他元数据 | 整个 worker context | 约 1 KiB |
| KSW query/profile/DP | 单次 mate rescue | 常见数 KiB，上限约 64 KiB |
| dedup key/scratch | 单次 dedup | 上限约 64 KiB |

KSW 临时区和 dedup scratch 不会在同一调用点同时活跃，但它们都与 context 常驻的
约 41 KiB 共存。峰值接近 100 KiB 以上，还要给栈和运行库留余量。

因此后续增加 LDM 缓冲时必须遵守：

1. 先画出生命周期重叠，而不是只把各函数的局部大小分别相加；
2. 优先复用互斥阶段的 scratch；
3. 长 read 和大候选集合必须保留 heap 路径；
4. 改变 `-cache_size` 后重新核算手动 LDM；
5. 不能假设所有 `ldm_malloc()` 失败点都已有自动 heap 重试。

## 7. 优化效果汇总

### 7.1 第一阶段热点优化的大数据结果

以下数据来自
`correctness_results/cpe_ldm_opt_20260826/analysis/bigdata_stage2.tsv`，是 context、
LDM、紧凑排序等第一阶段优化完成后的历史检查点，不包含后续所有双状态/i16
实验，因此不能和后面的百分比简单相乘。

| 数据 | stage 2 基线 | stage 2 优化后 | 加速 | part 3 基线 | part 3 优化后 | 加速 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 75 bp PE | 69.469 s | 61.445 s | 1.131x | 65.233 s | 56.911 s | 1.146x |
| small SRR7963242 150 bp PE | 216.244 s | 182.751 s | 1.183x | 212.954 s | 179.141 s | 1.189x |
| small SRR7963242 150 bp SE | 54.837 s | 52.008 s | 1.054x | - | - | - |
| SRR2496709 PE | 67.978 s | 60.089 s | 1.131x | - | - | - |

PE 收益高于 SE，符合热点来源：mate rescue、mate dedup 和 KSW LDM 优化主要只在
PE 路径活跃。

### 7.2 当前默认组合

当前生产默认应理解为：

| 项目 | 默认 |
| --- | --- |
| u8 KSW | `int32_16`，严格 16-lane 语义 |
| i16 KSW | `int32_8` |
| same-PE dual forward | 开启，但仅兼容 150 bp u8 候选实际进入 |
| CPE profiling | 关闭 |
| CPE 调度 | 每次原子领取一条 SE read 或一对 PE reads |
| context/SMEM/chain arena | 启用，长 read 保留原路径 |
| compact mate dedup | 启用，仅 mate-rescue 专用调用进入 |
| FP16 后端 | 保留但不默认启用 |

## 8. 已尝试但未保留的优化

### 8.1 BWT SA batch

尝试过一次批量恢复 4/8/16 个 SA 位置，希望提高访存并行度。局部 cycle 有时下降，
但端到端收益不稳定，并与 chain arena 等改动存在交互，最终撤销。

BWT/Occ 查询的地址具有数据依赖：下一步位置依赖当前查询结果。单纯扩大 batch 会
增加活跃状态、寄存器和随机访存压力，不一定提高实际吞吐。

实验记录：

- `correctness_results/bwt_sub_ablation_20260826/README.md`；
- `correctness_results/bwt_sa_arena_interaction_20260826/README.md`。

### 8.2 BWT prefetch、RID offset LDM 和 prefetch distance

这些方案在个别 region 的 cycle 上有改善，但没有形成稳定、可叠加的 stage 2
收益。BWT 约数 GB，访问随机且下一地址相关；软件预取容易过早、过晚或污染有限
cache。`bns_intv2rid()` 的小型 offset 缓存也不是主导热点。最终均删除。

### 8.3 用 DMA/RMA 搬 BWT 热点

DMA 适合连续、可提前确定的大块数据。BWT Occ/SA 恢复是细粒度随机访问，地址链
依赖明显，DMA 启动和等待成本通常高于一次普通缓存访问。参考窗口 fetch 相对只占
约 0.3% 的热点，也不足以支撑复杂 DMA 双缓冲。因此当前没有为 BWT 引入 DMA/RMA。

### 8.4 整个 ksw.c 使用 O3

全文件 O3 使 SE 路径回退，可能来自代码尺寸、寄存器压力或 `ksw_extend2()` 指令
调度变化。现只保留 `ksw_u8()`/`ksw_i16()` 的函数级 O3。

### 8.5 专用小数组排序

曾考虑对 mate dedup 使用插入排序或固定小网络。profiling 显示候选数并不总是小：
150 bp 数据中 65 个以上候选占比接近一半。最终选择“紧凑 key + 通用排序”，而不
采用只对极小 n 有利的专用 sorter。

### 8.6 FP16 32-lane 单 alignment

该方案确实能降低部分 cycle，也在个别数据上更快，但改变 stripe 和 Lazy-F/tie
语义，严格 MD5 不稳定。代码作为实验后端保留，默认禁用。

## 9. 正确性边界与审查重点

逐项审查时，建议特别检查下面几类风险：

1. **tie 语义**：KSW 相同分数下选哪个位置会影响 CIGAR、XA/SA 和最终 MD5；
2. **lane 边界**：双状态 16+16 的 shift 不能让低半状态进入高半；
3. **饱和语义**：int32 模拟 u8/i16 时必须在每个原本会饱和的操作后 clamp；
4. **对象来源**：LDM、arena、CPE pool、系统 heap 的 free 路径必须匹配；
5. **长 read 回退**：固定 256 项 SMEM 和 24 KiB arena 不能成为输入上限；
6. **普通 dedup 路径**：紧凑排序只能用于无 patch 的 mate-rescue 调用；
7. **PE 双任务匹配**：不兼容的 qlen/tlen 或非 u8 候选必须走串行原实现；
8. **LDM 峰值**：新增临时区要和 context 常驻区、栈以及 cache 划分一起计算。

## 10. 后续优化建议

基于现有 profiling 和消融结果，后续更值得投入的方向是：

1. 继续细化 same-PE 双任务覆盖率，寻找安全合并 reverse DP 的方法；
2. 对 i16 的真实长 read 数据做 profile，而不只依靠 `-A 2` 强制路径；
3. 在不改变 16-lane tie 语义的前提下改进 int32 模拟指令调度；
4. 复用 KSW 与 dedup 互斥 scratch，进一步释放 LDM 余量；
5. 分离“局部 cycle 改善”和“stage 2 端到端改善”，所有新方案至少做两次大数据
   消融并保留严格 MD5 检查。

不建议优先重做已经被否定的 BWT 大块 DMA、盲目预取、全文件 O3 或 FP16x32 单
alignment。这些方向已经证明局部指标好看不等于端到端、严格正确的收益。

## 11. 相关提交与实验目录

当前优化演进的主要提交：

| 提交 | 内容 |
| --- | --- |
| `7c3ebe6` | CPE 热点 profiling、mate rescue 复用、去重和 SAM 初始化优化 |
| `b68a473` | KSW LDM、SMEM/context、chain arena 和紧凑排序 |
| `cdd53fe` | 显式 u8 后端选择与 FP16 实验实现 |
| `05b705b` | same-PE forward 双状态 KSW 和 i16 SIMD |

建议配合阅读的实验记录：

- [`correctness_results/cpe_hotspot_opt_20260825/README.md`](correctness_results/cpe_hotspot_opt_20260825/README.md)
- [`correctness_results/cpe_ldm_opt_20260826/README.md`](correctness_results/cpe_ldm_opt_20260826/README.md)
- [`correctness_results/bwt_sub_ablation_20260826/README.md`](correctness_results/bwt_sub_ablation_20260826/README.md)
- [`correctness_results/bwt_sa_arena_interaction_20260826/README.md`](correctness_results/bwt_sa_arena_interaction_20260826/README.md)

这些目录中的数据来自不同优化检查点。做横向比较时应确认可执行文件提交、编译
开关、输入数据和 `-K/-A/-1` 参数一致，避免把多轮实验的比例直接相乘。
