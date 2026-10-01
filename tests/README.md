# Standalone Tests

这里保存不参与 SWBWA 默认构建的独立测试程序：

- `test_mpi_env.c`：MPI 线程级别和基础环境检查。
- `test_mpi_rma.c`：MPI RMA 原子 offset、progress 和双窗口诊断。
- `test_usleep_progress.c`：Sunway CPE progress/sleep 行为实验。

这些文件需要单独用平台 MPI 编译；它们不会被 `make` 自动链接进 `SWBWA`。

## 本机回归检查

`bash tests/run_host_checks.sh` 使用系统 C 编译器和 ASan/UBSan，测试真实的
CPE allocator（仅 SDK 的 LDM 接口由 `stubs/slave.h` 模拟），覆盖预算拒绝、
heap 回退前提、按批统计、pool 增长/扩容及错误路径；同时检查 discard 哈希解析。
不提交远端作业，不生成数据集的 SAM 结果；输出测试只在临时目录生成小型样本并清理。
这些检查不能代替 Sunway 的 cross 两遍编译和比对回归。

测试不链接到正式程序。生产路径保留错误/边界检查与可选性能统计；不要为了
精简代码删除 SAM 校验或 allocator 所有权检查。已移除实验分支对应的独立测试
也一并移除；`test_optimization_defaults.py` 额外检查三个执行模式的构建对象清单
及 `libswbwa.a` 链接名称，防止重新引入重复 CPE CLI/示例。

新增的输出检查覆盖非 MPI discard、完整哈希、正常 SAM 输出，以及
`CPE_DISCARD_DIGEST` 的 SE/PE 元数据生命周期、边界和错误状态。启用 CPE
digest 时检查的是生成的 SAM 字节，不包含被跳过的最终 SAM 拷贝。
`test_optimization_defaults.py` 检查 Make 默认开关、显式关闭及不支持的
配置组合，不需要 Sunway 工具链。

`bash tests/run_cpe_kernel_checks.sh` 运行 CPE 核心优化的差分检查：

- `test_ksw_fused_gap*.py`、`test_ksw_xor_select_native.py`：整数 SIMD 基元的等价性。
- `test_ksw_u8_native.py`、`test_ksw_i16_native.py`：真实 DP kernel 的返回值、状态和边界条件。
- `test_ksw_extend2_qp_native.py`：五符号 query profile 的生成顺序和结果。
- `test_ksw_extend2_ldm_native.py`：extension LDM 与 heap 路径、预算拒绝和释放。
- `test_chain_reuse_smem_ldm.py`：SMEM/chain arena 复用的容量和所有权。
- `test_matesw_dedup.py`：与主核原版方向循环差分比较 rescue/去重调用顺序，
  覆盖未新增命中、无效区间、多方向及双状态长度回退；不代替真实 SAM 回归。
- `test_ksw_modes_compile.py`：不同整数、FP16 及双状态配置的条件编译。

这些测试使用本机编译器向量适配层及 SDK stub，不执行神威指令；FP16
组合只检查编译，不能据此宣称 FP16 数值正确。运行前需要 Python 3 和
支持 ASan/UBSan 的 C 编译器，测试临时产物由各测试自行清理。

## LDM mode checks

`test_ldm_modes.py` checks the four production `CPE_LDM_MODE` presets,
SDK allocation gating and rejection of retired/incompatible options.
`test_ldm_allocator.py` retains isolated historical-policy regression coverage
and actual global-DP score/CIGAR checks under sanitizers. Both are included in
`run_cpe_kernel_checks.sh`. These native checks do not replace two-pass Sunway
builds and complete output validation.

`test_ldm_persistence.c` is a separate Sunway-only MPE/CPE probe for allocation,
data/TLS persistence and stack interference across ordinary CGS spawn/join
calls. It is not part of the native test runner or SWBWA build, and does not
validate the custom cross-segment runtime's PC/SP handling. Keep it as a
hardware diagnostic, not as evidence that arbitrary cross-runtime lifetimes
or large LDM arenas are safe.

Unified-pool checks (native SDK stubs, not a replacement for Sunway runs):

- `test_unified_ldm.py`: profile/growth switches, scratch coexistence, capacity,
  fallback and actual global-DP/CIGAR differential tests under sanitizers.
- `test_pool_bitmap.py`: lowest-free-slot agreement with the original tree,
  metadata-cache overflow, address lookup, SDK refusal and suspend/resume.
- `test_ldm_policy_allocator.py`: production B and optional C policies using
  the real allocator, including random realloc and data-retention checks.
- `test_ldm_policy_config.py`: built-in B defaults agree with the offline tool;
  malformed/unsafe environment settings fail instead of silently changing mode.
- `test_ldm_policy.py`: offline selection, node normalization and refusal of
  incomplete, incorrect, profiled or mixed-node measurements.
- `test_sam_stream_md5.py`: exact header-free bytes, arbitrary write boundaries,
  malformed/truncated input, >4 GiB counters and no-file output integration.

`python3 tests/test_output_rma.py` uses a native MPI installation (`MPICC` and
`MPIEXEC` may select it) and three local ranks. It compares normal and
`OUTPUT_RMA_ONLY=1` writer buffering/reservations, verifies complete-SAM samples
for the first 100 reads per chunk, oversized writes, an idle rank, and absence
of a disk output. It is separate from the non-MPI host checks.

## Cross Runtime Stack Regression

`python3 tests/test_cross_stack.py RUN_LOG` checks the `-v 4` cross-dispatch
trace: consecutive task IDs, three completed phases per batch, and stable
post-return stack ranges across all CPEs. With no arguments it tests the
log parser locally. Run the hardware check with at least 300 batches (a small
`-K` is useful), in addition to the normal chunk-size run. Compare full SAM
hashes against a `cgs` build using the same input and chunk configuration;
stable stack traces alone are not a correctness check. Cross binaries and
relocation files must come from the complete two-pass `build.sh` workflow.
Use `--min-batches 300` for the stress acceptance check. The parser requires
identical per-core extrema and an aligned LDM stack address (`bsub -b`); it
checks post-return stack accumulation, not peak stack usage inside a kernel.
# Chunk Output And MPE Helpers

`test_parallel_input.py` checks exact positioned input with 1--6 readers,
unaligned slices, empty input, short reads, `EINTR`, EOF, errors and offsets
beyond 4 GiB using ASan/UBSan. It does not measure Sunway I/O performance.

Six-MPE builds default to six `pread` readers; one-MPE builds default to serial
`fread`. `SWBWA_INPUT_READERS=0` explicitly selects serial `fread`. Runtime values
1--6 select `pread` with that many readers, limited by the build's
`HOST_MPE_THREADS` (default 6 for `cgs`/`cgs_cross`, 1 for `single`).
`test_optimization_defaults.py` checks Make/header defaults and explicit overrides;
`test_parallel_input.py` exercises reader selection and input under both defaults
and overrides. Multiple readers require `HOST_MPE_THREADS=6` and a
full-chip allocation. Byte slices are reassembled before FASTQ parsing;
MPI chunk IDs, boundaries and ticket reservations are unchanged.

`test_pipeline_queue.cpp` exercises the bounded, condition-variable queues
with empty input, repeated reuse, slow readers and slow writers. The test
checks FIFO delivery, one final flush and SAM ring ownership under ASan/UBSan.
Output capacity is reserved before Stage 2 starts using a ring slot.
`test_host_workers.c` also checks concurrent callers and a held input task:
short SAM preparation/packing tasks can run on their caller instead of waiting
for that input task. Serial fallback preserves the six logical byte slices so
measurement and packing can safely take different dispatch paths.

At `-v 4`, Stage 1 reports thread user/system CPU, faults and context switches
when `RUSAGE_THREAD` is available. `FASTQ read wall time` includes dispatch
and waiting; `pread calls (worker time sum)` sums concurrent calls and must
not be added to wall time. Changing defaults does not establish performance
gains for every storage environment; retain explicit overrides for comparisons.
In pipelined mode, queue wait times are reported separately from stage callback
times; concurrent waits overlap and must not be summed as elapsed time.

See [CHUNK_OUTPUT.md](CHUNK_OUTPUT.md) for the batch output contract,
`HOST_MPE_THREADS=1|6`, and local/target-platform validation.
## Ordered Output

`python3 tests/test_output_ordered.py` runs real local-MPI tests of the
experimental `single_ordered` writer. See [ORDERED_OUTPUT.md](ORDERED_OUTPUT.md)
for prefix dependencies, monotonic scheduling, bounded buffering and limitations.

`check_ordered_sam.py --fastq reads.fq --reference reference.sam --candidate output.sam`
compares **all** SAM bytes after excluding headers, without sorting. Add `--paired`
for PE input. It also checks FASTQ name order, primary alignment multiplicity,
record count, byte count and MD5; supplementary records are retained. This small
fixture checker expects unique FASTQ names (both mates share one name). Run it
where the data reside; SAM files need not be transferred for validation.
`test_ordered_sam_check.py` verifies that reordered, truncated, duplicated and
modified records are rejected. A non-MPI SWBWA reference checks output-path
consistency, not independent agreement with upstream BWA.
