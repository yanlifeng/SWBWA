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
- `test_ksw_modes_compile.py`：不同整数、FP16 及双状态配置的条件编译。

这些测试使用本机编译器向量适配层及 SDK stub，不执行神威指令；FP16
组合只检查编译，不能据此宣称 FP16 数值正确。运行前需要 Python 3 和
支持 ASan/UBSan 的 C 编译器，测试临时产物由各测试自行清理。

## LDM mode checks

`test_ldm_modes.py` checks the three production `CPE_LDM_MODE` presets,
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
