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
不提交远端作业，不产生 SAM，也不能代替 Sunway 的 cross 两遍编译和比对回归。
