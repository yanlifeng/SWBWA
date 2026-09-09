# Standalone Tests

这里保存不参与 SWBWA 默认构建的独立测试程序：

- `test_mpi_env.c`：MPI 线程级别和基础环境检查。
- `test_mpi_rma.c`：MPI RMA 原子 offset、progress 和双窗口诊断。
- `test_usleep_progress.c`：Sunway CPE progress/sleep 行为实验。

这些文件需要单独用平台 MPI 编译；它们不会被 `make` 自动链接进 `SWBWA`。
