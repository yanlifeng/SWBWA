# Source Layout

- `host/`：主核 C/C++ 实现，包括 FASTQ 驱动、BWA-MEM 流程、MPI 输入/输出和主核辅助模块。
- `slave/`：CPE 侧实现，包括 CPE kernel、从核辅助算法和从核专用头文件。

主核和从核使用各自目录中的同名头文件；`include/` 中只放跨目标共享的主核接口和
配置头文件。构建产物（`.o`、`.d`）保留在对应源代码目录中，由 `make clean` 清理。
