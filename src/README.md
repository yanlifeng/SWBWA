# Source Layout

- `host/`：主核 C/C++ 实现，包括 FASTQ 驱动、SWBWA 比对流程、MPI 输入/输出和主核辅助模块。
- `slave/`：CPE 侧实现，包括 CPE kernel、从核辅助算法和从核专用头文件。

主核和从核使用各自目录中的同名头文件；`include/` 中只放跨目标共享的主核接口和
配置头文件。构建产物（`.o`、`.d`）保留在对应源代码目录中，由 `make clean` 清理。

程序入口只有 `host/main.c`，主核算法库为 `libswbwa.a`。CPE 的入口由
`slave/slave.c` 提供，不编入独立 CLI 或示例程序。继承的 `bwa_*` / `mem_*`
内部接口仍保留，避免无关的跨主从核符号变更；来源与许可见根目录 `NOTICE.md`。
