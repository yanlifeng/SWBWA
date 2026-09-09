# Source Layout

- `host/`：主核 C/C++ 实现，包括 FASTQ 驱动、BWA-MEM 流程、MPI 输入/输出和主核辅助模块。
- `../slave/`：CPE 侧实现。该目录保持原位置，因为神威交叉编译和跨段链接规则直接依赖它。
- `../include`：当前项目的公共头文件仍位于仓库根目录，以兼容主核/从核同名头文件和现有 include 写法。
