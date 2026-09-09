# Shared Headers

这里保存主核编译使用的公共头文件、配置头文件和跨模块接口。

- `swbwa_config.h`：编译期开关和模式定义。
- `swbwa_cpe_layout.h`：由 `build.sh` 在跨段构建时生成的 CPE 布局信息。
- 其余头文件：BWA-MEM、MPI、输出、内存包装和基础数据结构接口。

CPE 专用的同名头文件放在 `src/slave/`，避免主核和从核的编译接口互相覆盖。
