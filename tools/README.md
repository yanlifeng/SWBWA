# Build Tools

- `xlink.py`：从链接后的可执行文件提取跨段执行所需的符号和地址信息。
- `get_tls.py`：提取构建阶段需要的线程局部存储布局信息。

上述两个脚本由 `build.sh` 自动调用，通常不需要单独执行。

- `ldm_policy.py`：可选的离线LDM准入参数选择工具（C），不由构建或程序启动自动调用。
  默认运行使用内置B策略。输入格式和限制见 [LDM_POLICY.md](../docs/LDM_POLICY.md)。
