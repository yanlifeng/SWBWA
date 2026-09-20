# Scripts

- `correctness.sh`：编译、运行、断点续跑和 MD5 校验入口。
- `md5_concat_sam.sh`、`check_discard_hash.py`：split/discard 输出的校验辅助工具。
- `analyze_*.py`：读取日志/TSV 并生成 profiling、MPI 调度和 CPE 优化分析。
- `sorted_sam.cpp`：用于正确性比较的 SAM 排序工具源码。
- `benchmark_ksw_u8_modes.sh`、`run_mpi_discard_profile.sh`：专项性能实验驱动脚本。
- `ab_profile/`：stage2 kernel 级 old/new A/B 剖析工具（构建、`bsub` 驱动、限速回传、解析与渲染）。见 `ab_profile/README.md`。
- `deploy_ab/`：A/B 部署、discard 性能和尾块实验的驱动/解析工具，见该目录 README。
- `prep_data.sh`：经批准后在计算节点生成 tiny/big4 输入，需要显式数据目录与空间确认。
- `startup_driver3.sh`：A/B 启动耗时实验，使用显式 `WORK_ROOT`、`DATA_ROOT`，默认一个节点。
- `bigdata_expected_md5.tsv`：大数据集的标准 MD5 对照表。

脚本默认从仓库根目录执行，结果写入 `correctness_results/`。平台相关的作业参数集中在脚本开头的环境变量中，不要把节点生成的 SAM 或缓存文件提交进 Git。
