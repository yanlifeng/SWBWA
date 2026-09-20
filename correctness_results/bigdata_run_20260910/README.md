# Big-data 正确性与性能回归 — 2026-09-10

## 测试环境

- 远端：`ylf@121.250.210.19`（跳板机）→ `swls-transfer`（Sunway 登陆节点 `sw_hpc_64.qd.sw`）
- 部署目录：`guoshi/ylf/SWBWA_20260910/`（源码打包上传，非 git clone；本地基线为 `master@a3e7beb` + `scripts/correctness.sh` 的 `bigdata --configs/--io-mode` 增量补丁）
- 数据：`guoshi/ylf/data/bwa_test_big_data/{ERR1203383,small_SRR7963242,SRR2496709}_{1,2}.fastq`
- 参考基因组：`guoshi/ylf/data/GRCh38.d1.vd1.fa`
- 队列：`q_share`，全程单节点、单点占用；本地始终只通过嵌套 SSH 发命令并每 60s~8min（按历史耗时估算）probe 一次，不在登陆节点跑 agent
- 测试范围：仅 `has1`（`-1` 单线程流水线路径），未测 `no1`
- I/O 限速：全部传输使用 `rsync --bwlimit=500`；上传源码包 464KB，回传结果包 435KB（远小于约定的 10MB/20MB 上限）

## 测试矩阵与结果（24/24 全部 PASS）

| 配置 | 数据 | PE | SE |
| --- | --- | --- | --- |
| single + system | ERR1203383 | PASS | PASS |
| single + system | small_SRR7963242 | PASS | PASS |
| single + system | SRR2496709 | PASS | PASS |
| cgs_cross + pool | ERR1203383 | PASS | PASS |
| cgs_cross + pool | small_SRR7963242 | PASS | PASS |
| cgs_cross + pool | SRR2496709 | PASS | PASS |
| mpi + dynamic + split | ERR1203383 | PASS | PASS |
| mpi + dynamic + split | small_SRR7963242 | PASS | PASS |
| mpi + dynamic + split | SRR2496709 | PASS | PASS |
| mpi + dynamic + single_unordered | ERR1203383 | PASS | PASS |
| mpi + dynamic + single_unordered | small_SRR7963242 | PASS | PASS |
| mpi + dynamic + single_unordered | SRR2496709 | PASS | PASS |

所有 MD5 均与 `scripts/bigdata_expected_md5.tsv` 的标准值逐位匹配（split/single_unordered 均已按脚本规则合并/排序后再比较）。详见 `bigdata_results/md5_results.tsv`、`bigdata_mpi_results/md5_results.tsv`。校验通过的 SAM 均已在远端删除，仅保留日志与 `.status` checkpoint。

## 性能对比（stage2 = CPE 比对与 SAM 生成，单位秒）

### 非 MPI（单 rank 完整时间）

| 配置 | 数据 | PE stage2 (本次) | PE stage2 (历史 08-21) | SE stage2 (本次) | SE stage2 (历史) |
| --- | --- | ---: | ---: | ---: | ---: |
| single_system | ERR1203383 | 303.29 | 351.03 (**-13.6%**) | 117.10 | 129.11 (**-9.3%**) |
| single_system | small_SRR7963242 | 968.25 | 1275.18 (**-24.1%**) | 261.96 | 277.62 (**-5.6%**) |
| single_system | SRR2496709 | 300.17 | 348.65 (**-13.9%**) | 119.19 | 130.56 (**-8.7%**) |
| cgs_cross_pool | ERR1203383 | 61.40 | 71.09 (**-13.6%**) | 27.40 | 29.28 (**-6.4%**) |
| cgs_cross_pool | small_SRR7963242 | 172.52 | 220.99 (**-21.9%**) | 49.95 | 55.64 (**-10.2%**) |
| cgs_cross_pool | SRR2496709 | 59.08 | 69.62 (**-15.1%**) | 26.33 | 29.02 (**-9.3%**) |

历史参照：`correctness_results/bigdata_results/timing_summary.tsv`（2026-08-21 采集）。本次实测的 CPE 比对阶段（stage2）相较历史基线普遍快 5%~24%，与 `DEVELOPMENT_HANDOFF_ZH.md` 中记录的 LDM/worker context/chain arena/KSW query profile 等优化的累计效果一致。完整 stage1/stage2/stage3 明细见 `nonmpi_timing_20260910.tsv`。

### MPI dynamic（6 rank，min-max 范围，单位秒）

| 配置 | 数据 | PE total | PE stage2 | SE total | SE stage2 |
| --- | --- | --- | --- | --- | --- |
| mpi_dynamic_split | ERR1203383 | 148.53-148.54 | 56.16-56.90 | 68.75-68.76 | 22.41-23.68 |
| mpi_dynamic_split | small_SRR7963242 | 240.50-240.50 | 165.57-172.93 | 87.08-87.09 | 47.04-50.22 |
| mpi_dynamic_split | SRR2496709 | 145.76-145.77 | 53.50-58.12 | 72.36-72.37 | 22.27-25.07 |
| mpi_dynamic_single | ERR1203383 | 168.60-168.61 | 51.55-60.03 | 78.25-78.25 | 22.14-24.22 |
| mpi_dynamic_single | small_SRR7963242 | 225.81-225.82 | 166.92-172.43 | 87.91-87.92 | 44.87-51.02 |
| mpi_dynamic_single | SRR2496709 | 162.41-162.42 | 51.09-57.36 | 78.17-78.17 | 21.02-24.28 |

与历史 `correctness_results/timing_summary.md`（小数据集 SRR7963242/SRR2496709/ERR1203383，非本次 bigdata 尺寸）不完全同源，仅作为量级参考；本次 MPI stage2 min-max 范围与 `mpi + dynamic` 的历史结论一致：`split` 输出通常比 `single_unordered` 更干净地反映 CPU 计算负载，`single_unordered` 的 stage3（RMA offset + pwrite）耗时更长且波动更大（如 ERR1203383 PE stage3 达 79.79-90.75s，对比 split 的 54.82-59.18s）。完整明细见 `mpi_timing_20260910.tsv`。

## 文件清单

- `bigdata_driver.log` / `bigdata_mpi_driver.log`：驱动脚本完整输出（含每个 case 的 dispatch/PASS 记录）
- `bigdata_results/`、`bigdata_mpi_results/`：各配置的 `build.log`、`*.run.log`、`*.sort.log`、`*.sorted.md5`、`.status/*.md5.ok`、`md5_results.tsv`
- `nonmpi_timing_20260910.tsv`：single_system / cgs_cross_pool 的 stage1/2/3/total 明细
- `mpi_timing_20260910.tsv`：mpi_dynamic_split / mpi_dynamic_single 的 6-rank min-max 明细

## 未覆盖范围

- 仅测 `has1`；`no1`（三阶段 pipeline_queue 路径）未在本次覆盖
- `mpi + static`（static input）未在本次覆盖，用户仅要求 dynamic
- 小数据集矩阵（SRR7963242/SRR2496709/ERR1203383 的 `../data` 版本，非 bigdata）未重新验证，沿用既有 `correctness_results/timing_summary.md`
