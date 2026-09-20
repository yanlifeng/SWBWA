# `ab_profile/`：stage2 kernel 级 old/new A/B 剖析工具

这套脚本用于回答一个具体问题：**从基线提交到当前 HEAD 的 CPE 优化，在 stage2 上究竟买到了什么，以及现在的时间花在哪里。**
它把“干净 stage2 墙钟时间”和“LWPF 每区域 cycle 计数”分成两组构建，在**同一集群、同一 `bsub` 窗口内背靠背**运行，用 ABBA 顺序抵消漂移。

## 实验设计

| | 提交 | 角色 |
| --- | --- | --- |
| `old` | `fc658b8` | 产生历史基线的提交 |
| `new` | `a3e7beb` | 当前 `HEAD` |

每个 variant 构建 4 个二进制：`{single(system), cgs_cross(pool)} × {CPE_PROFILE=0,1}`。
每个条件按 `old_np → new_np → new_p → old_p` 顺序跑 4 个阻塞式 `bsub`（`np` = 无 LWPF，量 stage2 墙钟；`p` = `CPE_PROFILE=1`，取 kernel cycle）。共 2 配置 × 3 数据 × PE/SE = 12 条件，全程同时最多占 1 个节点。

`-1`（has1）贯穿始终，因此 stage1/stage3 不与 stage2 重叠；stage2 包含主核准备、
CPE 格式化/比对/SAM 生成及等待，并非纯 CPE 指令执行时间。

## 前置条件

```text
$AB_ROOT/old    fc658b8 的源码 checkout
$AB_ROOT/new    当前 HEAD（a3e7beb）的源码 checkout
```

两个 checkout 都必须是能通过 Sunway 工具链构建的完整树（含 `tools/`）。`./ab_build.sh` 会自行创建 `logs/`、`runs/`、`scratch/`。
表中的 `new=a3e7beb` 是这次归档实验的固定版本，不代表以后执行时的最新 HEAD。
执行前显式设置 `AB_ROOT` 和 `DATA_ROOT`，不再默认访问旧账号目录；参考
[`../../TESTING_GUIDE_ZH.md`](../../TESTING_GUIDE_ZH.md) 检查 `gfsquota`。
复用已有构建目录前请手动确认版本；构建脚本不会覆盖已有 runs。
driver 和 MD5 检查共用文件锁；提交状态未知时会停止，必须先用 `bjobs` 确认，不能直接重启。

## 使用顺序

在 Sunway 登陆节点上执行（不要在登陆节点跑 agent；本目录脚本只发 `bsub`）：

```bash
# 1) 构建 8 套二进制 + 独立运行目录（每个目录自带 SWBWA + data.bin + data.bin2）
./ab_build.sh

# 2) 先跑一个最便宜的条件验证整条链路
./ab_driver.sh smoke

# 3) 跑完整 12 条件（可重复执行，已完成的条件会跳过）
nohup ./ab_driver.sh > logs/driver.out 2>&1 &

# 4) 正确性抽检：CPE_PROFILE=1 的构建必须复现参考 MD5
./md5_check.sh cross ERR1203383 PE dc5c0a6babd41641db22808caedb7a44
```

本地侧拉日志（新账号直连，限速 500 KB/s，预检最多 10 MiB，不含 SAM、不使用 `--delete`）：

```bash
AB_REMOTE=/path/to/completed/experiment bash fetch_logs.sh
```

解析与渲染：

```bash
python3 parse_ab.py <logs_root> <out_dir>     # 生成 4 个 TSV
python3 region_tree.py <out_dir>/kernel_cycles.tsv new_p   # 层级占比（不要逐行相加）
python3 report.py <out_dir> cross             # 渲染 markdown 表
```

## 关键约束（踩过的坑）

1. **`./SWBWA` 必须是 `bsub` 的直接子命令。** 神威 CPE/MPE launcher 只给直接子进程建立 CPE 映射；用 `bsub ... bash script.sh` 包一层会让 SWBWA 变成孙进程，启动即 SIGSEGV。这是之前整轮实验全部失败的真正原因，不是代码 bug。
2. **`SWBWA` 从当前目录读 `data.bin` / `data.bin2`**（见 `src/host/bwamem.c` 的 `load_cpe_relocation`）。所以每个运行目录都必须自带这两个文件，否则启动即死。
3. **裸跑 `./SWBWA`（无参数）在已知正常的构建上也会 segfault**，所以“无参数崩”不能作为构建损坏的证据。
4. **队列拥塞**：提交可能立刻失败并回 `No enough compute nodes`，脚本按 `RETRIES`/`RETRY_SLEEP` 退避重试；真正跑过但返回非零的作业不会被盲目重试。
5. **只回传日志**。计时作业成功后删除 SAM（这些不算 MD5 验证）；独立 MD5 作业只在 PASS 后
   删除 SAM。失败或完成状态不明时保留现场并停止。
6. **LWPF region 是嵌套的**，父区域的计数包含子区域。只有 `WORKER_ALIGNMENT` 和 `SAM_COPY` 互斥，两者之和才是 CPE 总时间；`SAM_FORMAT` 是 `WORKER_ALIGNMENT` 的子区域，`MATE_RESCUE`/`PAIRING` 又在 `SAM_FORMAT` 里面。不要把这些行相加。

## 结果去向

一次完整运行的产物已归档在 `correctness_results/profile_ab_20260913/`（含 `README.md`、`ab_stage2.tsv`、`kernel_diff.tsv`、`kernel_new_breakdown.tsv`、`kernel_cycles.tsv`、`report_cross.md`、`report_single.md`）。
