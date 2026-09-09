# p2 `ksw_extend2()` 分布调研

## 目的

这次实验针对“是否值得把主 chain extension 改成 8-lane `int32` SIMD”做
优化前调研。只在 p2 的独立副本
`/home/user_home/ylf/someGit/bwa_ksw_extend_profile_20260831` 中修改，原始
`/home/user_home/ylf/someGit/bwa` 未修改。

插桩位置是副本的 `ksw.c:ksw_extend2()`。每次调用只累计一次实际访问的 DP
band cell 数，并每 64 次调用抽取一条样本，输出到 TSV。没有改变 DP recurrence、
边界判断或返回值。

## 输入和运行方式

参考索引使用：

`/home/bigssd/ylf_data/BWA_index/GRCh38.d1.vd1.fa`

每个数据取前 400,000 行 FASTQ，即 100,000 条 read；PE 两端取相同前缀，
同时跑了 PE 和 R1-only 的 SE。程序使用 `-t 32 -K 5000000`，SAM 输出到
`/dev/null`，统计报告保存在远端副本的 `profile_runs/reports/`。

构建插桩版本：

```sh
make clean
make -j4 EXTENSION_PROFILE=1
```

离线分析：

```sh
python3 scripts/analyze_ksw_extend2_profile.py profile.tsv
```

## 结果摘要

| 数据集 | 模式 | `ksw_extend2` 调用数 | 平均 DP cells/调用 | DP cells 中位数 | p90 | p99 | 8-call 无重排利用率中位数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 | PE | 451,066 | 492 | 219 | 1,350 | 2,109 | 0.374 |
| ERR1203383 | SE | 214,440 | 487 | 222 | 1,362 | 2,077 | 0.374 |
| SRR2496709 | PE | 710,528 | 814 | 362 | 2,420 | 4,599 | 0.381 |
| SRR2496709 | SE | 339,599 | 757 | 355 | 2,199 | 4,502 | 0.374 |
| SRR7963242 | PE | 4,273,042 | 2,844 | 1,216 | 8,420 | 11,509 | 0.339 |
| SRR7963242 | SE | 1,682,797 | 3,172 | 1,314 | 8,950 | 11,169 | 0.358 |

长度分布也不是原始 read 长度的简单复制：

- ERR：`qlen` 中位数约 25，p99 约 55，最大 56；
- SRR2496709：`qlen` 中位数约 40--41，p99 约 80，最大 81；
- SRR7963242：`qlen` 中位数约 60--66，p90 约 110，p99 约 129，最大 131。

所有调用的请求 band width 都是 100，但实际有效 width 会被 query 长度和
gap 约束截短，因此 `effective_w` 大致跟 `qlen` 走。

## 对 8-lane SIMD 的判断

目前不适合直接把连续到达的 8 个 extension 调用塞进一个 8-lane kernel：

- 三个数据的无重排 8-call `max/min` 中位数约为 17--72 倍；
- 对应 lane 利用率中位数只有约 0.34--0.38；
- 只有 SRR2496709 的极少数样本组能达到 `max/min <= 1.5`，ERR 和
  SRR7963242 基本没有；
- 将样本按实际 DP work 离线排序后利用率接近 1，但这只是理论上界，不能
  直接作为在线调度方案。

因此，结论不是“8-lane SIMD 永远没有收益”，而是：**单纯在现有调用顺序
外面套一个 8-lane kernel，收益很可能被空 lane 和最长 lane 等待抵消。** 若继续
做，应先增加按 `dp_cells` 或粗粒度 work class 分桶的队列，并验证分桶、搬运和
结果回填成本；否则优先优化单调用的 scalar/SIMD 内核更稳妥。

## 正确性和开销校验

同一份 100,000-read ERR PE 前缀在基线和插桩二进制之间进行了比较。单线程下
SAM 正文完全一致，`cmp` 的差异只有 `@PG CL` 中二进制名称不同导致的命令行
字段变化。单线程运行时间为 11.096 s（基线）和 11.084 s（插桩），说明这次
统计的额外开销在该样本上不可见；多线程版本也只有同样的 `@PG` 头部差异。

## 文件

- `ERR1203383_PE.tsv`、`ERR1203383_SE.tsv`
- `SRR2496709_PE.tsv`、`SRR2496709_SE.tsv`
- `SRR7963242_PE.tsv`、`SRR7963242_SE.tsv`
- `summary.txt`：脚本生成的完整文本摘要
- `logs/`：六次采样运行和基线/插桩校验的轻量日志
- `../scripts/analyze_ksw_extend2_profile.py`：离线分析脚本
