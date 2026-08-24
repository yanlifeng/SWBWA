# SWBWA 从核性能剖析结果

## 实验设置

- 运行模式：`-t 1 -1`，仅分析 has1 的 Stage2。
- 数据集：75 bp 的 `ERR1203383 PE` 和 150 bp 的 `SRR7963242 PE`。
- 配置矩阵：`cgs/cgs_cross` 与 `system/pool` 的 2 x 2 组合。
- 采样范围：CG5 的 64 个从核。普通 CGS 和 cross CGS 使用相同采样位置。
- 每个 LWPF 单元格依次为 64 个从核的平均值、最小值和最大值。
- PMU 指标：周期、指令数、全局读写、D-cache 访问/缺失、内存屏障等待周期、取指缓冲为空周期。

各区域存在嵌套关系。例如 `WORKER_ALIGNMENT` 包含 `SAM_FORMAT`，后者又包含
`MEM_CHAIN`、`MATE_RESCUE` 等路径，因此区域周期不能相加。启用 PMU 插桩也会
带来额外开销，下面的绝对时间用于同一插桩版本之间的比较。

## Stage2 时间

| 配置 | ERR1203383 PE | 相对 cgs/system | SRR7963242 PE | 相对 cgs/system |
| --- | ---: | ---: | ---: | ---: |
| cgs + system | 8.738 s | 1.00x | 133.708 s | 1.00x |
| cgs + pool | 8.233 s | 0.94x | 58.186 s | 0.44x |
| cgs_cross + system | 5.848 s | 0.67x | 127.293 s | 0.95x |
| cgs_cross + pool | 4.873 s | 0.56x | 56.397 s | 0.42x |

## 关键 PMU 结果

以下数值均为 `WORKER_ALIGNMENT` 在采样 CG 上的每从核平均值。

| 数据集 | 配置 | 周期 | 指令 | 近似 CPI | 内存屏障等待 | 取指缓冲为空 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ERR1203383 | cgs + system | 17.0G | 3.04G | 5.60 | 1.39M | 12.0G |
| ERR1203383 | cgs + pool | 16.0G | 3.27G | 4.89 | 10 | 11.0G |
| ERR1203383 | cgs_cross + system | 10.0G | 3.08G | 3.24 | 1.63M | 1.72G |
| ERR1203383 | cgs_cross + pool | 8.23G | 3.27G | 2.52 | 11 | 2.50G |
| SRR7963242 | cgs + system | 289G | 73G | 3.96 | 826M | 184G |
| SRR7963242 | cgs + pool | 118G | 74G | 1.59 | 49 | 18G |
| SRR7963242 | cgs_cross + system | 273G | 72G | 3.79 | 830M | 170G |
| SRR7963242 | cgs_cross + pool | 114G | 74G | 1.54 | 49 | 14G |

## 结论

1. `SRR7963242` 上 pool 是主要收益来源。四种配置的动态指令数基本一致，
   但 system allocator 下总周期、内存屏障等待和取指缓冲为空周期大幅增加。
   这说明慢速主要来自等待和拥堵，而不是执行了更多比对工作。
2. 单独 cross 对长读帮助有限：Stage2 只改善约 4.8%，PMU 周期也仅从
   289G 降到 273G。pool 开启后再叠加 cross，可从 58.186 s 小幅降到
   56.397 s。
3. `ERR1203383` 的分配压力较低，pool 单独收益只有约 5.8%；cross 对它更
   有效，显著降低了取指缓冲为空周期。cross 与 pool 叠加后的 Stage2 最快。
4. 因此不同数据集表现不一致是合理的：短读更容易暴露代码放置/取指路径，
   长读的 mate rescue 和 SAM 后处理分配更重，更容易暴露 system allocator
   造成的跨核组访存与等待。

原始运行日志保存在四个配置子目录中。cross/system 的 `ERR1203383_PE.log`
还包含两次修复 cross 计数地址之前的零值诊断；应以文件最后一份非零报告为准。
