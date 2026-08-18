# SRR2496709 完整数据验证（空区间剪枝修复前）

> 本目录保留 2026-08-17 的历史结果。随后确认从核版
> `bwt_seed_strategy1()` 的空 FM-index 区间提前返回会跳过阈值前的
> 模糊碱基重启点，使完整 PE 输出比原版 BWA 多 2 字节。修复后的最终结果见
> `../bwt_zero_prune_mpi_full_20260818/README.md`。

## 输入与参数

- 参考基因组：`../data/GRCh38.d1.vd1.fa`
- PE：`../data/SRR2496709_1.fastq` 和 `../data/SRR2496709_2.fastq`
- SE：`../data/SRR2496709_1.fastq`
- 公共参数：`mem -v 4 -t 1 -I 170,80,500,1`
- 本轮运行均未使用 `-1`。

## 结果

| 配置 | 模式 | SAM 记录数 | 精确字节数 | MD5 |
| --- | --- | ---: | ---: | --- |
| MPI dynamic + single unordered，排序后 | PE | 25,396,635 | 8,119,600,202 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| MPI dynamic + single unordered，排序后 | SE | 12,684,037 | 3,685,790,019 | `322fe0c60d0ca365389307c06a21d77a` |
| single + system | PE | 25,396,635 | 8,119,600,202 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| single + system | SE | 12,684,037 | 3,685,790,019 | `322fe0c60d0ca365389307c06a21d77a` |
| cgs_cross + pool | PE | 25,396,635 | 8,119,600,202 | `f90eb3fae1e001d710e3bdbce5d1ffbf` |
| cgs_cross + pool | SE | 12,684,037 | 3,685,790,019 | `322fe0c60d0ca365389307c06a21d77a` |

三种配置得到的有序 PE 输出逐字节一致，有序 SE 输出也逐字节一致。

## 远端保留的 SAM

- MPI 原始无序输出：`SRR2496709_FULL_PE.sam`、`SRR2496709_FULL_SE.sam`
- MPI 排序后输出：`SRR2496709_FULL_PE.sorted.sam`、`SRR2496709_FULL_SE.sorted.sam`
- single + system：`SRR2496709_FULL_PE.single_system.sam`、`SRR2496709_FULL_SE.single_system.sam`
- cgs_cross + pool：`SRR2496709_FULL_PE.cgs_cross_pool.sam`、`SRR2496709_FULL_SE.cgs_cross_pool.sam`

## 排序说明

原来的全内存排序器在装入约 300 万条完整 SAM 记录后耗尽主核可用堆。
本结果目录中的 `external_sam_sort.cpp` 仅用于本轮结果处理：每批使用
128 MiB 内存生成有序 run，再执行稳定的多路归并。PE 共生成 70 个 run，
SE 共生成 32 个 run；成功结束后临时 run 已全部删除。相同 QNAME 的记录
保持原始相对顺序。

在当前神威批处理环境中，直接运行 `md5sum` 或 `wc` 会在打印完整结果后被
作业系统记录为退出码 1。上表仅采用对应 SWBWA/排序作业成功完成后的结果；
所有对比文件的精确字节数相同，MD5 也相同，因此可以确定内容逐字节一致。
