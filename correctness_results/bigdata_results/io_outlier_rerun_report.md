# Big-data I/O Outlier Rerun Report

## Scope

- Result set: five non-MPI configurations, `has1` and `no1`, three datasets,
  PE and SE.
- Correctness: all 60 cases pass their expected header-free SAM MD5.
- I/O metrics are taken from each run's `SWBWA Timing Report`:
  - input: `effective fread bandwidth`;
  - output: `POSIX write effective bandwidth`.

## Screening Result

The original result set contained one clear I/O outlier:

| Configuration | I/O mode | Case | Read bandwidth | Slowest fread |
| --- | --- | --- | ---: | ---: |
| `cgs_system` | `has1` | `ERR1203383 PE` | 28.072 MiB/s | 28.138 s |

All other original input measurements were at least 44.509 MiB/s. Output
measurements were 36.040-72.964 MiB/s. Values below 40 MiB/s occurred only in
the `no1` pipeline, where FASTQ input and SAM output overlap, and remained in
the expected tens-of-MiB/s range. They were therefore not treated as clear
filesystem stalls.

## Valid Rerun

- Scheduler job: `7870790`
- Queue: `q_share`
- Nodes: one
- Final status: `DONE`
- SAM MD5: `dc5c0a6babd41641db22808caedb7a44` (`PASS`)

| Metric | Original | Rerun |
| --- | ---: | ---: |
| Pipeline total | 503.262 s | 317.725 s |
| Stage 1 | 237.780 s | 59.192 s |
| Stage 2 | 121.190 s | 120.710 s |
| Stage 3 | 143.720 s | 137.229 s |
| FASTQ fread calls | 237.762 s | 59.187 s |
| Effective fread bandwidth | 28.072 MiB/s | 112.768 MiB/s |
| Slowest fread | 28.138 s | 1.694 s |
| POSIX write bandwidth | 68.404 MiB/s | 72.606 MiB/s |

The alignment time was nearly unchanged, while Stage 1 returned to the normal
range. This identifies the original slowdown as a transient input-filesystem
event rather than an SWBWA compute regression.

## Final I/O Range

After replacing the outlier log with the valid rerun:

| Metric | Minimum | Median | Maximum |
| --- | ---: | ---: | ---: |
| FASTQ read bandwidth | 44.509 MiB/s | 89.606 MiB/s | 139.187 MiB/s |
| SAM write bandwidth | 36.040 MiB/s | 58.951 MiB/s | 72.964 MiB/s |

No remaining case falls below 30 MiB/s. A second rerun was unnecessary because
the first valid rerun was normal and retained the expected MD5.

## Retained Files

- `cgs_system/has1/ERR1203383_PE.run.log`: valid rerun log used by the timing
  report.
- `cgs_system/has1/ERR1203383_PE.before_io_rerun.run.log`: original outlier log.
- `timing_summary.tsv` and `timing_report.html`: updated with the valid rerun.

Two earlier attempts (`7870735` and `7870757`) overlapped due to an orchestration
error and wrote the same output path. Both were terminated, their partial SAM
was removed, and neither attempt is included in this report.
