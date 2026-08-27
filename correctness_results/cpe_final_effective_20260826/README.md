# Final effective CPE optimizations

This directory records the final regression after removing the ineffective
BWT experiments. The retained production changes are:

- a per-read LDM arena for small MEM-chain seed arrays;
- the FP16-vector KSW backend with the original 16 logical lanes preserved;
- detailed LWPF regions used to verify CPE hotspots.

The SA batch, RID table copy, and BWT prefetch-trim experiments were removed
because their isolated gains did not survive the interaction test with the
chain arena.

## Configuration

```text
EXEC_MODE=cgs_cross
CPE_ALLOCATOR=pool
CPE_PROFILE=1
USE_MPI=0
FLOAT16_VECTOR=1
-t 1 -1 -v 4 -I 170,80,500,1
```

The FP16 backend is enabled by default and can be disabled with
`make FLOAT16_VECTOR=0 ...` for ablation or fallback builds.

## Full-output regression

| Dataset | Reads | SAM bytes | MD5 | Result | Stage 2 | CPE part 3 |
| --- | ---: | ---: | --- | --- | ---: | ---: |
| ERR1203383 PE (75 bp) | 2,000,259 | 523,048,834 | `c2af4bf0b057d5125ce2a9d770f13741` | PASS | 3.874 s | 3.450 s |
| SRR7963242 PE (150 bp) | 5,033,412 | 2,107,655,733 | `a799dc7268f389120ec92820e04b5118` | PASS | 45.030 s | 44.068 s |

Both checks used the complete SAM stream. The temporary SAM files were deleted
after the line count, byte count, and MD5 checks passed; build, timing, and LWPF
logs remain in this directory.

The default `CPE_PROFILE=0` production build was compiled separately and passed
an additional full ERR1203383 PE smoke run with output directed to `/dev/null`.
Its timing report recorded `stage2=3.997 s` and `part3=3.670 s`; this single run
is a build-path check, not part of the profiling-on comparison above.
