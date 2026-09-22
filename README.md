# SWBWA

SWBWA is a high-accuracy, high-performance short-read aligner optimized for the next-generation Sunway platform and based on BWA-MEM.

## Features

- A parallel framework redesigned for Sunway's heterogeneous architecture, with software prefetching and memory-access optimizations for bigshare mode.
- Up to **330x** speedup over the unoptimized single-threaded version, and **1.2-1.4x** faster than [bwamem](https://github.com/lh3/bwa) on a dual-socket 48-core x86 server, with nearly identical results.

## Directory layout

- `src/host/`: MPE C/C++ source code, including the alignment pipeline, MPI input/output, and host-side helpers.
- `src/slave/`: CPE source code and headers used only by the slave side.
- `include/`: public headers, configuration headers, and the generated CPE layout header.
- `tools/`: scripts for extracting addresses and TLS information during cross-segment builds.
- `tests/`: MPI/RMA and runtime diagnostic programs excluded from the default build.
- `scripts/`: correctness checks, performance tests, and result-analysis scripts.
- `correctness_results/`: run logs, correctness results, and experiment notes.

## Build

SWBWA supports only the next-generation Sunway platform.

### Dependencies

- `sw9gcc` (7.1.0 or newer)
- `zlib`

### Compilation

```bash
git clone https://github.com/RabbitBio/SWBWA.git
cd SWBWA
make -j4
```

### Build configuration

The default build uses single-CG execution, CPE FASTQ formatting, the system CPE allocator, dynamic MPI input, and unordered single-file MPI output:

```bash
make print-config
```

FASTQ formatting is always performed on the CPEs; there is no separate formatting-mode switch.

The supported build variables are:

| Variable | Values | Default |
| --- | --- | --- |
| `EXEC_MODE` | `single`, `cgs`, `cgs_cross` | `single` |
| `CPE_ALLOCATOR` | `system`, `pool` | `system` |
| `HOST_MALLOC_WRAPPER` | `0`, `1` | `1` |
| `HOST_MALLOC_STATS` | `0`, `1` | `0` |
| `CPE_KERNEL_OPT` | `0`, `1` | `1` for non-MPI `cgs_cross + pool`; `0` otherwise |
| `CPE_LDM_MODE` | `0` off, `1` tiered malloc pool, `2` manual | `2`; mode 1 requires non-MPI `cgs_cross + pool` |
| `CPE_DISCARD_DIGEST` | `0`, `1` | `1` for non-MPI `cgs_cross + pool` with `OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`; `0` otherwise |
| `USE_MPI` | `0`, `1` | `1` |
| `MPI_INPUT_MODE` | `static`, `dynamic` | `dynamic` with MPI |
| `OUTPUT_MODE` | `split`, `single_unordered`, `discard` | `single_unordered` with MPI; `split` otherwise |
| `MPI_EXACT_READ_INDEX` | `0`, `1` | `1` |

`MPI_INPUT_MODE` and `MPI_EXACT_READ_INDEX` are only effective when `USE_MPI=1`. `OUTPUT_MODE=discard` is also supported without MPI; non-MPI `split` retains the ordinary single-file writer. `single_unordered` requires MPI.

`CPE_LDM_MODE` controls explicit private-object placement, not runtime cache,
stack, SIMD or algorithm choices. Mode 0 disables both automatic and manual
LDM allocation. Mode 1 disables manual placement and routes wrapped allocation
through the existing 32 KiB lifetime/access-tier pool, with cross-segment heap
fallback. SAM and public/cached query profiles stay on the heap. Mode 2 disables
the new pool and retains the validated manual workspace/phase-reuse paths.
Pool metadata counts against the shared 40 KiB budget. The old Makefile options
`CPE_LDM_ALLOC`, `CPE_MANUAL_LDM`, `CPE_LDM_BYTES` and `LDM_SCRATCH_BUDGET`
now fail explicitly instead of silently mixing policies.
See [LDM configuration](docs/LDM_ALLOCATOR.md) for ownership and safety limits.

```bash
# Select 0/1/2; cross mode always requires the complete two-pass build.
bash build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 CPE_LDM_MODE=2
```

`MPI_EXACT_READ_INDEX=1` is intended for correctness checks. Rank 0 scans the complete FASTQ once before alignment to build exact record prefixes.

Dynamic MPI input uses the configured large chunks for the first 90% of the FASTQ. In the final 10% it uses chunks one quarter that size, and the final two rank waves use chunks one quarter of the medium-tail size. Set `SWBWA_MPI_TAIL_PERCENT=0` to disable tail refinement, or set `SWBWA_MPI_FINE_TAIL_WAVES=0` to keep the 10% medium tail without the final fine region.

`OUTPUT_MODE=discard` is a profiling mode: it creates no SAM file. By default it hashes each non-empty SAM blob submitted through the output interface and prints per-rank order-independent sum and XOR fingerprints. Set `SWBWA_DISCARD_HASH=0` only when measuring hash overhead itself.

`CPE_DISCARD_DIGEST` defaults to `1` only when all of
`EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`
are selected; it defaults to `0` otherwise. Explicitly set `CPE_DISCARD_DIGEST=0`
for core-only stage 2 timing comparisons and final-copy validation: the MPE
hashes the copied final SAM buffer in stage 3.
The eligible default, `CPE_DISCARD_DIGEST=1`, moves hashing to the CPE copy worker
in stage 2 part 5 and replaces the final SAM copy with per-read digest metadata
for overall discard throughput.
Keep part 3, part 5, and stage 3 separate when reporting this change.

The enabled path requires `EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0`
and `OUTPUT_MODE=discard DISCARD_HASH_BYTES=0`. It fingerprints every **generated SAM blob byte**
except the terminating NUL, preserving per-read/mate boundaries and
`calls/bytes/sum/xor` semantics. It does **not** execute or validate the skipped
final SAM copy; the log identifies `hash_scope=generated_sam final_sam_copy=0`.
Normal SAM output modes are unchanged. `SWBWA_DISCARD_HASH=0` also disables
CPE hashing, retaining calls/bytes but not FULL correctness evidence.
Set the Make variable `CPE_DISCARD_DIGEST`, not a duplicate definition in
`EXTRA_CPPFLAGS`, and repeat the complete two-pass build after changing it.

Core-only example for single-process, six-CG stage 2 measurements and FULL
final-copy validation (explicit `CPE_DISCARD_DIGEST=0`):

```bash
./build.sh 6 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 \
    OUTPUT_MODE=discard DISCARD_HASH_BYTES=0 CPE_DISCARD_DIGEST=0 CPE_PROFILE=0
SWBWA_DISCARD_HASH=1 bsub -I -b -q q_share -n 1 -cgsp 64 -mpecg 6 \
    -share_size 2000 -xmalloc -cross_size 42000 -cache_size 128 -priv_size 16 \
    ./SWBWA mem -v 4 -t 1 -1 -I 170,80,500,1 -o ignored.sam \
    ref.fa read1.fq read2.fq
```

`-1` runs the three stages serially. Compare stage 2 and its CPE part 3
separately from input and discard-hashing time. Validate all four output
fields (`calls`, `bytes`, `sum`, `xor`) against the same-input baseline, with
`enabled=1 hash_prefix_bytes=0`; a prefix-only or disabled hash is not a
full-output correctness check. These fingerprints are order-independent
regression checks, not a proof of equality or a replacement for SAM MD5.

For the default discard-throughput configuration, rebuild with the following
command and use the same run command above. Both `CPE_KERNEL_OPT` and
`CPE_DISCARD_DIGEST` default to `1` here; FULL checking covers generated SAM
(`hash_scope=generated_sam`), not the skipped final copy.

```bash
./build.sh 6 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 \
    OUTPUT_MODE=discard DISCARD_HASH_BYTES=0 CPE_PROFILE=0
```

### Exact CPE Kernel Optimizations

`CPE_KERNEL_OPT=1` groups the following changes without changing SIMD lane
counts, scoring, lazy-F termination, tie comparisons, or read scheduling:

- Reuse the completed SMEM collection's 16 KiB LDM scratch for chain seeds.
- Place extension EH/QP scratch in LDM when its combined size is at most
  4096 bytes; retain the existing reusable heap scratch as fallback.
- Build the five-symbol extension query profile with shared query loads.
- Fuse integer SIMD gap clamps and use predicate-based XOR selection.
- Reuse MPE SAM pointer/length scratch and leave terminators to the CPE copy.

The existing 40 KiB tracked LDM budget is unchanged. This switch does not
enable experimental FP16 lane layouts. Use `CPE_KERNEL_OPT=0` for the kernel
control build and always repeat the complete `build.sh` workflow for cross
execution. The single-process discard writer and batch-level output timing
are independent of this switch.

Portable differential checks use compiler-vector adapters and sanitizers;
they supplement, but do not replace, full Sunway output validation:

```bash
bash tests/run_host_checks.sh
bash tests/run_cpe_kernel_checks.sh
```

For example:

```bash
make clean
make -j4 EXEC_MODE=cgs CPE_ALLOCATOR=system
```

Run `make clean` before changing build modes because Make does not track compiler-flag changes in existing object files.

After source changes, regenerate the linked CPE layout and relocation data for cross-segment execution:

```bash
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=system
```

The first argument to `build.sh` is the number of parallel compilation jobs; subsequent arguments are passed directly as Make variables. The script performs the placeholder build, extracts ELF segment addresses, generates the layout header, performs the second build, and generates cross-segment relocation data.

Common configuration examples:

```bash
./build.sh 8 EXEC_MODE=single CPE_ALLOCATOR=system USE_MPI=0
./build.sh 8 EXEC_MODE=cgs CPE_ALLOCATOR=pool USE_MPI=0
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0
./build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=1 \
    MPI_INPUT_MODE=dynamic OUTPUT_MODE=single_unordered
```

## Basic usage

### Index the reference

SWBWA is compatible with index files generated by [bwamem](https://github.com/lh3/bwa).

```bash
./SWBWA index ref.fa
```

### Align reads

Single-end reads:

```bash
./SWBWA mem ref.fa reads.fq -o aln.sam
```

Paired-end reads:

```bash
./SWBWA mem ref.fa read1.fq read2.fq -o aln.sam
```

In an MPI job, each rank reads its own FASTQ range; rank 0 does not load the entire input into memory:

```bash
bsub -I -b -q q_share -N 1 -np 6 -cgsp 64 \
  -share_size 12000 -cache_size 128 -priv_size 16 \
  ./SWBWA mem -t 1 -1 -K 5000000 -o out.sam \
  ref.fa read1.fq read2.fq
```

`OUTPUT_MODE=split` creates one output file per rank. `single_unordered` uses MPI RMA to atomically reserve ranges in one output file, but does not guarantee record order. `discard` writes no SAM and is useful for measuring stage 2 performance.

The unified run entry point only submits jobs; it does not rebuild the program:

```bash
./run.sh single -- ./SWBWA mem -t 1 -o out.sam ref.fa reads.fq
./run.sh cgs_cross -- ./SWBWA mem -t 1 -o out.sam ref.fa reads.fq
./run.sh mpi --nodes 1 --ranks 6 -- \
  ./SWBWA mem -t 1 -K 5000000 -o out.sam ref.fa reads.fq
```

See [`scripts/README.md`](scripts/README.md) and [`correctness_results/README.md`](correctness_results/README.md) for testing scripts and result-directory conventions.

## Getting help

```bash
./SWBWA mem
```

This command displays the complete parameter reference.
