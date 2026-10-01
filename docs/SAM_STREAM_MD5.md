# Ordered SAM validation without a file

Build ordinary non-MPI output, then set `SWBWA_OUTPUT_MD5=1` at runtime:

```bash
bash build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 OUTPUT_MODE=split
SWBWA_OUTPUT_MD5=1 ./run.sh cgs_cross -- ./SWBWA mem -v 4 -t 1 \
    -I 170,80,500,1 -o unused.sam ref.fa reads_1.fastq reads_2.fastq
```

The example deliberately omits `-1`; it exercises the overlapping pipeline.
Ensure the job launcher passes the environment variable to the application.
The log must name `(ordered SAM streaming MD5; no file)` and end with one
`SAM_STREAM_MD5` JSON record with `status=COMPLETE`. `unused.sam` is not created.
No completion record, an unsuccessful job exit, or a count mismatch is a failure.

This mode generates and copies the complete SAM as normal. The final buffered
writer consumes the same bytes into an incremental MD5 instead of writing them.
Only leading `@HD`, `@SQ`, `@RG`, `@PG`, `@CO` header lines are removed. Records
are not sorted, re-encoded or sampled. Compare `md5`, `bytes` and `records`
against an independent header-free ordered reference. Records count SAM lines,
not reads; supplementary alignments can make them differ.

The optional standalone consumer uses the same implementation:

```bash
c++ -std=c++11 -O2 -Iinclude -DSWBWA_SAM_MD5_STANDALONE \
    src/host/swbwa_sam_md5.cpp -o sam_stream_md5
set -o pipefail
./bwa mem -t 64 -I 170,80,500,1 ref.fa reads_1.fastq reads_2.fastq | ./sam_stream_md5
```

This is a regression checksum, not cryptographic authentication. It cannot be
combined with MPI, `OUTPUT_MODE=discard` or `OUTPUT_RMA_ONLY=1`. Invalid values
are rejected. Unset the variable (or set it to 0) for ordinary output.
Empty or truncated streams and embedded NUL bytes are rejected. MD5/counting
cost belongs to Stage3 and must not be reported as real filesystem throughput.
See `tests/test_sam_stream_md5.py` for boundary and >4 GiB tests.
