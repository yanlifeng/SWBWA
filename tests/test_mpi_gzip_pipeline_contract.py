#!/usr/bin/env python3
"""Host-only fastmap syntax and pipeline contracts, without MPI or CPE hardware.

Compile the real fastmap translation unit and reader ABI with SDK/MPI stubs.
Run its actual process callback with the actual single/threaded kthread driver;
extract only main_mem's post-pipeline teardown for ordered lifecycle checks.
The decoder, alignment and output are controlled stubs, not a distributed
deadlock, codec correctness, PE layout or native resource-placement test.
All generated sources and executables live in a temporary directory.
"""
import argparse
import itertools
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

ATHREAD_H = r"""
#pragma once
#define SLAVE_FUN(x) slave_##x
#define __uncached
void athread_init(void); void athread_init_cgs(void);
void athread_join(void); void athread_join_cgs(void);
void __real_athread_spawn(void *, void *, int);
void __real_athread_spawn_cgs(void *, void *, int);
"""
MPI_H = r"""
#pragma once
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#define MPI_SUCCESS 0
#define MPI_INT 1
#define MPI_MIN 1
#define MPI_MAX 2
int MPI_Allreduce(const void *, void *, int, int, int, MPI_Comm);
"""

COMMON_C = r"""
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static const char *scenario;
static void *batch_calloc(size_t count, size_t size) {
    if (scenario && !strcmp(scenario, "allocate")) return NULL;
    return calloc(count, size);
}
#define calloc batch_calloc
#include "fastmap.c"
#undef calloc

int bwa_verbose = 0;
void err_fatal(const char *header, const char *fmt, ...) {
    va_list args;
    fprintf(stderr, "[LOCAL_FATAL][%s] ", header);
    va_start(args, fmt); vfprintf(stderr, fmt, args); va_end(args);
    fputc('\n', stderr);
    _Exit(87);
}
void swbwa_mpi_abort(const char *message) {
    fprintf(stderr, "[MPI_ABORT] %s\n", message ? message : "");
    _Exit(86);
}
"""

FATAL_C = r"""
int main(int argc, char **argv) {
    assert(argc == 3);
#if SWBWA_ENABLE_GZIP_INPUT
    if (!strcmp(argv[1], "magic")) {
        is_gzip_input(argv[2]);
        return 99;
    }
#endif
    swbwa_pipeline_fatal("probe", "detail=%s code=%d tail", argv[2], 17);
    return 99;
}
"""

PIPELINE_C = r"""
struct swbwa_mpi_gzip_reader { int unused; };
static int paired, next_calls, align_calls, write_calls, flush_calls, phase;
static int chunk_count_calls, finish_calls;
static int is_case(const char *name) { return !strcmp(scenario, name); }
static int records(int round) {
    static const int counts[] = {0, 3, 0, 2, 0};
    return is_case("empty") ? 0 : counts[round];
}
static int64_t first_record(int round) {
    if (is_case("boundary") && records(round))
        return INT64_MAX / (paired ? 2 : 1) - records(round);
    return 100 + round * 7;
}
static int expected_calls(void) { return is_case("empty") ? 0 : 2; }
/* The last nonempty local ID is 10, not the final READY round's ID 13. */
static int64_t expected_chunks(void) { return is_case("empty") ? 0 : 11; }
static void unexpected(void) {
    fprintf(stderr, "unexpected raw/local-reader path\n");
    abort();
}
int swbwa_mpi_rank(void) { return 1; }
int swbwa_mpi_size(void) { return 3; }
int swbwa_mpi_is_root(void) { return 0; }
int swbwa_mpi_gzip_next(swbwa_mpi_gzip_reader *reader, swbwa_mpi_gzip_batch *batch) {
    assert(reader && phase == 0);
    int round = next_calls++;
    if (is_case("read")) return -1;
    memset(batch, 0, sizeof(*batch));
    if (round == 5) return 0;
    assert(round < 5);
    batch->records = records(round);
    batch->first_record = first_record(round);
    batch->chunk_id = round * 3 + 1;
    if (is_case("overflow")) {
        batch->first_record = INT64_MAX / (paired ? 2 : 1);
        batch->records = 1;
        return 1;
    }
    if (is_case("negative")) { batch->first_record = -1; return 1; }
    if (is_case("capacity")) {
        batch->size[0] = (int64_t)SWBWA_CPE_FORMAT_BUFFER_BYTES + 1;
        return 1;
    }
    /* Empty READY may still transfer independently owned empty allocations. */
    if (batch->records) {
        batch->partitions = calloc(SWBWA_CPE_COUNT, sizeof(*batch->partitions));
        assert(batch->partitions);
        batch->partition_count = SWBWA_CPE_COUNT;
        batch->partitions[0].records = batch->records;
    }
    for (int mate = 0; mate < (paired ? 2 : 1); ++mate) {
        batch->data[mate] = malloc(2);
        assert(batch->data[mate]);
        batch->data[mate][0] = (char)round;
        batch->data[mate][1] = 0;
        batch->size[mate] = batch->records ? 1 : 0;
    }
    return 1;
}
const char *swbwa_mpi_gzip_error(const swbwa_mpi_gzip_reader *reader) {
    (void)reader; return "injected reader detail";
}
int64_t swbwa_mpi_gzip_chunk_count(const swbwa_mpi_gzip_reader *reader) {
    assert(reader && phase == 0);
    assert(next_calls == 6 && flush_calls == 1);
    assert(align_calls == expected_calls() && write_calls == expected_calls());
    assert(chunk_count_calls++ == 0 && finish_calls == 0);
    return expected_chunks();
}
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
int swbwa_output_ordered_finish(int64_t chunks) {
    assert(phase == 0 && next_calls == 6 && flush_calls == 1);
    assert(chunk_count_calls == 1 && finish_calls++ == 0);
    assert(chunks == expected_chunks());
    assert(align_calls == expected_calls() && write_calls == expected_calls());
    errno = EIO;
    return is_case("finish") ? -1 : 0;
}
#endif
int swbwa_mpi_gzip_report(swbwa_mpi_gzip_reader *reader, FILE *stream) {
    assert(reader && stream && phase == 0);
    assert(next_calls == 6 && flush_calls == 1);
    assert(align_calls == expected_calls() && write_calls == expected_calls());
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    assert(chunk_count_calls == 1 && finish_calls == 1);
#else
    assert(chunk_count_calls == 0 && finish_calls == 0);
#endif
    phase = 1;
    return is_case("report") ? -1 : 0;
}
void swbwa_mpi_gzip_close(swbwa_mpi_gzip_reader *reader) {
    assert(reader && phase == 1); phase = 2;
}
void swbwa_host_workers_destroy(void) { assert(phase == 2); phase = 3; }
int swbwa_output_close(void) {
    assert(phase == 3); phase = 4;
    errno = EIO;
    return is_case("close") ? -1 : 0;
}
int swbwa_mpi_progress_thread_stop(void) {
    assert(phase == 4); phase = 5;
    errno = EIO;
    return is_case("stop") ? -1 : 0;
}
int swbwa_output_flush(void) {
    ++flush_calls;
    errno = EIO;
    return is_case("flush") ? -1 : 0;
}
void mem_process_seqs_merge2(const mem_opt_t *opt, const bwt_t *bwt,
        const bntseq_t *bns, const uint8_t *pac, int64_t processed, int *n,
        bseq1_t **seqs, char *r1, char *r2, long long size1, long long size2,
        const mem_pestat_t *pes, int **lengths, const swbwa_fastq_partition_t *partitions) {
    (void)opt; (void)bwt; (void)bns; (void)pac; (void)pes;
    int round = (unsigned char)r1[0];
    assert(round < 5 && records(round) > 0);
    assert(size1 == 1 && size2 == paired);
    assert(partitions && partitions[0].records == records(round));
    assert(processed == first_record(round) * (paired ? 2 : 1));
    assert(!paired || (r2 && r2[0] == r1[0]));
    usleep(2000);
    *n = records(round) * (paired ? 2 : 1);
    if (is_case("format")) --*n;
    *seqs = calloc((size_t)*n, sizeof(**seqs));
    assert(*seqs);
    if (lengths) {
        *lengths = calloc((size_t)*n, sizeof(**lengths));
        assert(*lengths);
    }
    ++align_calls;
    free(r1); free(r2);
}
int swbwa_output_write_chunk(int64_t id, int64_t start, int64_t end,
        const bseq1_t *seqs, const int *lengths, int reads, int is_paired) {
    int round = (int)((id - 1) / 3);
    assert(round == (write_calls == 0 ? 1 : 3));
    assert(reads == records(round) * (paired ? 2 : 1));
    assert(start == first_record(round) && end == start + records(round));
    assert(seqs && is_paired == paired);
#if SWBWA_OUTPUT_SINGLE_FILE
    assert(lengths);
#else
    assert(!lengths);
#endif
    ++write_calls;
    errno = EIO;
    return is_case("write") ? -1 : 0;
}
int swbwa_gzip_read_next(swbwa_gzip_reader *r, char **a, int64_t *b,
                         char **c, int64_t *d) {
    (void)r; (void)a; (void)b; (void)c; (void)d; unexpected(); return -1;
}
const char *swbwa_gzip_error(const swbwa_gzip_reader *r) {
    (void)r; unexpected(); return "unexpected";
}
int swbwa_gzip_check_pe_layout(swbwa_gzip_reader *r, const char *a, int64_t b,
                               const char *c, int64_t d, unsigned workers) {
    (void)r; (void)a; (void)b; (void)c; (void)d; (void)workers;
    unexpected(); return -1;
}
int swbwa_gzip_report(const swbwa_gzip_reader *r, FILE *f) {
    (void)r; (void)f; unexpected(); return -1;
}
void swbwa_gzip_close(swbwa_gzip_reader *r) { (void)r; unexpected(); }
void swbwa_input_usage_sample(swbwa_input_usage_t *u) { (void)u; unexpected(); }
void swbwa_input_usage_add(swbwa_input_usage_t *a, const swbwa_input_usage_t *b,
                           const swbwa_input_usage_t *c) {
    (void)a; (void)b; (void)c; unexpected();
}
int swbwa_input_read(const swbwa_input_slice_t *a, int b, int c, int d,
                      swbwa_input_result_t *e) {
    (void)a; (void)b; (void)c; (void)d; (void)e; unexpected(); return -1;
}
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
int swbwa_mpi_fastq_scheduler_next(swbwa_fastq_range_t *r) {
    (void)r; unexpected(); return -1;
}
void swbwa_mpi_fastq_scheduler_record_stage2(int64_t id, int64_t n, double s) {
    (void)id; (void)n; (void)s; unexpected();
}
void swbwa_mpi_fastq_scheduler_stats(int64_t *a, int64_t *b, int64_t *c) {
    (void)a; (void)b; (void)c; unexpected();
}
void swbwa_mpi_fastq_scheduler_close(void) { unexpected(); }
#endif

static void actual_teardown(ktp_aux_t aux) {
    int gzip_input = 1;
    double t0 = GetTime();
#include "actual_teardown.inc"
}
int main(int argc, char **argv) {
    assert(argc == 4);
    scenario = argv[1]; paired = atoi(argv[2]);
    struct swbwa_mpi_gzip_reader reader = {0};
    mem_opt_t opt = {0}; bwaidx_t idx = {0}; ktp_aux_t aux = {0};
    opt.flag = paired ? MEM_F_PE : 0;
    aux.opt = &opt; aux.idx = &idx; aux.is_paired = paired;
    aux.mpi_gzip_reader = &reader;
    aux.n_processed = 777; /* Must never replace the batch's exact ordinal. */
    if (atoi(argv[3])) kt_pipeline_queue(3, process, &aux, 3);
    else kt_pipeline_single(1, process, &aux, 3);
    assert(stage1_timing.completed_batches == (uint64_t)expected_calls());
    actual_teardown(aux);
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC || SWBWA_OUTPUT_SINGLE_FILE
    assert(phase == 5);
#else
    assert(phase == 4);
#endif
    puts("pipeline contract PASS");
    return 0;
}
"""


def run(command, *, expected=0, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, timeout=60, **kwargs)
    if result.returncode != expected:
        raise AssertionError(f"{shlex.join(map(str, command))}\n"
                             f"exit={result.returncode}, expected={expected}\n"
                             f"{result.stdout}\n{result.stderr}")
    return result


def function_body(source, name):
    """Match one definition, counting braces with comments/literals masked."""
    masked = re.sub(r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                    lambda match: ' ' * len(match[0]), source, flags=re.S)
    matches = list(re.finditer(r'\b' + re.escape(name) + r'\s*\([^;{}]*\)\s*\{', masked))
    if len(matches) != 1:
        raise AssertionError(f"expected one definition of {name}, found {len(matches)}")
    start = matches[0].end()
    depth = 1
    for pos in range(start, len(masked)):
        depth += (masked[pos] == '{') - (masked[pos] == '}')
        if depth == 0:
            return source[start:pos]
    raise AssertionError(f"unterminated function {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize', action='store_true', help='ASan/UBSan for runtime harnesses')
    args = parser.parse_args()
    cc = shlex.split(os.environ.get('CC', 'cc'))
    cxx = shlex.split(os.environ.get('CXX', 'c++'))
    source = (ROOT / 'src/host/fastmap.c').read_text()
    body = function_body(source, 'main_mem')
    for name in ('process', 'main_mem', 'is_gzip_input', 'init_cpe_allocator',
                 'read_fastq_block', 'read_fastq_parallel', 'open_fastq_input'):
        assert not re.search(r'\berr_(?:fatal|fflush)\s*\(', function_body(source, name)), name
    # Keep extraction inside the function; never capture a neighboring helper.
    start = body.index('/* Stage 3 has flushed;')
    end = body.index('swbwa_mpi_progress_thread_report()', start)
    teardown = body[start:end]
    assert body.index('kt_pipeline_queue(') < start
    order = ['swbwa_output_ordered_finish(', 'swbwa_mpi_gzip_report(', 'swbwa_mpi_gzip_close(',
             'swbwa_host_workers_destroy(', 'swbwa_output_close(',
             'swbwa_mpi_progress_thread_stop(']
    positions = [teardown.index(name) for name in order]
    assert positions == sorted(positions)
    assert 'online MPI gzip does not yet support single_ordered' not in body
    assert body.index('swbwa_output_open_online(') < body.index('kt_pipeline_queue(')
    assert teardown.index('swbwa_mpi_gzip_chunk_count(') < teardown.index('swbwa_mpi_gzip_report(')

    with tempfile.TemporaryDirectory(prefix='swbwa-mpi-gzip-contract-') as temporary:
        work = Path(temporary)
        (work / 'athread.h').write_text(ATHREAD_H)
        (work / 'mpi.h').write_text(MPI_H)
        (work / 'actual_teardown.inc').write_text(teardown)
        (work / 'fatal.c').write_text(COMMON_C + FATAL_C)
        (work / 'pipeline.c').write_text(COMMON_C + PIPELINE_C)
        disabled = work / 'disabled'
        disabled.mkdir()
        for name in ('swbwa_gzip.h', 'swbwa_mpi_gzip.h', 'swbwa_codec_runtime.h', 'codec_job.h'):
            (disabled / name).write_text('#error gzip dependency reached with GZIP_INPUT=0\n')
        includes = ['-I' + str(work), '-I' + str(ROOT / 'include'), '-I' + str(ROOT / 'src/host')]
        common = ['-D_GNU_SOURCE', '-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0',
                  '-DSWBWA_EXEC_MODE=SWBWA_EXEC_CGS_CROSS',
                  '-DSWBWA_CPE_ALLOC_MODE=SWBWA_CPE_ALLOC_POOL',
                  '-Wall', '-Wextra', '-Werror=implicit-function-declaration',
                  '-Werror=incompatible-pointer-types', '-Werror=format',
                  '-Wno-unused-parameter', '-Wno-unused-variable', '-Wno-unused-function']
        sanitizer = (['-fsanitize=address,undefined', '-fno-sanitize-recover=all',
                      '-fno-omit-frame-pointer', '-fno-pie'] if args.sanitize else [])

        def defines(mpi, gzip, mode, output):
            return [f'-DSWBWA_USE_MPI={mpi}', f'-DSWBWA_ENABLE_GZIP_INPUT={gzip}',
                    f'-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_{mode}',
                    f'-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_{output}']

        syntax_count = 0
        for mpi, gzip, mode, output in itertools.product(
                (0, 1), (0, 1), ('STATIC', 'DYNAMIC'),
                ('SPLIT', 'SINGLE_UNORDERED', 'DISCARD', 'SINGLE_ORDERED')):
            if not mpi and (mode != 'STATIC' or output.startswith('SINGLE_')):
                continue
            if output == 'SINGLE_ORDERED' and mode != 'DYNAMIC' and not gzip:
                continue
            flags = defines(mpi, gzip, mode, output)
            poison = ['-I' + str(disabled)] if not gzip else []
            run(cc + ['-std=gnu11', *common, *poison, *includes, *flags,
                      '-fsyntax-only', str(ROOT / 'src/host/fastmap.c')])
            syntax_count += 1
            if mpi and gzip and output == 'SINGLE_ORDERED':
                run(cc + ['-std=gnu11', *common, *includes, *flags,
                          '-DSWBWA_OUTPUT_RMA_ONLY=1', '-fsyntax-only',
                          str(ROOT / 'src/host/fastmap.c')])
                syntax_count += 1
        print(f'PASS {syntax_count} full fastmap syntax configurations (gzip on/off, MPI on/off)', flush=True)

        def build(name, flags, source_name, threaded=False):
            obj = work / (name + '.o')
            run(cc + ['-std=gnu11', '-O1', '-g', *common, *includes, *flags,
                      *sanitizer, '-ffunction-sections', '-fdata-sections', '-c',
                      str(work / source_name), '-o', str(obj)])
            objects = [str(obj)]
            if threaded:
                kt = work / (name + '-kthread.o')
                run(cxx + ['-std=c++11', '-O1', '-g', *includes, *flags, *sanitizer,
                           '-include', 'swbwa_config.h', '-ffunction-sections', '-pthread',
                           '-c', str(ROOT / 'src/host/kthread.cpp'), '-o', str(kt)])
                objects.append(str(kt))
            executable = work / name
            run(cxx + objects + sanitizer + (['-no-pie'] if args.sanitize else []) +
                ['-Wl,--gc-sections', '-pthread', '-o', str(executable)])
            return str(executable)

        long_detail = 'x' * 16000
        for mpi, gzip in itertools.product((0, 1), repeat=2):
            executable = build(f'fatal-{mpi}-{gzip}', defines(mpi, gzip, 'STATIC', 'SPLIT'), 'fatal.c')
            result = run([executable, 'fatal', long_detail], expected=86 if mpi else 87)
            assert f'detail={long_detail} code=17 tail' in result.stderr
            assert ('[MPI_ABORT]' in result.stderr) == bool(mpi)
            if gzip:
                missing = str(work / 'missing.fq.gz')
                result = run([executable, 'magic', missing], expected=86 if mpi else 87)
                assert missing in result.stderr
        print('PASS detailed fatal routing and rank-local magic-open failure', flush=True)

        error_messages = {
            'allocate': 'failed to allocate pipeline batch', 'read': 'injected reader detail',
            'overflow': 'MPI gzip batch exceeds read/index limits',
            'negative': 'MPI gzip batch exceeds read/index limits',
            'capacity': 'MPI gzip batch exceeds CPE formatting capacity',
            'format': 'MPI gzip formatter changed the decoded record count',
            'write': 'failed to write output chunk', 'flush': 'failed to flush SAM output',
            'report': 'failed to report MPI gzip input metrics',
            'close': 'failed to close SAM output', 'stop': 'MPI progress thread failed',
        }
        cases = 0
        for mode, output in itertools.product(
                ('STATIC', 'DYNAMIC'), ('SPLIT', 'SINGLE_UNORDERED', 'DISCARD', 'SINGLE_ORDERED')):
            executable = build(f'pipeline-{mode}-{output}', defines(1, 1, mode, output), 'pipeline.c', True)
            for paired, threaded, scenario in itertools.product((0, 1), (0, 1), ('mixed', 'empty', 'boundary')):
                result = run([executable, scenario, str(paired), str(threaded)],
                             env=dict(os.environ, ASAN_OPTIONS='detect_leaks=1:halt_on_error=1'))
                assert 'pipeline contract PASS' in result.stdout
                cases += 1
            errors = dict(error_messages)
            if output == 'SINGLE_ORDERED':
                errors['finish'] = 'failed to finalize gzip ordered chunk count'
            for threaded, (scenario, message) in itertools.product((0, 1), errors.items()):
                if scenario == 'stop' and mode == 'STATIC' and not output.startswith('SINGLE_'):
                    continue
                result = run([executable, scenario, '1', str(threaded)], expected=86,
                             env=dict(os.environ, ASAN_OPTIONS='detect_leaks=0:halt_on_error=1'))
                assert message in result.stderr, result.stderr
                assert '[MPI_ABORT]' in result.stderr and '[LOCAL_FATAL]' not in result.stderr
                cases += 1
            print(f'PASS {mode}/{output}: SE/PE, single/queue, empty READY, EOF, close, injected errors', flush=True)
        print(f'PASS {cases} pipeline cases; no MPI launcher, remote commands or native CPE execution', flush=True)


if __name__ == '__main__':
    main()
