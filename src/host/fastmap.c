/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include "swbwa_config.h"
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "swbwa_cpe.h"
#include "swbwa_runtime.h"
#include "swbwa_mpi.h"
#include "swbwa_output.h"
#include "swbwa_host_workers.h"
#include "swbwa_input.h"
#include "swbwa_discard_digest.h"
#include "swbwa_cpe_profile.h"

#include <athread.h>
//KSEQ_DECLARE(gzFile)
KSEQ_DECLARE(int)

extern unsigned char nst_nt4_table[256];
extern void kt_pipeline_wait_report(void);

extern void SLAVE_FUN(state_init());

double t_malloc = 0;
double t_free = 0;

double t_tot = 0;
double t_step1 = 0;
double t_step1_1 = 0;
double t_step1_1_1 = 0;
double t_step1_read = 0;
double t_step2 = 0;
double t_step3 = 0;
double t_step3_1 = 0;

double t_work1 = 0;
double t_work2 = 0;

double t_work1_1 = 0;
double t_work1_2 = 0;
double t_work1_3 = 0;
double t_work1_4 = 0;
double t_work1_5 = 0;
double t_work1_6 = 0;


double t_work2_1 = 0;
double t_work2_2 = 0;
double t_work2_3 = 0;
double t_work2_4 = 0;
double t_work2_5 = 0;

double t_extend = 0;
double t_bwt_sa = 0;

typedef struct {
    uint64_t stage_calls;
    uint64_t completed_batches;
    uint64_t allocation_calls;
    uint64_t allocation_bytes;
    uint64_t seek_calls;
    uint64_t read_calls;
    uint64_t read_bytes;
    uint64_t boundary_calls;
    double allocation_seconds;
    double allocation_max_seconds;
    double seek_seconds;
    double seek_max_seconds;
    double read_seconds;
    double read_worker_seconds;
    double read_max_seconds;
    double boundary_seconds;
    double boundary_max_seconds;
    uint64_t slowest_allocation_bytes;
    int slowest_read_input;
    int64_t slowest_read_offset;
    uint64_t slowest_read_bytes;
    swbwa_input_usage_t usage;
} stage1_timing_t;

/* Stage 1 has exactly one producer, even in the threaded pipeline. */
static stage1_timing_t stage1_timing;
static int input_readers;

long long s_reg_sum = 0;
long long c_px2 = 0;
long long s_px2 = 0;

static void print_timing_line(const char *label, double seconds, int indent)
{
    fprintf(stderr, "%*s%-52s %10.3f s\n", indent, "", label, seconds);
}

static void print_timing_stat(const char *label, double seconds,
                              uint64_t calls, double max_seconds, int indent)
{
    fprintf(stderr,
            "%*s%-52s %10.3f s  (%" PRIu64 " calls, max %.3f s)\n",
            indent, "", label, seconds, calls, max_seconds);
}

static void record_stage1_allocation(uint64_t bytes, double seconds)
{
    ++stage1_timing.allocation_calls;
    stage1_timing.allocation_bytes += bytes;
    stage1_timing.allocation_seconds += seconds;
    if (seconds > stage1_timing.allocation_max_seconds) {
        stage1_timing.allocation_max_seconds = seconds;
        stage1_timing.slowest_allocation_bytes = bytes;
    }
    if (bwa_verbose >= 4 && seconds >= 1.0)
        fprintf(stderr,
                "[DBG stage1] slow allocation: bytes=%" PRIu64
                " elapsed=%.3f s\n",
                bytes, seconds);
}

static void print_timing_report_body(void)
{
    double stage1_accounted = stage1_timing.allocation_seconds +
                              stage1_timing.seek_seconds +
                              stage1_timing.read_seconds +
                              stage1_timing.boundary_seconds;
    double stage1_other = t_step1 - stage1_accounted;
    double stage1_read_mib = (double)stage1_timing.read_bytes / (1 << 20);
    double stage3_cleanup = t_step3 - t_step3_1;

    if (stage1_other < 0.0) stage1_other = 0.0;
    if (stage3_cleanup < 0.0) stage3_cleanup = 0.0;

    fprintf(stderr,
            "\n"
            "======================= SWBWA Timing Report =======================\n");
#if SWBWA_USE_MPI
    fprintf(stderr, "  MPI rank: %06d / %06d\n",
            swbwa_mpi_rank(), swbwa_mpi_size());
#endif
    fprintf(stderr,
            "  Accumulated wall-clock time; indented rows are included in\n"
            "  their parent row and must not be added to it.\n"
            "\n"
            "  Pipeline\n");
    print_timing_line("total - complete three-stage pipeline", t_tot, 4);
    print_timing_line("stage 1 - allocate and read raw FASTQ blocks", t_step1, 4);
    print_timing_line("stage 2 - align reads and generate SAM records", t_step2, 4);
    print_timing_line("stage 3 - write SAM and release batch data", t_step3, 4);
    print_timing_line("SAM output writes", t_step3_1, 6);
    print_timing_line("batch cleanup and loop overhead (derived)", stage3_cleanup, 6);

    fprintf(stderr, "\n  Stage 1 input details\n");
    fprintf(stderr,
            "    batches completed / stage calls              %10" PRIu64
            " / %" PRIu64 "\n",
            stage1_timing.completed_batches, stage1_timing.stage_calls);
    print_timing_stat("input-buffer and batch allocations",
                      stage1_timing.allocation_seconds,
                      stage1_timing.allocation_calls,
                      stage1_timing.allocation_max_seconds, 4);
    fprintf(stderr,
            "    aggregate requested allocation bytes           %10" PRIu64
            "\n",
            stage1_timing.allocation_bytes);
    print_timing_stat("FASTQ fseeko calls", stage1_timing.seek_seconds,
                      stage1_timing.seek_calls,
                      stage1_timing.seek_max_seconds, 4);
    print_timing_line("FASTQ read wall time", stage1_timing.read_seconds, 4);
    print_timing_stat(input_readers ? "pread calls (worker time sum)" : "FASTQ fread calls",
                      stage1_timing.read_worker_seconds,
                      stage1_timing.read_calls,
                      stage1_timing.read_max_seconds, 6);
    fprintf(stderr, "    read API / readers                            %s / %d\n",
            input_readers ? "pread" : "fread", input_readers ? input_readers : 1);
    if (bwa_verbose >= 4) {
        const swbwa_input_usage_t *u = &stage1_timing.usage;
        fprintf(stderr, "    thread resource samples                       %u\n", u->samples);
        if (u->samples) {
            print_timing_line("reader user CPU (worker sum)", u->user_seconds, 6);
            print_timing_line("reader system CPU (worker sum)", u->system_seconds, 6);
            fprintf(stderr,
                    "      reader minor / major faults                 %" PRId64 " / %" PRId64 "\n"
                    "      reader voluntary / involuntary switches     %" PRId64 " / %" PRId64 "\n",
                    u->minor_faults, u->major_faults,
                    u->voluntary_switches, u->involuntary_switches);
        } else {
            fprintf(stderr, "      thread resource accounting unavailable\n");
        }
    }
    print_timing_stat("align chunk ends to FASTQ record boundaries",
                      stage1_timing.boundary_seconds,
                      stage1_timing.boundary_calls,
                      stage1_timing.boundary_max_seconds, 4);
    print_timing_line("scheduler, checks and bookkeeping (derived)",
                      stage1_other, 4);
    fprintf(stderr,
            "    raw FASTQ bytes read                           %10" PRIu64
            "  (%.3f MiB)\n",
            stage1_timing.read_bytes, stage1_read_mib);
    if (stage1_timing.read_seconds > 0.0)
        fprintf(stderr,
                "    effective read bandwidth                      %10.3f MiB/s\n",
                stage1_read_mib / stage1_timing.read_seconds);
    if (stage1_timing.allocation_calls > 0)
        fprintf(stderr,
                "    slowest allocation                            %10.3f s"
                "  (%" PRIu64 " bytes)\n",
                stage1_timing.allocation_max_seconds,
                stage1_timing.slowest_allocation_bytes);
    if (stage1_timing.read_calls > 0)
        fprintf(stderr,
                "    slowest read call                             %10.3f s"
                "  (R%d offset=%" PRId64 ", bytes=%" PRIu64 ")\n",
                stage1_timing.read_max_seconds,
                stage1_timing.slowest_read_input,
                stage1_timing.slowest_read_offset,
                stage1_timing.slowest_read_bytes);

    fprintf(stderr, "\n  Stage 2 worker details\n");
    print_timing_line("total - active merge worker", t_work1, 4);
    print_timing_line("part 1 - prepare CPE task and reusable buffers", t_work1_1, 6);
    print_timing_line("part 2 - CPE FASTQ formatting and input release", t_work1_2, 6);
    print_timing_line("part 3 - CPE alignment and SAM length pass", t_work1_3, 6);
    print_timing_line("part 4 - assign slices in the shared SAM buffer", t_work1_4, 6);
    print_timing_line("part 5 - CPE SAM record generation", t_work1_5, 6);
    print_timing_line("part 6 - release temporary worker data", t_work1_6, 6);

    if (t_work2 != 0.0 || t_work2_1 != 0.0 || t_work2_2 != 0.0 ||
        t_work2_3 != 0.0 || t_work2_4 != 0.0 || t_work2_5 != 0.0) {
        fprintf(stderr, "\n  Legacy split-SAM worker details\n");
        print_timing_line("total - legacy SAM worker", t_work2, 4);
        print_timing_line("part 1 - CPE SAM length pass", t_work2_1, 6);
        print_timing_line("part 2 - allocate per-read SAM buffers", t_work2_2, 6);
        print_timing_line("part 3 - CPE SAM record generation", t_work2_3, 6);
        print_timing_line("part 4 - release alignment and SAM metadata", t_work2_4, 6);
        print_timing_line("part 5 - release remaining worker data", t_work2_5, 6);
    }

    fprintf(stderr, "\n  Reserved allocator timers (no active timing sites)\n");
    print_timing_line("malloc timer", t_malloc, 4);
    print_timing_line("free timer", t_free, 4);
    fprintf(stderr,
            "===================================================================\n"
            "\n");
}

static void print_timing_report(void)
{
    swbwa_mpi_print_rank_ordered(print_timing_report_body);
}

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_single(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_thread(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_queue(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);


typedef struct {
	mem_opt_t *opt;
	mem_pestat_t *pes;
	int64_t n_processed;
	int64_t fastq_bytes_per_cg;
	int is_paired;
	bwaidx_t *idx;
	int64_t fastq_chunk_bytes;
	int64_t input_position[2];
	int64_t input_end[2];
} ktp_aux_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	int *sam_lengths;
	int64_t chunk_id;
	int64_t n_processed;
	char *fastq_buffer[2];
	long long fastq_size[2];
	int64_t chunk_start, chunk_end;
} ktp_data_t;

static void skip_to_line_end(char *data_, long long *pos_, const long long size_) {
    while (*pos_ < size_ && data_[*pos_] != '\n') {
        ++(*pos_);
    }
}

static int64_t get_next_fastq(char *data_, long long pos_, const long long size_) {
    if(pos_ < 0) pos_ = 0;
    skip_to_line_end(data_, &pos_, size_);
    if (pos_ >= size_) return size_;
    ++pos_;

    while (pos_ < size_ && data_[pos_] != '@') {
        skip_to_line_end(data_, &pos_, size_);
        if (pos_ >= size_) return size_;
        ++pos_;
    }
    if (pos_ >= size_) return size_;
    int64_t pos0 = pos_;

    skip_to_line_end(data_, &pos_, size_);
    if (pos_ >= size_) return pos0;
    ++pos_;

    if (pos_ < size_ && data_[pos_] == '@')
        return pos_;
    skip_to_line_end(data_, &pos_, size_);
    if (pos_ >= size_) return pos0;
    ++pos_;
    if (pos_ >= size_ || data_[pos_] != '+')
        err_fatal(__func__, "invalid FASTQ record near byte %lld", pos0);
    return pos0;
}


static FILE *file1_ptr;
static FILE *file2_ptr;

static int64_t read_fastq_block(FILE *file, int64_t *file_offset,
                                int64_t file_end, char *buffer,
                                int64_t capacity)
{
    const int input_index = file == file2_ptr ? 2 : 1;
    const int64_t read_offset = *file_offset;
    int64_t bytes_to_read;
    int64_t bytes_read;
    int64_t record_bytes;
    int seek_result;
    int read_error;
    double start;
    double elapsed;
    swbwa_input_usage_t before = {0}, after;

    if (*file_offset >= file_end) return 0;
    bytes_to_read = capacity;
    if (file_end - *file_offset < bytes_to_read)
        bytes_to_read = file_end - *file_offset;
    start = GetTime();
    seek_result = fseeko(file, (off_t)*file_offset, SEEK_SET);
    elapsed = GetTime() - start;
    ++stage1_timing.seek_calls;
    stage1_timing.seek_seconds += elapsed;
    if (elapsed > stage1_timing.seek_max_seconds)
        stage1_timing.seek_max_seconds = elapsed;
    if (bwa_verbose >= 4 && elapsed >= 1.0)
        fprintf(stderr,
                "[DBG stage1] slow fseeko: R%d offset=%" PRId64
                " elapsed=%.3f s\n",
                input_index, read_offset, elapsed);
    if (seek_result != 0)
        err_fatal(__func__, "failed to seek FASTQ input: %s", strerror(errno));

    if (bwa_verbose >= 4) swbwa_input_usage_sample(&before);
    start = GetTime();
    bytes_read = (int64_t)fread(buffer, 1, (size_t)bytes_to_read, file);
    read_error = errno;
    elapsed = GetTime() - start;
    if (bwa_verbose >= 4) {
        swbwa_input_usage_sample(&after);
        swbwa_input_usage_add(&stage1_timing.usage, &before, &after);
    }
    ++stage1_timing.read_calls;
    stage1_timing.read_seconds += elapsed;
    stage1_timing.read_worker_seconds += elapsed;
    if (bytes_read > 0)
        stage1_timing.read_bytes += (uint64_t)bytes_read;
    if (elapsed > stage1_timing.read_max_seconds) {
        stage1_timing.read_max_seconds = elapsed;
        stage1_timing.slowest_read_input = input_index;
        stage1_timing.slowest_read_offset = read_offset;
        stage1_timing.slowest_read_bytes = bytes_read > 0
                                         ? (uint64_t)bytes_read : 0;
    }
    if (bwa_verbose >= 4 && elapsed >= 1.0)
        fprintf(stderr,
                "[DBG stage1] slow fread: R%d offset=%" PRId64
                " requested=%" PRId64 " read=%" PRId64
                " elapsed=%.3f s\n",
                input_index, read_offset, bytes_to_read, bytes_read, elapsed);
    if (ferror(file))
        err_fatal(__func__, "failed to read FASTQ input: %s", strerror(read_error));

    record_bytes = bytes_read;
    if (bytes_read == capacity && *file_offset + bytes_read < file_end) {
        start = GetTime();
        record_bytes = get_next_fastq(buffer, bytes_read - (1 << 10), bytes_read);
        elapsed = GetTime() - start;
        ++stage1_timing.boundary_calls;
        stage1_timing.boundary_seconds += elapsed;
        if (elapsed > stage1_timing.boundary_max_seconds)
            stage1_timing.boundary_max_seconds = elapsed;
        if (bwa_verbose >= 4 && elapsed >= 1.0)
            fprintf(stderr,
                    "[DBG stage1] slow FASTQ boundary alignment:"
                    " R%d offset=%" PRId64 " elapsed=%.3f s\n",
                    input_index, read_offset, elapsed);
    }
    *file_offset += record_bytes;
    return record_bytes;
}

static void read_fastq_parallel(ktp_aux_t *aux, char *buffer, char *buffer2,
                                int64_t capacity, long long *size,
                                long long *size2)
{
    swbwa_input_slice_t inputs[2];
    swbwa_input_result_t result;
    int count = aux->is_paired ? 2 : 1;
    int i;
    char *buffers[2] = {buffer, buffer2};
    FILE *files[2] = {file1_ptr, file2_ptr};
    long long *sizes[2] = {size, size2};

    for (i = 0; i < count; ++i) {
        int64_t remaining = aux->input_end[i] - aux->input_position[i];
        inputs[i].fd = fileno(files[i]);
        inputs[i].buffer = buffers[i];
        inputs[i].offset = aux->input_position[i];
        inputs[i].length = remaining <= 0 ? 0 :
            (size_t)(remaining < capacity ? remaining : capacity);
    }
    if (swbwa_input_read(inputs, count, input_readers, bwa_verbose >= 4, &result) != 0)
        err_fatal(__func__, "positioned FASTQ read failed: %s", strerror(errno));
    stage1_timing.read_calls += result.calls;
    stage1_timing.read_bytes += result.bytes;
    stage1_timing.read_seconds += result.wall_seconds;
    stage1_timing.read_worker_seconds += result.syscall_seconds;
    if (result.max_seconds > stage1_timing.read_max_seconds) {
        stage1_timing.read_max_seconds = result.max_seconds;
        stage1_timing.slowest_read_input = result.max_input;
        stage1_timing.slowest_read_offset = result.max_offset;
        stage1_timing.slowest_read_bytes = result.max_bytes;
    }
    stage1_timing.usage.samples += result.usage.samples;
    stage1_timing.usage.user_seconds += result.usage.user_seconds;
    stage1_timing.usage.system_seconds += result.usage.system_seconds;
    stage1_timing.usage.minor_faults += result.usage.minor_faults;
    stage1_timing.usage.major_faults += result.usage.major_faults;
    stage1_timing.usage.voluntary_switches += result.usage.voluntary_switches;
    stage1_timing.usage.involuntary_switches += result.usage.involuntary_switches;
    for (i = 0; i < count; ++i) {
        int64_t bytes = (int64_t)inputs[i].length;
        if (bytes == capacity && inputs[i].offset + bytes < aux->input_end[i]) {
            double start = GetTime(), elapsed;
            bytes = get_next_fastq(buffers[i], bytes - (1 << 10), bytes);
            elapsed = GetTime() - start;
            ++stage1_timing.boundary_calls;
            stage1_timing.boundary_seconds += elapsed;
            if (elapsed > stage1_timing.boundary_max_seconds)
                stage1_timing.boundary_max_seconds = elapsed;
        }
        aux->input_position[i] += bytes;
        *sizes[i] = bytes;
    }
}


static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	if (step == 0) {
		double t0 = GetTime();
		double allocation_start;
		ktp_data_t *ret;
		long long block_capacity;

		++stage1_timing.stage_calls;
		allocation_start = GetTime();
		ret = calloc(1, sizeof(*ret));
		record_stage1_allocation(sizeof(*ret),
		                         GetTime() - allocation_start);
		if (ret == NULL)
			err_fatal(__func__, "failed to allocate pipeline batch");

#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
		swbwa_fastq_range_t chunk;
		int next_chunk;

		next_chunk = swbwa_mpi_fastq_scheduler_next(&chunk);
		if (next_chunk < 0)
			err_fatal(__func__, "failed to claim FASTQ chunk: %s", strerror(errno));
		if (next_chunk == 0) {
			free(ret);
			t_step1 += GetTime() - t0;
			return 0;
		}
		aux->input_position[0] = aux->input_position[1] = chunk.start;
		aux->input_end[0] = aux->input_end[1] = chunk.end;
		ret->chunk_id = chunk.chunk_id;
		{
			int64_t read_index;
			int index_scale = aux->is_paired ? 2 : 1;

#if SWBWA_MPI_EXACT_READ_INDEX
			read_index = chunk.first_record;
#else
			/* Fast mode needs only a stable, non-overlapping hash seed. */
			read_index = chunk.start;
#endif
			if (read_index > INT64_MAX / index_scale)
				err_fatal(__func__, "FASTQ read index overflow");
			ret->n_processed = read_index * index_scale;
		}
		block_capacity = chunk.end - chunk.start;
#else
		int64_t remaining_bytes =
			aux->input_end[0] - aux->input_position[0];

		if (remaining_bytes <= 0) {
			free(ret);
			t_step1 += GetTime() - t0;
			return 0;
		}
		block_capacity = aux->fastq_chunk_bytes;
		if (remaining_bytes < block_capacity)
			block_capacity = remaining_bytes;
#endif
		const int is_paired = aux->is_paired;
		char *block_buffer;
		char *block_buffer2;
		long long real_size;
		long long real_size2 = 0;

		if (block_capacity <= 0 ||
		    (uint64_t)block_capacity > SIZE_MAX - 1024)
			err_fatal(__func__, "FASTQ chunk is too large: %lld bytes",
			          block_capacity);
		allocation_start = GetTime();
		block_buffer = malloc((size_t)block_capacity + 1024);
		record_stage1_allocation((uint64_t)block_capacity + 1024,
		                         GetTime() - allocation_start);
		if (is_paired) {
			allocation_start = GetTime();
			block_buffer2 = malloc((size_t)block_capacity + 1024);
			record_stage1_allocation((uint64_t)block_capacity + 1024,
			                         GetTime() - allocation_start);
		} else {
			block_buffer2 = NULL;
		}

		if (block_buffer == NULL || (is_paired && block_buffer2 == NULL))
			err_fatal(__func__, "failed to allocate FASTQ input buffers");
		if (bwa_verbose >= 3)
			fprintf(stderr,
			        "[M::%s] read chunk (%lld bytes/CG x %d CGs,"
			        " %lld bytes)...\n",
			        __func__, (long long)aux->fastq_bytes_per_cg,
			        SWBWA_CG_COUNT, block_capacity);

		ret->chunk_start = aux->input_position[0];
#if !SWBWA_USE_MPI || SWBWA_MPI_INPUT_MODE != SWBWA_MPI_INPUT_DYNAMIC
        ret->chunk_id = ret->chunk_start;
#endif
		if (input_readers) {
			read_fastq_parallel(aux, block_buffer, block_buffer2,
			                    block_capacity, &real_size, &real_size2);
		} else {
			real_size = read_fastq_block(file1_ptr, &aux->input_position[0],
		                             aux->input_end[0], block_buffer,
		                             block_capacity);
			if (is_paired)
				real_size2 = read_fastq_block(
					file2_ptr, &aux->input_position[1], aux->input_end[1],
					block_buffer2, block_capacity);
		}
#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
		if (real_size != block_capacity)
			err_fatal(__func__, "short FASTQ chunk read: %lld != %lld",
			          real_size, block_capacity);
#endif
		if (is_paired) {
			if (real_size != real_size2)
				err_fatal(__func__,
				          "paired FASTQ blocks have different sizes:"
				          " %lld != %lld",
				          real_size, real_size2);
		}
		if (real_size == 0) {
			free(block_buffer);
			free(block_buffer2);
			free(ret);
			t_step1 += GetTime() - t0;
			return 0;
		}
		ret->fastq_buffer[0] = block_buffer;
		ret->fastq_buffer[1] = block_buffer2;
		ret->fastq_size[0] = real_size;
		ret->fastq_size[1] = real_size2;
		ret->chunk_end = aux->input_position[0];
		++stage1_timing.completed_batches;
		t_step1 += GetTime() - t0;
		return ret;
	} else if (step == 1) {
	        double t0 = GetTime();
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
		int64_t n_processed = data->n_processed;
#else
		int64_t n_processed = aux->n_processed;
#endif
		if (opt->flag & MEM_F_SMARTPE) {
			err_fatal(__func__, "smart pairing is not supported by SWBWA");
		} else {
            int **output_lengths = NULL;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
            output_lengths = &data->sam_lengths;
#endif
			mem_process_seqs_merge2(opt, idx->bwt, idx->bns, idx->pac, n_processed, &(data->n_seqs), &(data->seqs), data->fastq_buffer[0], data->fastq_buffer[1], data->fastq_size[0], data->fastq_size[1], aux->pes, output_lengths);
		}
#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
		{
			double stage2_seconds = GetTime() - t0;

			if (aux->is_paired && (data->n_seqs & 1) != 0)
				err_fatal(__func__,
				          "paired FASTQ chunk produced an odd read count");
			swbwa_mpi_fastq_scheduler_record_stage2(
				data->chunk_id,
				data->n_seqs / (aux->is_paired ? 2 : 1),
				stage2_seconds);
			t_step2 += stage2_seconds;
		}
#else
		aux->n_processed += data->n_seqs;
		t_step2 += GetTime() - t0;
#endif
		return data;
	} else if (step == 2) {
        double t0 = GetTime();
		double write_start = GetTime();
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
		for (int i = 0; i < data->n_seqs; ++i) {
			if (data->seqs[i].sam) {
                if (swbwa_output_write_digest(data->seqs[i].sam) != 0)
                    err_fatal(__func__, "failed to write SAM output: %s",
                              strerror(errno));
            }
		}
#else
        if (swbwa_output_write_chunk(data->chunk_id, data->chunk_start,
                                     data->chunk_end, data->seqs, data->sam_lengths,
                                     data->n_seqs, aux->is_paired) != 0)
            err_fatal(__func__, "failed to write output chunk: %s", strerror(errno));
#endif
		t_step3_1 += GetTime() - write_start;
		free(data->sam_lengths);
		free(data->seqs); free(data);
        t_step3 += GetTime() - t0;
		return 0;
	} else if(step == 3) {
        if (swbwa_output_flush() != 0)
            err_fatal(__func__, "failed to flush SAM output: %s", strerror(errno));
    }
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}


static void open_fastq_input(const char *path, FILE **file_ptr, int64_t start)
{
    char resolved_path[PATH_MAX];
    int n = snprintf(resolved_path, sizeof(resolved_path), "%s", path);
    if (n < 0 || (size_t)n >= sizeof(resolved_path))
        err_fatal(__func__, "input path is too long: '%s'", path);
    if (swbwa_mpi_is_root())
        fprintf(stderr, "the input file is %s\n", resolved_path);

    *file_ptr = fopen(resolved_path, "rb");
    if (*file_ptr == NULL)
        err_fatal(__func__, "failed to open input file '%s': %s", resolved_path, strerror(errno));
    if (fseeko(*file_ptr, (off_t)start, SEEK_SET) != 0)
        err_fatal(__func__, "failed to seek input file '%s': %s",
	                  resolved_path, strerror(errno));
}

static int prepare_sequential_fastq_input(const char *read1_path,
                                          const char *read2_path,
                                          ktp_aux_t *aux,
                                          int64_t *file_size)
{
	if (swbwa_fastq_chunk_bytes(
	        read1_path, read2_path, aux->fastq_bytes_per_cg,
	        file_size, &aux->fastq_chunk_bytes) != 0)
		return -1;

	if (swbwa_mpi_is_root() && bwa_verbose >= 3)
		fprintf(stderr,
		        "[CPE input] bytes_per_cg=%lld cg_count=%d"
		        " chunk_bytes=%lld file_bytes=%lld\n",
		        (long long)aux->fastq_bytes_per_cg, SWBWA_CG_COUNT,
		        (long long)aux->fastq_chunk_bytes,
		        (long long)*file_size);
	return 0;
}

static void init_cpe_allocator(void)
{
#if SWBWA_CPE_ALLOC_MODE == SWBWA_CPE_ALLOC_POOL
    static char *pool_buffer;
    swbwa_cpe_pool_params_t params;

    pool_buffer = (char*)malloc(SWBWA_CPE_POOL_TOTAL_BYTES);
    if (pool_buffer == NULL)
        err_fatal(__func__, "failed to allocate CPE malloc buffer: bytes=%lld", (long long)SWBWA_CPE_POOL_TOTAL_BYTES);
    memset(pool_buffer, 0, SWBWA_CPE_POOL_TOTAL_BYTES);
    params.buffer = pool_buffer;
    params.bytes_per_cpe = SWBWA_CPE_POOL_BYTES_PER_CPE;
    swbwa_cpe_run((void*)slave_state_init, &params);
#endif
}

int main_mem(int argc, char *argv[])
{
	mem_opt_t *opt, opt0;
    int i, c, ignore_alt = 0, no_mt_io = 0;
#if SWBWA_OUTPUT_RMA_ONLY
    int saved_report_stderr = -1;
#endif
	int64_t fastq_bytes_per_cg = SWBWA_DEFAULT_FASTQ_BYTES_PER_CG;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	const char *output_path = 0;
	const char *read2_path;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "51qpaMSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:F:z:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'u') opt->flag |= MEM_F_XB;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'o' || c == 'f') output_path = optarg;
        else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
        else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'K') {
			char *end = NULL;
			long long value;

			errno = 0;
			value = strtoll(optarg, &end, 10);
			if (errno == ERANGE || end == optarg || *end != '\0' || value <= 0)
				err_fatal(__func__, "invalid -K value: '%s'", optarg);
			fastq_bytes_per_cg = (int64_t)value;
		}
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'F') bwa_dbg = atoi(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'z') opt->XA_drop_ratio = atof(optarg);
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		} else if (c == 'I') { // specify the insert size distribution
			aux.pes = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3 && swbwa_mpi_is_root())
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}

#if SWBWA_USE_MPI
	if (!swbwa_mpi_is_root() && bwa_verbose > 1 && bwa_verbose < 4)
		bwa_verbose = 1;
#endif

	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		if (!swbwa_mpi_is_root()) {
			free(opt);
			return 1;
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: SWBWA mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -o FILE       sam file to output results to [stdout]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "       -5            for split alignment, take the alignment with the smallest query (not genomic) coordinate as primary\n");
		fprintf(stderr, "       -q            don't modify mapQ of supplementary alignments\n");
		fprintf(stderr,
		        "       -K INT        target raw FASTQ bytes per CG and input file"
		        " [%lld]\n",
		        (long long)SWBWA_DEFAULT_FASTQ_BYTES_PER_CG);
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >%.2f%% of the max score, output all in XA [%d,%d]\n", 
				opt->XA_drop_ratio * 100.0,
				opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "                     A second value may be given for alternate sequences.\n");
		fprintf(stderr, "       -z FLOAT      The fraction of the max score to use with -h [%f].\n", opt->XA_drop_ratio);
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "       -u            output XB instead of XA; XB is XA with the alignment score and mapping quality added.\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	aux.fastq_bytes_per_cg = fastq_bytes_per_cg;
	read2_path = (optind + 2 < argc && !(opt->flag & MEM_F_PE))
		? argv[optind + 2] : NULL;

	aux.input_position[0] = aux.input_position[1] = 0;
	aux.input_end[0] = aux.input_end[1] = INT64_MAX;
#if SWBWA_USE_MPI
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
	{
		swbwa_fastq_range_t assigned_range;
		int64_t chunk_count;

		if (swbwa_mpi_fastq_scheduler_open(
		        argv[optind + 1], read2_path, aux.fastq_bytes_per_cg,
		        bwa_verbose >= 4,
		        &assigned_range, &chunk_count) != 0)
			return 1;
		aux.input_position[0] = aux.input_position[1] = assigned_range.start;
		aux.input_end[0] = aux.input_end[1] = assigned_range.end;
		if (swbwa_mpi_is_root())
			fprintf(stderr,
			        "[MPI input] mode=dynamic chunks=%lld bytes_per_cg=%lld"
			        " cg_count=%d chunk_bytes=%lld micro_chunk_bytes=%lld"
			        " fine_chunk_bytes=%lld tail_percent=%d"
			        " fine_tail_waves=%d file_bytes=%lld"
			        " scheduler=%s"
#if SWBWA_MPI_EXACT_READ_INDEX
			        " read_index=exact\n",
#else
			        " read_index=byte_offset\n",
#endif
			        (long long)chunk_count,
			        (long long)aux.fastq_bytes_per_cg, SWBWA_CG_COUNT,
			        (long long)swbwa_mpi_fastq_scheduler_chunk_bytes(),
			        (long long)swbwa_mpi_fastq_scheduler_micro_chunk_bytes(),
			        (long long)swbwa_mpi_fastq_scheduler_fine_chunk_bytes(),
			        swbwa_mpi_fastq_scheduler_tail_percent(),
			        swbwa_mpi_fastq_scheduler_fine_tail_waves(),
			        (long long)assigned_range.file_size,
			        swbwa_mpi_fastq_scheduler_ticket_mode());
	}
#else
	{
		int64_t file_size;
		swbwa_fastq_range_t range;

		if (prepare_sequential_fastq_input(
		        argv[optind + 1], read2_path, &aux, &file_size) != 0)
			return 1;
		if (swbwa_mpi_fastq_range(argv[optind + 1], read2_path, &range) != 0)
			return 1;
		aux.input_position[0] = aux.input_position[1] = range.start;
		aux.input_end[0] = aux.input_end[1] = range.end;
#if SWBWA_MPI_EXACT_READ_INDEX
		aux.n_processed =
			range.first_record * (read2_path != NULL ? 2 : 1);
		fprintf(stderr,
		        "[MPI rank %06d/%06d] input range: records=[%lld, %lld)"
		        " bytes=[%lld, %lld)\n",
		        swbwa_mpi_rank(), swbwa_mpi_size(),
		        (long long)range.first_record,
		        (long long)(range.first_record + range.record_count),
		        (long long)range.start, (long long)range.end);
#else
		fprintf(stderr,
		        "[MPI rank %06d/%06d] input range: bytes=[%lld, %lld)\n",
		        swbwa_mpi_rank(), swbwa_mpi_size(),
		        (long long)range.start, (long long)range.end);
#endif
	}
#endif
#else
	{
		int64_t file_size;

		if (prepare_sequential_fastq_input(
		        argv[optind + 1], read2_path, &aux, &file_size) != 0)
			return 1;
		aux.input_end[0] = aux.input_end[1] = file_size;
	}
#endif

	swbwa_runtime_init();

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the SWBWA index from shared memory\n", __func__);
	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;

	open_fastq_input(argv[optind + 1], &file1_ptr,
	                 aux.input_position[0]);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		} else {
			open_fastq_input(argv[optind + 2], &file2_ptr,
			                 aux.input_position[1]);
			aux.is_paired = 1;
			opt->flag |= MEM_F_PE;
		}
	}
#if SWBWA_USE_MPI
	if (output_path == NULL || strcmp(output_path, "-") == 0)
		err_fatal(__func__, "MPI output requires an explicit -o FILE path");
#endif
	if (swbwa_output_open(output_path, bwa_verbose >= 4) != 0)
		err_fatal(__func__, "failed to open SAM output: %s", strerror(errno));
	if (output_path != NULL || SWBWA_USE_MPI) {
#if SWBWA_USE_MPI
		fprintf(stderr, "[MPI rank %06d/%06d] output: %s\n",
		        swbwa_mpi_rank(), swbwa_mpi_size(), swbwa_output_name());
#else
		fprintf(stderr, "the output file is %s\n", swbwa_output_name());
#endif
	}

	{
		const char *value = getenv("SWBWA_INPUT_READERS");
		input_readers = swbwa_input_reader_count(value);
		if (input_readers < 0)
			err_fatal(__func__, "SWBWA_INPUT_READERS must be 0..%d (0=fread)",
			          SWBWA_HOST_MPE_THREADS);
	}
	if (swbwa_host_workers_init() != 0)
        err_fatal(__func__, "failed to initialize MPE helpers: %s", strerror(errno));
	init_cpe_allocator();
	swbwa_cpe_profile_init();


#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
	if (swbwa_mpi_barrier() != 0)
		err_fatal(__func__, "failed to synchronize dynamic MPI workers");
#endif

	swbwa_cpe_progress_debug_enable(bwa_verbose >= 4);
#if SWBWA_USE_MPI && \
    (SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC || \
     SWBWA_OUTPUT_SINGLE_FILE)
	if (swbwa_mpi_progress_thread_start(bwa_verbose >= 4) != 0)
		err_fatal(__func__, "failed to start MPI progress thread: %s",
		          strerror(errno));
#endif
	double t0 = GetTime();
	if (no_mt_io) kt_pipeline_single(1, process, &aux, 3);
	else kt_pipeline_queue(3, process, &aux, 3);
    /* Stage 3 has flushed; exclude rank-ordered reports and MPI teardown. */
    t_tot += GetTime() - t0;
    swbwa_host_workers_destroy();
#if SWBWA_OUTPUT_RMA_ONLY
    {
        const char *prefix = getenv("SWBWA_RANK_REPORT_PREFIX");
        if (prefix != NULL && *prefix != '\0') {
            char path[PATH_MAX];
            int report_fd;
            int length = snprintf(path, sizeof(path), "%s.rank%06d.log",
                                  prefix, swbwa_mpi_rank());
            if (length < 0 || (size_t)length >= sizeof(path))
                err_fatal(__func__, "rank report path is too long");
            err_fflush(stderr);
            /* Keep the launcher's stderr pipe alive until normal shutdown. */
            saved_report_stderr = dup(STDERR_FILENO);
            if (saved_report_stderr < 0)
                err_fatal(__func__, "failed to save stderr: %s", strerror(errno));
            report_fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0666);
            if (report_fd < 0)
                err_fatal(__func__, "failed to open rank report: %s", strerror(errno));
            if (dup2(report_fd, STDERR_FILENO) < 0)
                err_fatal(__func__, "failed to redirect rank report: %s", strerror(errno));
            close(report_fd);
        }
    }
#endif
#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
	{
		int64_t chunks, records, bytes;

		swbwa_mpi_fastq_scheduler_stats(&chunks, &records, &bytes);
		fprintf(stderr,
		        "[MPI rank %06d/%06d] input work: chunks=%lld records=%lld"
		        " bytes=%lld\n",
		        swbwa_mpi_rank(), swbwa_mpi_size(),
		        (long long)chunks, (long long)records,
		        (long long)bytes);
	}
	swbwa_mpi_fastq_scheduler_close();
#endif
	if (swbwa_output_close() != 0)
		err_fatal(__func__, "failed to close SAM output: %s", strerror(errno));
#if SWBWA_USE_MPI && \
    (SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC || \
     SWBWA_OUTPUT_SINGLE_FILE)
	if (swbwa_mpi_progress_thread_stop() != 0)
		err_fatal(__func__, "MPI progress thread failed: %s", strerror(errno));
#endif

	swbwa_mpi_progress_thread_report();
	swbwa_cpe_progress_debug_report();
	swbwa_cpe_profile_report(stderr);
    print_timing_report();
    if (bwa_verbose >= 4 && !no_mt_io) kt_pipeline_wait_report();
#if SWBWA_OUTPUT_RMA_ONLY
    if (saved_report_stderr >= 0) {
        err_fflush(stderr);
        if (dup2(saved_report_stderr, STDERR_FILENO) < 0)
            err_fatal(__func__, "failed to restore stderr: %s", strerror(errno));
        close(saved_report_stderr);
    }
#endif

	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
    if (file1_ptr) fclose(file1_ptr);
	if (aux.is_paired) {
        if (file2_ptr) fclose(file2_ptr);
	}
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX;
	uint64_t max_intv = 0;
	kseq_t *seq;
	bwtint_t k;
	int fd;
	void *ko;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:pi:I:L:")) >= 0) {
		switch (c) {
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'i': min_intv = atoi(optarg); break;
			case 'I': max_intv = atol(optarg); break;
			case 'L': max_len  = atoi(optarg); break;
		    default: return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   SWBWA fastmap [options] <idxbase> <in.fq>\n\n");
		fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
		fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
		fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
		fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
		fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
		fprintf(stderr, "\n");
		return 1;
	}

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	seq = kseq_init(fd);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	itr = smem_itr_init(idx->bwt);
	smem_config(itr, min_intv, max_len, max_intv);
	while (kseq_read(seq) >= 0) {
		err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			err_putchar('\t');
			err_puts(seq->seq.s);
		} else err_putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				err_printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_puts("\t*");
				err_putchar('\n');
			}
		}
		err_puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	kclose(ko);
	return 0;
}
