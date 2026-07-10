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
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include "swbwa_config.h"
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "swbwa_cpe.h"
#include "swbwa_runtime.h"

#include <athread.h>
//KSEQ_DECLARE(gzFile)
KSEQ_DECLARE(int)

extern unsigned char nst_nt4_table[256];

#ifndef O_DIRECT
#define O_DIRECT 0
#endif

#if SWBWA_ENABLE_LWPF
#define LWPF_UNITS U(TEST)
#include "lwpf.h"
#endif

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

long long s_reg_sum = 0;
long long c_px2 = 0;
long long s_px2 = 0;

static void print_timing_line(const char *label, double seconds, int indent)
{
    fprintf(stderr, "%*s%-52s %10.3f s\n", indent, "", label, seconds);
}

static void print_timing_report(void)
{
    double stage3_cleanup = t_step3 - t_step3_1;

    if (stage3_cleanup < 0.0) stage3_cleanup = 0.0;

    fprintf(stderr,
            "\n"
            "======================= SWBWA Timing Report =======================\n"
            "  Accumulated wall-clock time; indented rows are included in\n"
            "  their parent row and must not be added to it.\n"
            "\n"
            "  Pipeline\n");
    print_timing_line("total - complete three-stage pipeline", t_tot, 4);
#if SWBWA_ENABLE_CPE_FORMAT
    print_timing_line("stage 1 - allocate and read raw FASTQ blocks", t_step1, 4);
#else
    print_timing_line("stage 1 - read and format FASTQ records", t_step1, 4);
    print_timing_line("bseq_read - parse and materialize records", t_step1_1, 6);
    print_timing_line("file read system calls", t_step1_read, 8);
    print_timing_line("allocate/copy bseq record fields", t_step1_1_1, 8);
#endif
    print_timing_line("stage 2 - align reads and generate SAM records", t_step2, 4);
    print_timing_line("stage 3 - write SAM and release batch data", t_step3, 4);
    print_timing_line("SAM output writes", t_step3_1, 6);
    print_timing_line("batch cleanup and loop overhead (derived)", stage3_cleanup, 6);

    fprintf(stderr, "\n  Stage 2 worker details\n");
    print_timing_line("total - active merge worker", t_work1, 4);
#if SWBWA_ENABLE_CPE_FORMAT
    print_timing_line("part 1 - prepare CPE task and reusable buffers", t_work1_1, 6);
    print_timing_line("part 2 - CPE FASTQ formatting and input release", t_work1_2, 6);
    print_timing_line("part 3 - CPE alignment and SAM length pass", t_work1_3, 6);
    print_timing_line("part 4 - assign slices in the shared SAM buffer", t_work1_4, 6);
    print_timing_line("part 5 - CPE SAM record generation", t_work1_5, 6);
    print_timing_line("part 6 - release temporary worker data", t_work1_6, 6);
#else
    print_timing_line("part 1 - prepare CPE task and runtime", t_work1_1, 6);
    print_timing_line("part 2 - CPE alignment and SAM length pass", t_work1_2, 6);
    print_timing_line("part 3 - allocate and assign packed SAM blocks", t_work1_3, 6);
    print_timing_line("part 4 - CPE SAM record generation", t_work1_4, 6);
    print_timing_line("part 5 - release temporary worker data", t_work1_5, 6);
#endif

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

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_single(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_thread(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void kt_pipeline_queue(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);


static int file_out_fd = -1;
static char sam_out_file_name[256];

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
    char *fastq_buffer[2];
    long long fastq_size[2];
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

static long long read_fastq_block(FILE *file, long long *file_offset,
                                  char *buffer, long long capacity)
{
    long long bytes_read;
    long long record_bytes;

    if (fseek(file, *file_offset, SEEK_SET) != 0)
        err_fatal(__func__, "failed to seek FASTQ input: %s", strerror(errno));

    bytes_read = fread(buffer, 1, capacity, file);
    if (ferror(file))
        err_fatal(__func__, "failed to read FASTQ input: %s", strerror(errno));

    record_bytes = bytes_read;
    if (bytes_read == capacity)
        record_bytes = get_next_fastq(buffer, bytes_read - (1 << 10), bytes_read);
    *file_offset += record_bytes;
    return record_bytes;
}


static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
        double t0 = GetTime();
#if SWBWA_ENABLE_CPE_FORMAT
        static long long now_file1_pos = 0;
        static long long now_file2_pos = 0;
        const long long total_block_size =
            SWBWA_CPE_COUNT * SWBWA_READS_PER_CPE_BLOCK *
            SWBWA_ESTIMATED_FASTQ_RECORD_BYTES;
        const int is_paired = aux->ks2 != NULL;
        char *block_buffer = malloc(total_block_size + 1024);
        char *block_buffer2 = is_paired ? malloc(total_block_size + 1024) : NULL;
        long long real_size;
        long long real_size2 = 0;
        ktp_data_t *ret = calloc(1, sizeof(*ret));

        if (ret == NULL || block_buffer == NULL || (is_paired && block_buffer2 == NULL))
            err_fatal(__func__, "failed to allocate FASTQ input buffers");
        if (bwa_verbose >= 3)
            fprintf(stderr, "[M::%s] read block (%lld)...\n", __func__, total_block_size);

        real_size = read_fastq_block(file1_ptr, &now_file1_pos,
                                     block_buffer, total_block_size);
        if (is_paired) {
            real_size2 = read_fastq_block(file2_ptr, &now_file2_pos,
                                          block_buffer2, total_block_size);
            if (real_size != real_size2)
                err_fatal(__func__, "paired FASTQ blocks have different sizes: %lld != %lld",
                          real_size, real_size2);
        }
        if (real_size == 0) {
            free(block_buffer);
            free(block_buffer2);
            free(ret);
            return 0;
        }
        ret->fastq_buffer[0] = block_buffer;
        ret->fastq_buffer[1] = block_buffer2;
        ret->fastq_size[0] = real_size;
        ret->fastq_size[1] = real_size2;
#else
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
        double t1 = GetTime();
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
        t_step1_1 += GetTime() - t1;
		if (ret->seqs == 0) {
			free(ret);
            t_step1 += GetTime() - t0;
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
#endif
        t_step1 += GetTime() - t0;
		return ret;
	} else if (step == 1) {
        double t0 = GetTime();
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		if (opt->flag & MEM_F_SMARTPE) {
			err_fatal(__func__, "smart pairing is not supported by SWBWA");
		} else {

#if SWBWA_ENABLE_CPE_FORMAT
            mem_process_seqs_merge2(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, &(data->n_seqs), &(data->seqs), data->fastq_buffer[0], data->fastq_buffer[1], data->fastq_size[0], data->fastq_size[1], aux->pes);
#else
            mem_process_seqs_merge(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes);
#endif
        }
		aux->n_processed += data->n_seqs;
        t_step2 += GetTime() - t0;
		return data;
	} else if (step == 2) {
        double t0 = GetTime();
		for (i = 0; i < data->n_seqs; ++i) {
            double t1 = GetTime();
            if (data->seqs[i].sam) {
                if (file_out_fd >= 0) my_align_write(data->seqs[i].sam, file_out_fd, 0, NULL);
                else err_fputs(data->seqs[i].sam, stdout);
            }
            t_step3_1 += GetTime() - t1;
#if !SWBWA_ENABLE_CPE_FORMAT
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual);
            if (data->seqs[i].is_new_addr == 1) {
				free(data->seqs[i].sam);
			}
#endif
		}
		free(data->seqs); free(data);
        t_step3 += GetTime() - t0;
		return 0;
	} else if(step == 3) {
        if (file_out_fd >= 0) my_align_write(NULL, file_out_fd, 1, sam_out_file_name);
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


static void open_sam_output(const char *path)
{
    int n = snprintf(sam_out_file_name, sizeof(sam_out_file_name), "%s", path);
    if (n < 0 || (size_t)n >= sizeof(sam_out_file_name))
        err_fatal(__func__, "output path is too long: '%s'", path);
    file_out_fd = open(sam_out_file_name, O_CREAT | O_RDWR | O_DIRECT, 0666);
    if (file_out_fd == -1)
        err_fatal(__func__, "failed to open output file '%s': %s", sam_out_file_name, strerror(errno));
    fprintf(stderr, "the output file is %s\n", sam_out_file_name);
}

static void open_fastq_input(const char *path, FILE **file_ptr, void **ko, int *fd)
{
    char resolved_path[PATH_MAX];
    int n = snprintf(resolved_path, sizeof(resolved_path), "%s", path);
    if (n < 0 || (size_t)n >= sizeof(resolved_path))
        err_fatal(__func__, "input path is too long: '%s'", path);
    fprintf(stderr, "the input file is %s\n", resolved_path);

    *file_ptr = fopen(resolved_path, "r");
    if (*file_ptr == NULL)
        err_fatal(__func__, "failed to open input file '%s': %s", resolved_path, strerror(errno));

    *ko = kopen(resolved_path, fd);
    if (*ko == NULL)
        err_fatal(__func__, "failed to open input stream '%s'", resolved_path);
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
	int fd, fd2;
    int i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	//gzFile fp, fp2 = 0;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:F:z:")) >= 0) {
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
		else if (c == 'o' || c == 'f') {
			open_sam_output(optarg);
        }
        //xreopen(optarg, "wb", stdout);
        else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
        else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
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
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}

	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: SWBWA mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
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
		fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
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
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
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

	swbwa_runtime_init();

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the SWBWA index from shared memory\n", __func__);
	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;

	open_fastq_input(argv[optind + 1], &file1_ptr, &ko, &fd);
	//fp = gzdopen(fd, "r");
	//aux.ks = kseq_init(fp);
	aux.ks = kseq_init(fd);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		} else {
			open_fastq_input(argv[optind + 2], &file2_ptr, &ko2, &fd2);
			//fp2 = gzdopen(fd2, "r");
			//aux.ks2 = kseq_init(fp2);
			aux.ks2 = kseq_init(fd2);
			opt->flag |= MEM_F_PE;
		}
	}
	//bwa_print_sam_hdr(aux.idx->bns, hdr_line);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;

	init_cpe_allocator();

#if SWBWA_ENABLE_LWPF
    lwpf_init(NULL);
#endif

    double t0 = GetTime();
    if(no_mt_io) kt_pipeline_single(1, process, &aux, 3);
    else kt_pipeline_queue(3, process, &aux, 3);

#if SWBWA_ENABLE_LWPF
    char filename[50];
    snprintf(filename, sizeof(filename), "lwpf_%dCG.log", SWBWA_CG_COUNT);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open lwpf.log");
        return 1;
    }
    lwpf_report_summary(file);
    lwpf_report_detail(file);
    fclose(file);
#endif

    t_tot += GetTime() - t0;


    print_timing_report();


	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	//err_gzclose(fp);
    kclose(ko);
    if (file1_ptr) fclose(file1_ptr);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		//err_gzclose(fp2);
        kclose(ko2);
        if (file2_ptr) fclose(file2_ptr);
	}
    //if(file_out_fd != -1) close(file_out_fd);
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
