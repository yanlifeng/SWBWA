#include <slave.h>
#include <crts.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "bwt.h"
#include "bwamem.h"
#include "../swbwa_config.h"
#include "../swbwa_cpe.h"

#include "lwpf3_my_cpe.h"
#include "malloc_wrap.h"



/***********************************************************/
extern unsigned long swbwa_cpe_text_start;
extern unsigned long swbwa_cpe_text_size;
extern unsigned long swbwa_cpe_data_start;
extern unsigned long swbwa_cpe_data_size;
/***********************************************************/


__thread swbwa_cpe_task_t *swbwa_task;

static inline void swbwa_enter_cross_runtime(void)
{
    unsigned long csr_value;

    asm volatile("mov %0, $29\n\t"::"r"(swbwa_task->relocated_gp):);
    asm volatile("rcsr %0, 0xc4" : "=r"(csr_value));
    if (csr_value < 0x500000000000) {
        unsigned long relocated_csr =
            (unsigned long)swbwa_task->private_segment_copies +
            SWBWA_CPE_CSR_COPY_BYTES * _MYID + csr_value -
            SWBWA_CPE_PRIVATE_BASE;
        asm volatile("wcsr %0, 0xc4\n\t"::"r"(relocated_csr):);
    }
}

static inline void swbwa_finish_cross_task(void)
{
    swbwa_task->completion_flags[_MYID] = 1;
    flush_slave_cache();
    while (1) { }
}

void state_init(swbwa_cpe_pool_params_t *para) {
    set_big_buffer(para->buffer, para->bytes_per_cpe);
}

void pass_para(swbwa_cpe_task_t *ptr)
{
    swbwa_task = ptr;
}

void copy_priv_segment(swbwa_cpe_task_t *ptr) {
    void* my_priv_addr = ptr->private_segment_copies + SWBWA_CPE_CSR_COPY_BYTES * _MYID;
    memcpy(my_priv_addr, (void*)SWBWA_CPE_PRIVATE_BASE, SWBWA_CPE_CSR_COPY_BYTES);
}

void change_priv_segment(swbwa_cpe_task_t *ptr) {
    int tls_size = ptr->tls_relocation_count;
    unsigned long *tls_content = ptr->tls_relocations;
    void* my_priv_addr = ptr->private_segment_copies + SWBWA_CPE_CSR_COPY_BYTES * _MYID;
    int pre_instr_size = 8;
    for(unsigned long i = 0; i < (SWBWA_CPE_CSR_COPY_BYTES) / pre_instr_size; ++i) {
        unsigned long disp = i * pre_instr_size;
        unsigned long *addr = (unsigned long*)((unsigned long)my_priv_addr + disp);
        unsigned long content = *addr;
        if(content >= swbwa_cpe_text_start + swbwa_cpe_text_size) continue;
        if(content < swbwa_cpe_text_start) continue;

        /* Relocate CPE TLS references that point into the copied text. */
        for(int j = 0; j < tls_size; j++) {
            unsigned long content2 = *((unsigned long *)tls_content[j]);
            if(content == content2) {
                *addr = content + ptr->segment_offset;
                break;
            }
        }
    }
}

#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
__uncached long work_counter;
int task_list[SWBWA_CPE_COUNT][SWBWA_MAX_TASKS_PER_CPE];
int task_num[SWBWA_CPE_COUNT];
long cur_id[SWBWA_CPE_COUNT];

int acquire_task(int block_num) {
    asm volatile("faal %0, 0(%1)\n\t"
                 : "=r"(cur_id[_MYID])
                 : "r"(&work_counter)
                 : "memory");
    if (cur_id[_MYID] < block_num) {
        assert(task_num[_MYID] < SWBWA_MAX_TASKS_PER_CPE);
        task_list[_MYID][task_num[_MYID]++] = cur_id[_MYID];
    }
    return (int)cur_id[_MYID];
}
#endif

void worker1_s_pre_fast_cross(void) {
    swbwa_enter_cross_runtime();
    asm volatile("memb\n\t":::);
    swbwa_cpe_task_t *para = swbwa_task;

    lwpf_enter(TEST);
    lwpf_start(l_worker1_1);

#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    if(_MYID == 0) work_counter = 0;
# if SWBWA_USE_CGS
    athread_ssync_node();
# else
    athread_ssync_array();
# endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int range_l = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker1_pre_fast(para->worker_data, j, _MYID, para->alignment_regions);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker1_pre_fast(para->worker_data, i, _MYID, para->alignment_regions);
    }
#endif

    lwpf_stop(l_worker1_1);
    lwpf_exit(TEST);

    swbwa_finish_cross_task();
}

void worker1_s_fast_cross(void) {
    swbwa_enter_cross_runtime();
    swbwa_cpe_task_t *para = swbwa_task;

    lwpf_enter(TEST);
    lwpf_start(l_worker1_2);

#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int range_l = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker1_fast(para->worker_data, j, _MYID, para->alignment_regions);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker1_fast(para->worker_data, i, _MYID, para->alignment_regions);
    }
#endif
    lwpf_stop(l_worker1_2);
    lwpf_exit(TEST);

    swbwa_finish_cross_task();
}


void worker2_s_pre_fast_cross(void) {
    swbwa_enter_cross_runtime();
    swbwa_cpe_task_t *para = swbwa_task;

    lwpf_enter(TEST);
    lwpf_start(l_worker2_1);

#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    if(_MYID == 0) work_counter = 0;
#if SWBWA_USE_CGS
    athread_ssync_node();
#else
    athread_ssync_array();
#endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int range_l = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker2_pre_fast(para->worker_data, j, _MYID, para->sam_lengths, para->sam_records);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker2_pre_fast(para->worker_data, i, _MYID, para->sam_lengths, para->sam_records);
    }
#endif
    lwpf_stop(l_worker2_1);
    lwpf_exit(TEST);


    swbwa_finish_cross_task();
}


void worker2_s_fast_cross(void) {
    swbwa_enter_cross_runtime();
    swbwa_cpe_task_t *para = swbwa_task;
 
    lwpf_enter(TEST);
    lwpf_start(l_worker2_2);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int range_l = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker2_fast(para->worker_data, j, _MYID, para->sam_lengths, para->sam_records);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker2_fast(para->worker_data, i, _MYID, para->sam_lengths, para->sam_records);
    }
#endif
    lwpf_stop(l_worker2_2);
    lwpf_exit(TEST);

 
    swbwa_finish_cross_task();
}

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
    assert(pos_ < size_ && data_[pos_] == '+');
    return pos0;
}

static int format_fastq_partition(swbwa_cpe_task_t *params)
{
    long long input_chunk_size;
    long long input_begin;
    long long input_end;
    long long output_chunk_size;
    long long output_begin;
    long long output_end;
    int is_paired = params->fastq_size[1] > 0;

    if (params->fastq_size[0] < SWBWA_CPE_COUNT * (1 << 10)) {
        if (_MYID == 0) {
            return format_seqs(params->fastq_buffer[0], params->fastq_size[0],
                               params->fastq_buffer[1], params->fastq_size[1],
                               params->formatted_buffer[0], params->formatted_buffer[1],
                               params->formatted_buffer_size, params->worker_data);
        }
        return format_seqs(NULL, 0, NULL, 0, NULL, NULL, 0,
                           params->worker_data);
    }

    input_chunk_size = (params->fastq_size[0] + SWBWA_CPE_COUNT - 1) /
                       SWBWA_CPE_COUNT;
    input_begin = _MYID * input_chunk_size;
    input_end = input_begin + input_chunk_size;
    if (input_end > params->fastq_size[0]) input_end = params->fastq_size[0];

    if (_MYID > 0)
        input_begin = get_next_fastq(params->fastq_buffer[0],
                                     input_begin - (1 << 10),
                                     params->fastq_size[0]);
    if (_MYID < SWBWA_CPE_COUNT - 1)
        input_end = get_next_fastq(params->fastq_buffer[0],
                                   input_end - (1 << 10),
                                   params->fastq_size[0]);

    output_chunk_size = (params->formatted_buffer_size + SWBWA_CPE_COUNT - 1) /
                        SWBWA_CPE_COUNT;
    output_begin = _MYID * output_chunk_size;
    output_end = output_begin + output_chunk_size;
    if (output_end > params->formatted_buffer_size)
        output_end = params->formatted_buffer_size;

    return format_seqs(params->fastq_buffer[0] + input_begin,
                       input_end - input_begin,
                       is_paired ? params->fastq_buffer[1] + input_begin : NULL,
                       is_paired ? input_end - input_begin : 0,
                       params->formatted_buffer[0] + output_begin,
                       is_paired ? params->formatted_buffer[1] + output_begin : NULL,
                       output_end - output_begin, params->worker_data);
}

void cpe_format_pre_cross(void)
{
    swbwa_enter_cross_runtime();
    swbwa_task->formatted_read_counts[_MYID] =
        format_fastq_partition(swbwa_task);
    swbwa_finish_cross_task();
}

void worker12_s_pre_fast_cross(void) {
    swbwa_enter_cross_runtime();
    swbwa_cpe_task_t *para = swbwa_task;
    lwpf_enter(TEST);
    lwpf_start(l_worker12_1);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    if(_MYID == 0) work_counter = 0;
#if SWBWA_USE_CGS
    athread_ssync_node();
#else
    athread_ssync_array();
#endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int l_pos = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int r_pos = l_pos + SWBWA_READS_PER_DYNAMIC_TASK;
        if(r_pos > para->work_item_count) r_pos = para->work_item_count;
        lwpf_start(l_worker12i_1);
        worker12_pre_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->pes, para->sequence_ids);
        lwpf_stop(l_worker12i_1);
    }
#else
    int pre_n = ceil(1.0 * para->work_item_count / SWBWA_CPE_COUNT);
    int l_pos = _MYID * pre_n;
    int r_pos = l_pos + pre_n;
    if(r_pos > para->work_item_count) r_pos = para->work_item_count;
    worker12_pre_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->pes, para->sequence_ids);
#endif

    lwpf_stop(l_worker12_1);
    lwpf_exit(TEST);

    swbwa_finish_cross_task();
}


void worker12_s_fast_cross(void) {
    swbwa_enter_cross_runtime();
    swbwa_cpe_task_t *para = swbwa_task;
    lwpf_enter(TEST);
    lwpf_start(l_worker12_2);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int l_pos = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int r_pos = l_pos + SWBWA_READS_PER_DYNAMIC_TASK;
        if(r_pos > para->work_item_count) r_pos = para->work_item_count;
        worker12_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->sequence_ids);
    }
#else
    int pre_n = ceil(1.0 * para->work_item_count / SWBWA_CPE_COUNT);
    int l_pos = _MYID * pre_n;
    int r_pos = l_pos + pre_n;
    if(r_pos > para->work_item_count) r_pos = para->work_item_count;
    worker12_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->sequence_ids);
#endif

    lwpf_stop(l_worker12_2);
    lwpf_exit(TEST);

    swbwa_finish_cross_task();
}


void worker1_s_pre_fast(swbwa_cpe_task_t *para) {
 
    lwpf_enter(TEST);
    lwpf_start(l_worker1_1);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
# if SWBWA_USE_CGS
    if(_MYID == 0) work_counter = 0;
    athread_ssync_node();
# else
    if(_PEN == 0) work_counter = 0;
    athread_ssync_array();
# endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int range_l = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker1_pre_fast(para->worker_data, j, _MYID, para->alignment_regions);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker1_pre_fast(para->worker_data, i, _MYID, para->alignment_regions);
    }
#endif
    lwpf_stop(l_worker1_1);
    lwpf_exit(TEST);
}

void worker1_s_fast(swbwa_cpe_task_t *para) {
    lwpf_enter(TEST);
    lwpf_start(l_worker1_2);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int range_l = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker1_fast(para->worker_data, j, _MYID, para->alignment_regions);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker1_fast(para->worker_data, i, _MYID, para->alignment_regions);
    }
#endif
    lwpf_stop(l_worker1_2);
    lwpf_exit(TEST);

}



void worker2_s_pre_fast(swbwa_cpe_task_t *para) {
    lwpf_enter(TEST);
    lwpf_start(l_worker2_1);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    if(_MYID == 0) work_counter = 0;
#if SWBWA_USE_CGS
    athread_ssync_node();
#else
    athread_ssync_array();
#endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int range_l = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker2_pre_fast(para->worker_data, j, _MYID, para->sam_lengths, para->sam_records);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker2_pre_fast(para->worker_data, i, _MYID, para->sam_lengths, para->sam_records);
    }
#endif
    lwpf_stop(l_worker2_1);
    lwpf_exit(TEST);

}


void worker2_s_fast(swbwa_cpe_task_t *para) {
    lwpf_enter(TEST);
    lwpf_start(l_worker2_2);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int range_l = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int range_r = range_l + SWBWA_READS_PER_DYNAMIC_TASK;
        if(range_r > para->work_item_count) range_r = para->work_item_count;
        for(int j = range_l; j < range_r; j++) {
            worker2_fast(para->worker_data, j, _MYID, para->sam_lengths, para->sam_records);
        }
    }
#else
    for(long i = _MYID; i < para->work_item_count; i += SWBWA_CPE_COUNT) {
        worker2_fast(para->worker_data, i, _MYID, para->sam_lengths, para->sam_records);
    }
#endif
    lwpf_stop(l_worker2_2);
    lwpf_exit(TEST);

}

void cpe_format_pre(swbwa_cpe_task_t *para)
{
    para->formatted_read_counts[_MYID] =
        format_fastq_partition(para);
}


void worker12_s_pre_fast(swbwa_cpe_task_t *para) {
    lwpf_enter(TEST);
    lwpf_start(l_worker12_1);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    if(_MYID == 0) work_counter = 0;
#if SWBWA_USE_CGS
    athread_ssync_node();
#else
    athread_ssync_array();
#endif
    task_num[_MYID] = 0;
    int block_num = ceil(1.0 * para->work_item_count / SWBWA_READS_PER_DYNAMIC_TASK);
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        int l_pos = i * SWBWA_READS_PER_DYNAMIC_TASK;
        int r_pos = l_pos + SWBWA_READS_PER_DYNAMIC_TASK;
        if(r_pos > para->work_item_count) r_pos = para->work_item_count;
        lwpf_start(l_worker12i_1);
        worker12_pre_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->pes, para->sequence_ids);
        lwpf_stop(l_worker12i_1);
    }
#else
    int pre_n = ceil(1.0 * para->work_item_count / SWBWA_CPE_COUNT);
    int l_pos = _MYID * pre_n;
    int r_pos = l_pos + pre_n;
    if(r_pos > para->work_item_count) r_pos = para->work_item_count;
    worker12_pre_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->pes, para->sequence_ids);
#endif
    lwpf_stop(l_worker12_1);
    lwpf_exit(TEST);
}


void worker12_s_fast(swbwa_cpe_task_t *para) {
    lwpf_enter(TEST);
    lwpf_start(l_worker12_2);
#if SWBWA_ENABLE_DYNAMIC_SCHEDULING
    for(long i = 0; i < task_num[_MYID]; i++) {
        int l_pos = task_list[_MYID][i] * SWBWA_READS_PER_DYNAMIC_TASK;
        int r_pos = l_pos + SWBWA_READS_PER_DYNAMIC_TASK;
        if(r_pos > para->work_item_count) r_pos = para->work_item_count;
        worker12_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->sequence_ids);
    }
#else
    int pre_n = ceil(1.0 * para->work_item_count / SWBWA_CPE_COUNT);
    int l_pos = _MYID * pre_n;
    int r_pos = l_pos + pre_n;
    if(r_pos > para->work_item_count) r_pos = para->work_item_count;
    worker12_fast(para->worker_data, l_pos, r_pos, _MYID, para->sam_lengths, para->sam_records, para->sequence_ids);
#endif
    lwpf_stop(l_worker12_2);
    lwpf_exit(TEST);
}
