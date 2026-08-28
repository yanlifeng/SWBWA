#include <inttypes.h>
#include <string.h>

#include "swbwa_config.h"
#include "swbwa_cpe_profile.h"

#if SWBWA_ENABLE_CPE_PROFILE

#define LWPF_NOCOLOR
#define LWPF_UNITS U(SWBWA)
#include "lwpf.h"

enum { SWBWA_LWPF_CPE_COUNT = 64 };

static swbwa_matesw_profile_t matesw_profiles[SWBWA_CPE_COUNT];

static void print_ratio_bins(FILE *output, const char *label,
                             const uint64_t bins[SWBWA_MATESW_RATIO_BINS],
                             uint64_t total)
{
    fprintf(output,
            "  %-30s <=1.10=%" PRIu64 " (%.2f%%)"
            " <=1.25=%" PRIu64 " <=1.50=%" PRIu64
            " <=2.00=%" PRIu64 " >2.00=%" PRIu64 "\n",
            label, bins[0], total == 0 ? 0.0 : 100.0 * bins[0] / total,
            bins[1], bins[2], bins[3], bins[4]);
}

void swbwa_cpe_profile_init(void)
{
    evt_conf_t *config = &lwpf_evt_config_SWBWA;

    memset(lwpf_global_counter_SWBWA, 0,
           SWBWA_LWPF_CPE_COUNT * lwpf_kernel_count_SWBWA * NPCR *
               sizeof(*lwpf_global_counter_SWBWA));
    memset(config, 0, sizeof(*config));
    memset(matesw_profiles, 0, sizeof(matesw_profiles));
    config->pc_mask = 0xff;
    config->evt[0] = PC0_CYCLE;
    config->evt[1] = PC1_INST;
    config->evt[2] = PC2_GLD;
    config->evt[3] = PC3_GST;
    config->evt[4] = PC4_DCACHE_ACCESS;
    config->evt[5] = PC5_DCACHE_MISS;
    config->evt[6] = PC6_CYC_MEMB_WAIT;
    config->evt[7] = PC7_CYC_INSTBUFFEREMPTY;
}

swbwa_matesw_profile_t *swbwa_cpe_matesw_profile(void)
{
    return matesw_profiles;
}

long *swbwa_cpe_profile_counters(void)
{
    return lwpf_global_counter_SWBWA;
}

void swbwa_cpe_profile_report(FILE *output)
{
    swbwa_matesw_profile_t total;
    uint64_t pair_ratio_total;
    int i;

    if (output == NULL) output = stderr;
    fprintf(output,
            "\n================ SWBWA CPE Profile (LWPF) ================\n"
            "  sampled CG: %d (64 CPEs)\n"
            "  counters: cycle, instructions, global loads/stores,\n"
            "            D-cache accesses/misses, memory-barrier waits,\n"
            "            instruction-buffer-empty cycles\n"
            "  values: average | minimum | maximum across sampled CPEs\n\n",
            SWBWA_CPE_PROFILE_CG);
    lwpf_report_summary_one(output, "SWBWA",
                            lwpf_global_counter_SWBWA,
                            lwpf_kernel_names_SWBWA,
                            lwpf_kernel_count_SWBWA,
                            &lwpf_evt_config_SWBWA);

    memset(&total, 0, sizeof(total));
    for (i = 0; i < SWBWA_CPE_COUNT; ++i) {
        int j;

        total.invocations += matesw_profiles[i].invocations;
        total.candidates += matesw_profiles[i].candidates;
        total.pairs += matesw_profiles[i].pairs;
        total.paired_candidates += matesw_profiles[i].paired_candidates;
        total.same_orientation_pairs +=
            matesw_profiles[i].same_orientation_pairs;
        total.serial_work += matesw_profiles[i].serial_work;
        total.paired_work += matesw_profiles[i].paired_work;
        for (j = 0; j < SWBWA_MATESW_CANDIDATE_BINS; ++j)
            total.candidate_bins[j] +=
                matesw_profiles[i].candidate_bins[j];
        for (j = 0; j < SWBWA_MATESW_RATIO_BINS; ++j)
            total.ratio_bins[j] += matesw_profiles[i].ratio_bins[j];
        total.sam_pe_calls += matesw_profiles[i].sam_pe_calls;
        for (j = 0; j < SWBWA_MATESW_SAM_CANDIDATE_BINS; ++j)
            total.sam_pe_candidate_bins[j] +=
                matesw_profiles[i].sam_pe_candidate_bins[j];
        total.sam_pe_both_directions +=
            matesw_profiles[i].sam_pe_both_directions;
        total.same_read_pairs += matesw_profiles[i].same_read_pairs;
        total.same_read_paired_candidates +=
            matesw_profiles[i].same_read_paired_candidates;
        total.adjacent_read_groups += matesw_profiles[i].adjacent_read_groups;
        total.adjacent_read_pairs += matesw_profiles[i].adjacent_read_pairs;
        total.adjacent_read_paired_candidates +=
            matesw_profiles[i].adjacent_read_paired_candidates;
        total.sam_u8_candidates += matesw_profiles[i].sam_u8_candidates;
        total.sam_profiled_pairs += matesw_profiles[i].sam_profiled_pairs;
        total.sam_u8_pairs += matesw_profiles[i].sam_u8_pairs;
        total.sam_qlen_equal_pairs +=
            matesw_profiles[i].sam_qlen_equal_pairs;
        total.sam_tlen_equal_pairs +=
            matesw_profiles[i].sam_tlen_equal_pairs;
        total.sam_forward_dimension_equal_pairs +=
            matesw_profiles[i].sam_forward_dimension_equal_pairs;
        total.sam_reverse_dimension_equal_pairs +=
            matesw_profiles[i].sam_reverse_dimension_equal_pairs;
        total.sam_dimension_equal_pairs +=
            matesw_profiles[i].sam_dimension_equal_pairs;
        total.sam_profile_overflow +=
            matesw_profiles[i].sam_profile_overflow;
        for (j = 0; j < SWBWA_MATESW_RATIO_BINS; ++j) {
            total.sam_qlen_ratio_bins[j] +=
                matesw_profiles[i].sam_qlen_ratio_bins[j];
            total.sam_tlen_ratio_bins[j] +=
                matesw_profiles[i].sam_tlen_ratio_bins[j];
            total.sam_forward_work_ratio_bins[j] +=
                matesw_profiles[i].sam_forward_work_ratio_bins[j];
            total.sam_reverse_work_ratio_bins[j] +=
                matesw_profiles[i].sam_reverse_work_ratio_bins[j];
            total.sam_total_work_ratio_bins[j] +=
                matesw_profiles[i].sam_total_work_ratio_bins[j];
        }
        total.sam_total_forward_main_steps +=
            matesw_profiles[i].sam_total_forward_main_steps;
        total.sam_total_forward_lazy_steps +=
            matesw_profiles[i].sam_total_forward_lazy_steps;
        total.sam_total_reverse_main_steps +=
            matesw_profiles[i].sam_total_reverse_main_steps;
        total.sam_total_reverse_lazy_steps +=
            matesw_profiles[i].sam_total_reverse_lazy_steps;
        total.sam_paired_forward_serial_work +=
            matesw_profiles[i].sam_paired_forward_serial_work;
        total.sam_paired_forward_lockstep_work +=
            matesw_profiles[i].sam_paired_forward_lockstep_work;
        total.sam_paired_reverse_serial_work +=
            matesw_profiles[i].sam_paired_reverse_serial_work;
        total.sam_paired_reverse_lockstep_work +=
            matesw_profiles[i].sam_paired_reverse_lockstep_work;
    }

    pair_ratio_total = total.pairs;
    fprintf(output,
            "\nMate-SW dual-alignment feasibility (all %d CPEs)\n"
            "  mem_matesw invocations:       %" PRIu64 "\n"
            "  valid KSW candidates:         %" PRIu64 "\n"
            "  candidate count per call:     0=%" PRIu64
            " 1=%" PRIu64 " 2=%" PRIu64 " 3=%" PRIu64
            " 4=%" PRIu64 "\n"
            "  candidate pairs:              %" PRIu64 "\n"
            "  candidates covered by pairs:  %" PRIu64 " (%.2f%%)\n"
            "  same-orientation pairs:        %" PRIu64 " (%.2f%%)\n"
            "  pair tlen ratio max/min:       <=1.10=%" PRIu64
            " <=1.25=%" PRIu64 " <=1.50=%" PRIu64
            " <=2.00=%" PRIu64 " >2.00=%" PRIu64 "\n"
            "  forward-DP work proxy:         serial=%" PRIu64
            " paired=%" PRIu64 "\n"
            "  ideal paired-kernel speedup:   %.3fx\n"
            "  ideal physical-lane utility:   %.2f%%\n"
            "  note: candidates are paired by descending tlen within each\n"
            "        mem_matesw call; the work proxy excludes reverse-DP\n"
            "        divergence and implementation overhead.\n",
            SWBWA_CPE_COUNT, total.invocations, total.candidates,
            total.candidate_bins[0], total.candidate_bins[1],
            total.candidate_bins[2], total.candidate_bins[3],
            total.candidate_bins[4], total.pairs,
            total.paired_candidates,
            total.candidates == 0 ? 0.0
                                  : 100.0 * total.paired_candidates /
                                        total.candidates,
            total.same_orientation_pairs,
            pair_ratio_total == 0 ? 0.0
                                  : 100.0 * total.same_orientation_pairs /
                                        pair_ratio_total,
            total.ratio_bins[0], total.ratio_bins[1],
            total.ratio_bins[2], total.ratio_bins[3],
            total.ratio_bins[4], total.serial_work, total.paired_work,
            total.paired_work == 0 ? 0.0
                                   : (double)total.serial_work /
                                         total.paired_work,
            total.paired_work == 0 ? 0.0
                                   : 50.0 * total.serial_work /
                                         total.paired_work);
    fprintf(output,
            "\nCross-call dual-alignment opportunities\n"
            "  mem_sam_pe calls:              %" PRIu64 "\n"
            "  KSW candidates per PE:         0=%" PRIu64
            " 1=%" PRIu64 " 2=%" PRIu64 " 3=%" PRIu64
            " 4=%" PRIu64 " 5+=%" PRIu64 "\n"
            "  both mate directions active:   %" PRIu64 " (%.2f%%)\n"
            "  same-PE direction pairs:       %" PRIu64
            " (%" PRIu64 " candidates, %.2f%% coverage)\n"
            "  adjacent-PE groups:            %" PRIu64 "\n"
            "  adjacent-PE pairs:             %" PRIu64
            " (%" PRIu64 " candidates, %.2f%% coverage)\n"
            "  note: same-PE pairs preserve ordering independently in the\n"
            "        two mate directions. Adjacent-PE pairs group consecutive\n"
            "        mem_sam_pe calls executed by each CPE; both are count-only\n"
            "        upper bounds before qlen/tlen compatibility filtering.\n",
            total.sam_pe_calls, total.sam_pe_candidate_bins[0],
            total.sam_pe_candidate_bins[1],
            total.sam_pe_candidate_bins[2],
            total.sam_pe_candidate_bins[3],
            total.sam_pe_candidate_bins[4],
            total.sam_pe_candidate_bins[5],
            total.sam_pe_both_directions,
            total.sam_pe_calls == 0 ? 0.0
                                    : 100.0 * total.sam_pe_both_directions /
                                          total.sam_pe_calls,
            total.same_read_pairs,
            total.same_read_paired_candidates,
            total.candidates == 0 ? 0.0
                                  : 100.0 * total.same_read_paired_candidates /
                                        total.candidates,
            total.adjacent_read_groups, total.adjacent_read_pairs,
            total.adjacent_read_paired_candidates,
            total.candidates == 0 ? 0.0
                                  : 100.0 *
                                        total.adjacent_read_paired_candidates /
                                        total.candidates);
    {
        uint64_t forward_work = total.sam_total_forward_main_steps +
                                total.sam_total_forward_lazy_steps;
        uint64_t reverse_work = total.sam_total_reverse_main_steps +
                                total.sam_total_reverse_lazy_steps;
        uint64_t paired_serial = total.sam_paired_forward_serial_work +
                                 total.sam_paired_reverse_serial_work;
        uint64_t paired_lockstep = total.sam_paired_forward_lockstep_work +
                                   total.sam_paired_reverse_lockstep_work;
        uint64_t projected_work = forward_work + reverse_work -
                                  paired_serial + paired_lockstep;

        fprintf(output,
                "\nSame-PE length and measured DP-work compatibility\n"
                "  profiled ordered pairs:        %" PRIu64
                " (%.2f%% of count-based pairs)\n"
                "  u8 candidates / pairs:         %" PRIu64
                " / %" PRIu64 "\n"
                "  exact forward qlen pairs:      %" PRIu64 " (%.2f%%)\n"
                "  exact forward tlen pairs:      %" PRIu64 " (%.2f%%)\n"
                "  exact forward dimensions:      %" PRIu64 " (%.2f%%)\n"
                "  exact reverse dimensions:      %" PRIu64 " (%.2f%%)\n"
                "  exact both-pass dimensions:    %" PRIu64 " (%.2f%%)\n"
                "  dropped profile candidates:    %" PRIu64 "\n",
                total.sam_profiled_pairs,
                total.same_read_pairs == 0 ? 0.0
                    : 100.0 * total.sam_profiled_pairs /
                      total.same_read_pairs,
                total.sam_u8_candidates, total.sam_u8_pairs,
                total.sam_qlen_equal_pairs,
                total.sam_profiled_pairs == 0 ? 0.0
                    : 100.0 * total.sam_qlen_equal_pairs /
                      total.sam_profiled_pairs,
                total.sam_tlen_equal_pairs,
                total.sam_profiled_pairs == 0 ? 0.0
                    : 100.0 * total.sam_tlen_equal_pairs /
                      total.sam_profiled_pairs,
                total.sam_forward_dimension_equal_pairs,
                total.sam_profiled_pairs == 0 ? 0.0
                    : 100.0 * total.sam_forward_dimension_equal_pairs /
                      total.sam_profiled_pairs,
                total.sam_reverse_dimension_equal_pairs,
                total.sam_profiled_pairs == 0 ? 0.0
                    : 100.0 * total.sam_reverse_dimension_equal_pairs /
                      total.sam_profiled_pairs,
                total.sam_dimension_equal_pairs,
                total.sam_profiled_pairs == 0 ? 0.0
                    : 100.0 * total.sam_dimension_equal_pairs /
                      total.sam_profiled_pairs,
                total.sam_profile_overflow);
        print_ratio_bins(output, "qlen ratio max/min:",
                         total.sam_qlen_ratio_bins,
                         total.sam_profiled_pairs);
        print_ratio_bins(output, "tlen ratio max/min:",
                         total.sam_tlen_ratio_bins,
                         total.sam_profiled_pairs);
        print_ratio_bins(output, "forward work ratio max/min:",
                         total.sam_forward_work_ratio_bins,
                         total.sam_profiled_pairs);
        print_ratio_bins(output, "reverse work ratio max/min:",
                         total.sam_reverse_work_ratio_bins,
                         total.sam_profiled_pairs);
        print_ratio_bins(output, "total DP work ratio max/min:",
                         total.sam_total_work_ratio_bins,
                         total.sam_profiled_pairs);
        fprintf(output,
                "\n  Measured vector-step work (main + Lazy-F)\n"
                "    all forward:                 %" PRIu64
                " + %" PRIu64 " (Lazy-F %.2f%%)\n"
                "    all reverse:                 %" PRIu64
                " + %" PRIu64 " (Lazy-F %.2f%%)\n"
                "    paired forward:              serial=%" PRIu64
                " lockstep=%" PRIu64 " (%.3fx)\n"
                "    paired reverse:              serial=%" PRIu64
                " lockstep=%" PRIu64 " (%.3fx)\n"
                "    all candidate DP projected:  original=%" PRIu64
                " lockstep=%" PRIu64 " (%.3fx)\n"
                "  note: pairs are matched in original call order across the\n"
                "        two independent mate directions. Work counts striped\n"
                "        main-loop and actual Lazy-F vector steps; it excludes\n"
                "        query-profile construction and implementation overhead.\n",
                total.sam_total_forward_main_steps,
                total.sam_total_forward_lazy_steps,
                forward_work == 0 ? 0.0
                    : 100.0 * total.sam_total_forward_lazy_steps /
                      forward_work,
                total.sam_total_reverse_main_steps,
                total.sam_total_reverse_lazy_steps,
                reverse_work == 0 ? 0.0
                    : 100.0 * total.sam_total_reverse_lazy_steps /
                      reverse_work,
                total.sam_paired_forward_serial_work,
                total.sam_paired_forward_lockstep_work,
                total.sam_paired_forward_lockstep_work == 0 ? 0.0
                    : (double)total.sam_paired_forward_serial_work /
                      total.sam_paired_forward_lockstep_work,
                total.sam_paired_reverse_serial_work,
                total.sam_paired_reverse_lockstep_work,
                total.sam_paired_reverse_lockstep_work == 0 ? 0.0
                    : (double)total.sam_paired_reverse_serial_work /
                      total.sam_paired_reverse_lockstep_work,
                forward_work + reverse_work, projected_work,
                projected_work == 0 ? 0.0
                    : (double)(forward_work + reverse_work) /
                      projected_work);
    }
    fputs("=============================================================\n",
          output);
}

#else

void swbwa_cpe_profile_init(void)
{
}

long *swbwa_cpe_profile_counters(void)
{
    return NULL;
}

swbwa_matesw_profile_t *swbwa_cpe_matesw_profile(void)
{
    return NULL;
}

void swbwa_cpe_profile_report(FILE *output)
{
    (void)output;
}

#endif
