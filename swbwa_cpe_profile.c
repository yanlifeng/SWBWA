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
