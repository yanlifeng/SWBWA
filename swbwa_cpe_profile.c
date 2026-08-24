#include <string.h>

#include "swbwa_config.h"
#include "swbwa_cpe_profile.h"

#if SWBWA_ENABLE_CPE_PROFILE

#define LWPF_NOCOLOR
#define LWPF_UNITS U(SWBWA)
#include "lwpf.h"

enum { SWBWA_LWPF_CPE_COUNT = 64 };

void swbwa_cpe_profile_init(void)
{
    evt_conf_t *config = &lwpf_evt_config_SWBWA;

    memset(lwpf_global_counter_SWBWA, 0,
           SWBWA_LWPF_CPE_COUNT * lwpf_kernel_count_SWBWA * NPCR *
               sizeof(*lwpf_global_counter_SWBWA));
    memset(config, 0, sizeof(*config));
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

long *swbwa_cpe_profile_counters(void)
{
    return lwpf_global_counter_SWBWA;
}

void swbwa_cpe_profile_report(FILE *output)
{
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

void swbwa_cpe_profile_report(FILE *output)
{
    (void)output;
}

#endif
