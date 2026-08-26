#ifndef SWBWA_SLAVE_CPE_PROFILE_H
#define SWBWA_SLAVE_CPE_PROFILE_H

#include "../swbwa_config.h"

enum swbwa_cpe_profile_region {
    SWBWA_CPE_PROFILE_WORKER_ALIGNMENT,
    SWBWA_CPE_PROFILE_MEM_CHAIN,
    SWBWA_CPE_PROFILE_CHAIN_FILTER,
    SWBWA_CPE_PROFILE_CHAIN_EXTENSION,
    SWBWA_CPE_PROFILE_ALIGNMENT_FINALIZE,
    SWBWA_CPE_PROFILE_MATE_RESCUE,
    SWBWA_CPE_PROFILE_PAIRING,
    SWBWA_CPE_PROFILE_SAM_FORMAT,
    SWBWA_CPE_PROFILE_SAM_COPY,
    SWBWA_CPE_PROFILE_MEM_CHAIN_COLLECT,
    SWBWA_CPE_PROFILE_MEM_CHAIN_BUILD,
    SWBWA_CPE_PROFILE_CHAIN_EXTENSION_DP,
    SWBWA_CPE_PROFILE_MATE_REF_FETCH,
    SWBWA_CPE_PROFILE_MATE_KSW_ALIGN,
    SWBWA_CPE_PROFILE_KSW_QUERY_INIT_FORWARD,
    SWBWA_CPE_PROFILE_KSW_DP_FORWARD,
    SWBWA_CPE_PROFILE_KSW_QUERY_INIT_REVERSE,
    SWBWA_CPE_PROFILE_KSW_DP_REVERSE,
    SWBWA_CPE_PROFILE_MATE_DEDUP,
    SWBWA_CPE_PROFILE_DEDUP_SORT_END,
    SWBWA_CPE_PROFILE_DEDUP_REDUNDANCY,
    SWBWA_CPE_PROFILE_DEDUP_SORT_SCORE
};

#if SWBWA_ENABLE_CPE_PROFILE

#include <slave.h>

#define EVT_PC0 PC0_CYCLE
#define EVT_PC1 PC1_INST
#define EVT_PC2 PC2_GLD
#define EVT_PC3 PC3_GST
#define EVT_PC4 PC4_DCACHE_ACCESS
#define EVT_PC5 PC5_DCACHE_MISS
#define EVT_PC6 PC6_CYC_MEMB_WAIT
#define EVT_PC7 PC7_CYC_INSTBUFFEREMPTY

#define LWPF_KERNELS \
    K(WORKER_ALIGNMENT) \
    K(MEM_CHAIN) \
    K(CHAIN_FILTER) \
    K(CHAIN_EXTENSION) \
    K(ALIGNMENT_FINALIZE) \
    K(MATE_RESCUE) \
    K(PAIRING) \
    K(SAM_FORMAT) \
    K(SAM_COPY) \
    K(MEM_CHAIN_COLLECT) \
    K(MEM_CHAIN_BUILD) \
    K(CHAIN_EXTENSION_DP) \
    K(MATE_REF_FETCH) \
    K(MATE_KSW_ALIGN) \
    K(KSW_QUERY_INIT_FORWARD) \
    K(KSW_DP_FORWARD) \
    K(KSW_QUERY_INIT_REVERSE) \
    K(KSW_DP_REVERSE) \
    K(MATE_DEDUP) \
    K(DEDUP_SORT_END) \
    K(DEDUP_REDUNDANCY) \
    K(DEDUP_SORT_SCORE)
#define LWPF_UNIT U(SWBWA)
#include "lwpf.h"

static inline int swbwa_cpe_profile_active(void)
{
    return (int)_MYID / 64 == SWBWA_CPE_PROFILE_CG;
}

static inline void swbwa_cpe_profile_enter(long *counters)
{
    evt_conf_t local_config;
    volatile int reply = 0;
    int local_id;

    if (!swbwa_cpe_profile_active()) return;
    local_id = (int)_MYID & 63;
    lwpf_sync_counters_m2c(
        counters + local_id * lwpf_kernel_count_SWBWA * NPCR,
                           lwpf_kernel_count_SWBWA);
    athread_get(PE_MODE, &lwpf_evt_config_SWBWA, &local_config,
                sizeof(local_config), (void *)&reply, 0, 0, 0);
    while (reply != 1) { }
    config_pcrs(local_config.evt);
}

static inline void swbwa_cpe_profile_exit(long *counters)
{
    int local_id;

    if (!swbwa_cpe_profile_active()) return;
    local_id = (int)_MYID & 63;
    lwpf_sync_counters_c2m(
        counters + local_id * lwpf_kernel_count_SWBWA * NPCR,
                           lwpf_kernel_count_SWBWA);
}

static inline void swbwa_cpe_profile_start(
    enum swbwa_cpe_profile_region region)
{
    if (swbwa_cpe_profile_active()) lwpf_start(region);
}

static inline void swbwa_cpe_profile_stop(
    enum swbwa_cpe_profile_region region)
{
    if (swbwa_cpe_profile_active()) lwpf_stop(region);
}

#else

static inline void swbwa_cpe_profile_enter(long *counters)
{
    (void)counters;
}

static inline void swbwa_cpe_profile_exit(long *counters)
{
    (void)counters;
}

static inline void swbwa_cpe_profile_start(
    enum swbwa_cpe_profile_region region)
{
    (void)region;
}

static inline void swbwa_cpe_profile_stop(
    enum swbwa_cpe_profile_region region)
{
    (void)region;
}

#endif

#endif
