#ifndef SWBWA_CONFIG_H
#define SWBWA_CONFIG_H

/*
 * Build-time configuration shared by the MPE and CPE code.
 *
 * Override these from the compiler command line when needed, for example:
 *   -DSWBWA_USE_CGS_MODE=0 -DSWBWA_ENABLE_LWPF3=1
 *
 * The lower-case compatibility macros are kept so the existing code can be
 * migrated gradually instead of carrying scattered local #defines forever.
 */

#ifndef SWBWA_USE_CGS_MODE
#define SWBWA_USE_CGS_MODE 1
#endif

#if SWBWA_USE_CGS_MODE
#define SWBWA_CPE_NUM 384
#define SWBWA_CG_NUM 6
#ifndef use_cgs_mode
#define use_cgs_mode 1
#endif
#else
#define SWBWA_CPE_NUM 64
#define SWBWA_CG_NUM 1
#endif

#ifndef SWBWA_ENABLE_CPE_STEP0
#define SWBWA_ENABLE_CPE_STEP0 0
#endif
#if SWBWA_ENABLE_CPE_STEP0
#ifndef use_cpe_step0
#define use_cpe_step0 1
#endif
#endif

#ifndef SWBWA_ENABLE_MY_MPI
#define SWBWA_ENABLE_MY_MPI 0
#endif
#if SWBWA_ENABLE_MY_MPI
#ifndef use_my_mpi
#define use_my_mpi 1
#endif
#endif

#ifndef SWBWA_ENABLE_LWPF3
#define SWBWA_ENABLE_LWPF3 0
#endif
#if SWBWA_ENABLE_LWPF3
#ifndef use_lwpf3
#define use_lwpf3 1
#endif
#endif

#ifndef SWBWA_PRE_CPE_READ_NUM
#define SWBWA_PRE_CPE_READ_NUM 1024
#endif

#ifndef SWBWA_EVAL_READ_SIZE
#define SWBWA_EVAL_READ_SIZE 300
#endif

#ifndef SWBWA_CPE_MALLOC_PER_CPE
#define SWBWA_CPE_MALLOC_PER_CPE (24LL << 20)
#endif

#define SWBWA_CPE_MALLOC_TOTAL_SIZE (1LL * SWBWA_CPE_NUM * SWBWA_CPE_MALLOC_PER_CPE)

#ifndef cpe_num
#define cpe_num SWBWA_CPE_NUM
#endif

#ifndef cg_num
#define cg_num SWBWA_CG_NUM
#endif

#ifndef cpe_num_slave
#define cpe_num_slave SWBWA_CPE_NUM
#endif

#endif
