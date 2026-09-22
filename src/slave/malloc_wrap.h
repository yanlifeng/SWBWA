#ifndef SWBWA_CPE_MALLOC_WRAP_H
#define SWBWA_CPE_MALLOC_WRAP_H

#include <stdlib.h>
#include <string.h>
#include "swbwa_ldm_alloc.h"

#ifdef __cplusplus
extern "C" {
#endif

void swbwa_cpe_malloc_stats_init(void);
void swbwa_cpe_malloc_stats_print(void);
void set_big_buffer(char *buffer, long long bytes_per_cpe);

void *wrap_calloc(size_t nmemb, size_t size,
                  const char *file, unsigned int line, const char *func);
void *wrap_malloc(size_t size,
                  const char *file, unsigned int line, const char *func);
void *wrap_realloc(void *ptr, size_t size,
                   const char *file, unsigned int line, const char *func);
char *wrap_strdup(const char *s,
                  const char *file, unsigned int line, const char *func);
void wrap_free(void *ptr,
               const char *file, unsigned int line, const char *func);

#ifdef __cplusplus
}
#endif

#if SWBWA_ENABLE_CPE_MALLOC_WRAPPER && !defined(SWBWA_ALLOC_IMPLEMENTATION)
#ifdef calloc
#undef calloc
#endif
#define calloc(n, s) wrap_calloc((n), (s), __FILE__, __LINE__, __func__)

#ifdef malloc
#undef malloc
#endif
#define malloc(s) wrap_malloc((s), __FILE__, __LINE__, __func__)

#ifdef realloc
#undef realloc
#endif
#define realloc(p, s) wrap_realloc((p), (s), __FILE__, __LINE__, __func__)

#ifdef strdup
#undef strdup
#endif
#define strdup(s) wrap_strdup((s), __FILE__, __LINE__, __func__)

#ifdef free
#undef free
#endif
#define free(p) wrap_free((p), __FILE__, __LINE__, __func__)
#if SWBWA_CPE_LDM_ALLOC
#undef malloc
#undef calloc
#undef realloc
#undef free
#if SWBWA_CPE_LDM_ALLOC == 2 || SWBWA_CPE_LDM_ALLOC == 4
/* Unknown fresh objects cannot belong to LDM. Fold this branch at the call
 * site; only profile mode pays to count the ineligible allocations. */
#define malloc(s) (swbwa_ldm_alloc_site(__func__) == SWBWA_LDM_SITE_OTHER \
    ? wrap_malloc((s), __FILE__, __LINE__, __func__) \
    : swbwa_auto_malloc((s), swbwa_ldm_alloc_site(__func__)))
#define calloc(n, s) (swbwa_ldm_alloc_site(__func__) == SWBWA_LDM_SITE_OTHER \
    ? wrap_calloc((n), (s), __FILE__, __LINE__, __func__) \
    : swbwa_auto_calloc((n), (s), swbwa_ldm_alloc_site(__func__)))
#else
#define malloc(s) swbwa_auto_malloc((s), swbwa_ldm_alloc_site(__func__))
#define calloc(n, s) swbwa_auto_calloc((n), (s), swbwa_ldm_alloc_site(__func__))
#endif
#define realloc(p, s) swbwa_auto_realloc((p), (s), swbwa_ldm_alloc_site(__func__))
#define free(p) swbwa_auto_free(p)
#if SWBWA_CPE_LDM_ALLOC == 3
#undef strdup
#define strdup(s) swbwa_auto_strdup((s), swbwa_ldm_alloc_site(__func__))
#endif
#endif
#endif

#endif /* SWBWA_CPE_MALLOC_WRAP_H */
