#ifndef SWBWA_HOST_MALLOC_WRAP_H
#define SWBWA_HOST_MALLOC_WRAP_H

#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

void swbwa_host_malloc_stats_init(void);
void swbwa_host_malloc_stats_print(void);

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

#if SWBWA_ENABLE_HOST_MALLOC_WRAPPER
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
#endif

#endif /* SWBWA_HOST_MALLOC_WRAP_H */
