#include <stdlib.h>
#include <string.h>

#ifdef HOST_USE_MALLOC_WRAPPERS
#  undef HOST_USE_MALLOC_WRAPPERS
#endif

#include "malloc_wrap.h"

void htmalloccount_init(void)
{
}

void htmalloccount_print(void)
{
}

void *wrap_malloc(size_t size, const char *file, unsigned int line, const char *func)
{
    (void)file;
    (void)line;
    (void)func;
    return malloc(size);
}

void wrap_free(void *ptr, const char *file, unsigned int line, const char *func)
{
    (void)file;
    (void)line;
    (void)func;
    free(ptr);
}

void *wrap_realloc(void *ptr, size_t size, const char *file, unsigned int line, const char *func)
{
    (void)file;
    (void)line;
    (void)func;
    return realloc(ptr, size);
}

void *wrap_calloc(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func)
{
    size_t bytes;
    void *ptr;

    (void)file;
    (void)line;
    (void)func;

    if (size != 0 && nmemb > (size_t)-1 / size)
        return NULL;

    bytes = nmemb * size;
    ptr = malloc(bytes);
    if (ptr != NULL)
        memset(ptr, 0, bytes);
    return ptr;
}

char *wrap_strdup(const char *s, const char *file, unsigned int line, const char *func)
{
    size_t len;
    char *copy;

    (void)file;
    (void)line;
    (void)func;

    if (s == NULL)
        return NULL;

    len = strlen(s) + 1;
    copy = malloc(len);
    if (copy != NULL)
        memcpy(copy, s, len);
    return copy;
}
