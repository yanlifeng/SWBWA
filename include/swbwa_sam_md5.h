#ifndef SWBWA_SAM_MD5_H
#define SWBWA_SAM_MD5_H
#include <stddef.h>
#include <stdio.h>
#ifdef __cplusplus
extern "C" {
#endif
/* One ordered stream, one writer; chunks may split any header/record boundary. */
int swbwa_sam_md5_open(void);
int swbwa_sam_md5_update(const void *data, size_t size);
int swbwa_sam_md5_close(FILE *report);
void swbwa_sam_md5_abort(void);
#ifdef __cplusplus
}
#endif
#endif
