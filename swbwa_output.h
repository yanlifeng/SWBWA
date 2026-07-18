#ifndef SWBWA_OUTPUT_H
#define SWBWA_OUTPUT_H

#include <stddef.h>

int swbwa_output_open(const char *path);
int swbwa_output_write(const void *data, size_t length);
int swbwa_output_flush(void);
int swbwa_output_close(void);
const char *swbwa_output_name(void);

#endif /* SWBWA_OUTPUT_H */
