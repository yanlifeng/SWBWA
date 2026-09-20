#ifndef SWBWA_TEST_SLAVE_H
#define SWBWA_TEST_SLAVE_H
extern int swbwa_test_cpe_id;
#define _MYID swbwa_test_cpe_id
void *ldm_malloc(unsigned long bytes);
void ldm_free(void *ptr, unsigned long bytes);
#endif
