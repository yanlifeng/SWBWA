#ifndef SWBWA_CPE_PROFILE_H
#define SWBWA_CPE_PROFILE_H

#include <stdio.h>

void swbwa_cpe_profile_init(void);
long *swbwa_cpe_profile_counters(void);
void swbwa_cpe_profile_report(FILE *output);

#endif
