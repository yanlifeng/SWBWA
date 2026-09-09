#ifndef SWBWA_CPE_PROFILE_H
#define SWBWA_CPE_PROFILE_H

#include <stdio.h>

#include "swbwa_matesw_profile.h"

void swbwa_cpe_profile_init(void);
long *swbwa_cpe_profile_counters(void);
swbwa_matesw_profile_t *swbwa_cpe_matesw_profile(void);
void swbwa_cpe_profile_report(FILE *output);

#endif
