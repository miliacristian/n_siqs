#ifndef M_B_FUNC_H
#define M_B_FUNC_H

#include "parameters.h"
#include <gmp.h>
#include <math.h>
#include <mpfr.h>

float calculate_log_thresold(const mpz_t n,long M);
void calculate_news_M_and_B(long*M,long*B);
#endif