#ifndef M_B_FUNC_H
#define M_B_FUNC_H

#include "parameters.h"
#include <gmp.h>
#include <math.h>
#include <mpfr.h>

float calculate_log_thresold(const mpz_t n,long M);
void calculate_news_M_and_B(long*M,long*B);
void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B,long*num_elem_array_number);
#endif