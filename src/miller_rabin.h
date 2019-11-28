#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <gmp.h>
#include <stdio.h>
#include "parameters.h"

char get_and_check_n(mpz_t n,FILE*file_number);
void test_n(const mpz_t n,int num_test);

#endif
