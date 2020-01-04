
#ifndef ARRAY_NUMBERS_H
#define ARRAY_NUMBERS_H

#include "print.h"
#include "mpz_functions.h"

struct number {
	int sum_log;//somma del logaritmo
    int j;//indice j va da -M a M
    //number=a^2j^2+2*a*b*j+b^2-n,b^2-n=a*c
};
#define OFFSET_FROM_STRUCT_NUMBER_TO_SUM_LOG offsetof(struct number,sum_log)
#define NUM_INT_IN_STRUCT_NUMBER sizeof(struct number)/sizeof(int)
#endif //CRIVELLO_QUADRATICO__COPIA_1_ARRAY_NUMBERS_H
