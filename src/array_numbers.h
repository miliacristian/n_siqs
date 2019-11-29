
#ifndef ARRAY_NUMBERS_H
#define ARRAY_NUMBERS_H

#include "print.h"
#include "mpz_functions.h"

struct number {
    int j;//indice j va da -M a M
    int sum_log;//somma del logaritmo
    int first_index_f_base;//primo indice del primo della factor base che divide number,-1 non è considerato
    int last_index_f_base;//ultimo indice del primo della factor base che divide number,-1 non è considerato
    //number=a^2j^2+2*a*b*j+b^2-n,b^2-n=a*c
};
#endif //CRIVELLO_QUADRATICO__COPIA_1_ARRAY_NUMBERS_H
