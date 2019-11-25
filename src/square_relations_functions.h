#ifndef Y_H
#define Y_H
#include "factorization_functionsv2.h"
#include <gmp.h>
struct square_relation {
    mpz_t square;
    struct node_factorization*head_factorization;
    mpz_t num;
    mpz_t residuos;
};
struct node_square_relation {
    struct square_relation square_relation;
    struct node_square_relation*next;
    struct node_square_relation*prev;
};
#endif
