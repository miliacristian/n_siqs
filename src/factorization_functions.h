#ifndef FACTOR_H
#define FACTOR_H

#include "properties.h"
#include "mpz_functions.h"
#include "list_square_relation.h"
struct node_factorization {
    int number;
    int exp_of_number;
    int index;
    struct node_factorization*next;
    struct node_factorization*prev;
};
void print_factorization(const mpz_t num,struct node_factorization*head_factor);
#endif
