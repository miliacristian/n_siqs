#ifndef Y_H
#define Y_H
#include "factorization_functionsv2.h"
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "properties.h"
#include "matrix_function_v2.h"
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
void free_memory_list_square_relation(struct node_square_relation*head);
unsigned long**create_binary_linear_system(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth,int*num_col_binary_matrix);
#endif
