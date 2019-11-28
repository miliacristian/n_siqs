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
char verify_square_relation(struct square_relation square_relation,const mpz_t n);
void calculate_a_and_b_siqs(const int*solution,struct node_square_relation*head,int num_B_smooth,int card_f_base,mpz_t a,mpz_t b,const mpz_t n);
int find_factor_of_n_from_base_matrix_char(int **base_matrix,int num_row,int* num_column,char*matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,struct node_square_relation*head,int num_B_smooth,int card_f_base);
#endif
