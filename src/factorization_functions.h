#ifndef FACTOR_H
#define FACTOR_H

#include "properties.h"
#include "mpz_functions.h"
#include "factorization_functionsv2.h"
#include "square_relations_functions.h"
#include "a_b_c_BK_functions.h"
#include "factor_base_functions.h"
void print_factorization(const mpz_t num,struct node_factorization*head_factor);
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s);
#endif
