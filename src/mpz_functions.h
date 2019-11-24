#ifndef MPZ_H
#define MPZ_H
#include "print.h"

void print_array_mpz(mpz_t*array,int length);

void print_list_mpz(struct node*head);
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col);
#endif
