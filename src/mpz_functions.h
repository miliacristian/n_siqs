#ifndef MPZ_H
#define MPZ_H
#include "printv2.h"
#include "dynamic_list.h"

void print_array_mpz(mpz_t*array,int length);

void print_list_mpz(struct node*head);
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col);
void free_memory_array_mpz(mpz_t*array,long length);
int quadratic_residue(mpz_t x,const mpz_t q,const mpz_t n);
#endif
