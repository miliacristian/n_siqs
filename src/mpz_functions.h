#ifndef MPZ_H
#define MPZ_H
#include "print.h"
#include "dynamic_list.h"

void print_array_mpz(mpz_t*array,int length);
void print_list_mpz(struct node*head);
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col);
void free_memory_array_mpz(mpz_t*array,long length);
int quadratic_residue(mpz_t x,const mpz_t q,const mpz_t n);
mpz_t*alloc_array_mpz(int length);
void reduce_array_mpz_mod_n(mpz_t*array,int length,const mpz_t a);
void sum_elem_multiple_of_2_mpz(mpz_t*vector1,mpz_t*vector2,int length1,int length2);
void divide_elem_multiple_of_2_by_x(mpz_t*vector,int length,double x);
#endif
