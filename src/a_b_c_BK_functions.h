
#ifndef CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H
#define CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H

#include "mpz_functions.h"
struct a_struct {
    int number_prime_a;
    int index_prime_a;
};
int compare_a_struct( const void* a, const void* b);
void print_array_a_struct(struct a_struct*array_a_struct,int length);
char value_is_in_sorted_array(int index_of_prime,struct a_struct*array_a_struct,int length);
void calculate_root_poly_second_degree_mod_p_to_k(mpz_t j1t,const mpz_t p_to_k,long p,long k,const mpz_t a,const mpz_t b,const mpz_t n);
void create_num(mpz_t num,const mpz_t a,const mpz_t b,const mpz_t n,long j);
void calculate_square(mpz_t square,const mpz_t a,int index,const mpz_t b,const mpz_t n);
#endif //CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H
