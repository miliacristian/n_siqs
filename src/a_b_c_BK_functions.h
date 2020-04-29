
#ifndef A_B_C_FUNCTIONS_H
#define A_B_C_FUNCTIONS_H

#include "mpz_functions.h"
#include <mpfr.h>
#include "factor_base_functions.h"
#include "matrix_functions.h"

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
int find_factor_of_n_from_base_matrix(int **base_matrix,int num_row,int* num_column,int **matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base);
void adjust_array_bi(mpz_t *array_Bk,int s,const mpz_t a);
mpz_t* calculate_array_Bk(int*number_prime_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1);
void increment_M_and_B(long*M,long*B);
void calculate_thresold_a(mpfr_t thresold_a,const mpz_t radn,long M);
void calculate_x0(mpz_t x0,const mpz_t n,int k,char*factorized);
void multiply_n_for_k(mpz_t n,int*k,char*factorized);
//void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B);
mpz_t* calculate_bi(mpz_t *array_Bk,const mpz_t b1,int s);
void calculate_a(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,int**best_q,int**best_q_number);
void calculate_thresold_large_prime(mpz_t thresold_large_prime,int max_prime);
struct a_struct*create_array_a_struct(int*number_prime_a,int*index_number_a,int length);
#endif //CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H
