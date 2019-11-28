
#ifndef CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H
#define CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H

#include "mpz_functions.h"
#include "properties.h"
#include <mpfr.h>
#include "factor_base_functions.h"
#include "matrix_function_v2.h"
//#include "square_relations_functions.h"

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
void calculate_news_M_and_B(long*M,long*B);
void adjust_array_bi(mpz_t *array_Bk,int s,const mpz_t a);
mpz_t*calculate_array_Bk_f(int*number_prime_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1);
void increment_M_and_B(long*M,long*B);
void calculate_thresold_a(mpfr_t thresold_a,const mpz_t radn,long M);
void calculate_x0(mpz_t x0,const mpz_t n,int k,char*factorized);
void multiply_n_for_k(mpz_t n,int*k,char*factorized);
void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B);
mpz_t*calculate_bi(mpz_t *array_Bk,const mpz_t b1,int s);
void calculate_a_f2(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,int**best_q,int**best_q_number);
//int find_factor_of_n_from_base_matrix(int **base_matrix,int num_row,int* num_column,int **matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base);


/*void calculate_a_and_b(int*solution,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base,mpz_t a,mpz_t b,const mpz_t n);
int try_to_factor(const mpz_t a,const mpz_t b,const mpz_t n,mpz_t factor1,mpz_t factor2);
int factor_n(int **matrix_solution,int dim_sol,int num_B_smooth,mpz_t**matrix_B_smooth,int card_f_base,const mpz_t n);
int**generate_all_solution(int **base_matrix,int num_row,int* num_column,int**matrix_linear_system,int num_row_matrix,int num_col_matrix);
void calculate_best_M(const mpz_t n,long*M);
void calculate_best_B(const mpz_t n,long*B);
void calculate_thresold_a(mpfr_t thresold_a,const mpz_t radn,long M);
struct node* create_factor_base(int*cardinality_factor_base,long B,struct node**tail,const mpz_t n,mpz_t q,const mpz_t thresold_q);
long calculate_best_s();
void calculate_a(mpz_t a,const mpfr_t thresold_a,int* s,struct node*head_f_base,long cardinality_factor_base,long**array_of_prime_chosen_for_a);
mpz_t *calculate_array_Bk(long*array_of_prime_chosen_for_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1);
void multiply_n_for_k(mpz_t n,int*k,char*factorized);
mpz_t*calculate_bi(mpz_t *array_Bk,const mpz_t b1,int s);
mpz_t *create_array_of_number(const mpz_t a,const mpz_t b,long M,const mpz_t n);
void calculate_x0(mpz_t x0,const mpz_t n,int k,char*factorized);
void calculate_thresold_q(mpz_t thresold_q,const mpq_t thresold_a);
void calculate_b_mpqs(mpz_t b_mpqs,const mpz_t thresold_q,const mpz_t n,const mpz_t a);
void add_remainder_to_matrix_factorization(mpz_t **matrix_factorization,mpz_t*array_number,long M,int cardinality_factor_base);
void adjust_array_bi(mpz_t *array_Bk,int s,const mpz_t a);
void calculate_news_M_and_B(long*M,long*B);
void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B);
int find_factor_of_n_from_base_matrix(int **base_matrix,int num_row,int* num_column,int **matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base);
void create_factor_base_f(int*cardinality_factor_base,long B,struct node_factor_base**head,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base);
void calculate_a_f(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,
long**array_of_prime_chosen_for_a,long**Q,long**Q_number);
void calculate_target_a1_f(mpfr_t target_a1,const mpfr_t target_a,struct node_factor_base*head_f_base_f,
long p_min_i,long p_max_i,int 	cardinality_factor_base);
void calculate_p_min_p_max_i_f(long*p_min_i,long*p_max_i,struct node_factor_base*head_f_base_f,long cardinality_factor_base);
void calculate_a_f2(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,int**best_q,int**best_q_number);
mpz_t*calculate_array_Bk_f(int*number_prime_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1);
int find_factor_of_n_from_base_matrix_char(int **base_matrix,int num_row,int* num_column,char*matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,struct node_square_relation*head,int num_B_smooth,int card_f_base);
struct node_factor_base*initialize_factor_base(int*cardinality_factor_base,long B,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base);
void increment_M_and_B(long*M,long*B);*/
#endif //CRIVELLO_QUADRATICO__COPIA_1_A_B_C_FUNCTIONS_H
