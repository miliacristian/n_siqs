#include "dynamic_list.h"
#include "list_factor_base.h"
#include "basic.h"
#include <gmp.h>
void print_array_char(char*array,int lenght);
void print_array_long(long*array,int lenght);
void print_array_int(int*array,int lenght);
void print_list(struct node*head);
void print_vector_char(char*array,int lenght);
void print_vector_long(long*array,int lenght);
void print_vector_int(int*array,int lenght);
void print_matrix_int(int**matrix,int num_row,int num_col);
 void print_matrix_factorization(mpz_t**matrix_factorization,int M,int cardinality_factor_base);
void print_matrix_B_smooth(mpz_t**matrix_B_smooth,int num_of_B_smooth,int cardinality_factor_base);
void print_linear_system(int**linear_system,int cardinality_factor_base,int num_of_B_smooth);
void print_base_linear_system(int**linear_system,int num_B_smooth,int dim_sol);
void print_all_solution(int**matrix_solution,int dim_sol,int num_B_smooth);
void print_factor_base(struct node*head_f_base);
void print_array_number(mpz_t*array,int lenght);
void print_array_double(double*array,int lenght);
void print_matrix_double(double**matrix,int num_row,int num_col);
void print_matrix_long(long**matrix,int num_row,int num_col);
void print_int(int a,char*string);
void print_char(char a,char*string);
void print_double(double a,char*string);
void print_long(long a,char*string);
void print_array_Bk(mpz_t*array_Bk,long s);
void print_array_mpz(mpz_t*array,int length);
void print_array_bi(mpz_t* array_bi,long s);
void print_index_B_smooth(struct node*head);
void print_array_chosen_for_a(long*array_of_prime_chosen_for_a,int card_factor_base,int s);
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col);
void print_array_matrix_same_dimension(mpz_t***array_matrix_mpz,int length_array,int num_row,int num_col);
void print_matrix_factorization_f(struct matrix_factorization m);
void print_array_float(float*array,int length);
void print_list_f(struct node_factor_base*head);
