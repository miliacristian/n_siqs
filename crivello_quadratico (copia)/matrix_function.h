#include "math_function.h"
#include "list_factor_base.h"
#ifndef O_H
#define O_H

int* sum_vector(int*vector1,int*vector2,int lenght1,int lenght2);
long* sum_elem_multiple_of_2(long*vector1,long*vector2,int length1,int length2);
void sum_elem_multiple_of_2_mpz(mpz_t*vector1,mpz_t*vector2,int length1,int length2);
void prod_vector_and_scalar(int*vector,int scalar,int lenght);
int* prod_vector_and_scalar_v2(int*vector,int scalar,int lenght);
int*get_col(int matrix[7][3],int num_row,int num_col,int index_col);
void reduce_matrix_mod_n(int**linear_system,int cardinality_factor_base,int num_of_B_smooth,int n);
mpz_t*create_array_temp_factorization(int card_f_base,mpz_t*row_matrix_b_smooth);
char is_in_array_long(long*array,long length,long p_i);
mpz_t**create_matrix_factorization(long M,int cardinality_factor_base,struct node*head_f_base,long*array_of_prime_chosen_for_a,const mpz_t a,const mpz_t b);
int**alloc_linear_system(int cardinality_f_base,int num_B_smooth);//alloca matrice con la fattorizzazione dei numeri b_smmoth
int**create_linear_system(int cardinality_f_base,mpz_t** matrix_B_smooth,int num_B_smooth);
long*alloc_array_long(int lenght);
int*alloc_array_int(int lenght);
int**alloc_array_pointer_to_int(int lenght);
int**alloc_matrix_int(int num_row,int num_col);
void divide_elem_multiple_of_2_by_x(mpz_t*vector,int length,double x);//vector = num1 exp1 num2 exp2 ecc
int**alloc_matrix_int(int num_row,int num_col);
long**alloc_matrix_long(int num_row,int num_col);
mpz_t**alloc_matrix_mpz(int num_row,int num_col);
char*alloc_array_char(long lenght);
char find_pivot_not_null_in_column(int **matrix,int num_row,int num_col,int pivot[2]);
int*get_coli(int **matrix,int num_row,int num_col,int index_col);
int scan_array_to_find_element_not_null(int*array,int start,int lenght_array);
void swap_row(int**matrix,int num_row,int num_col,int row1,int row2);
void reduce_echelon_form(int**matrix,int num_row,int num_col);
int count_rows_not_null(int **matrix,int num_row,int num_col);
char* find_free_var(int**matrix,int num_row,int num_col);
int calculate_dim_sol(int**matrix,int num_row,int num_col);
int* alloc_vector_base(char *array_var,int lenght,int num_vector_base);
int**calculate_base_linear_system(int**matrix,int num_row,int num_col,int*dim_sol);
void calculate_vector_base(int **matrix,int num_row,int num_col,char*array_var,int*v);
mpz_t*alloc_array_mpz(int length);
void sort_residuos(long size,mpz_t**matrix_factorization,int cardinality_factor_base,long M);
void copy_array_mpz(mpz_t*array1,mpz_t*array2,int length);
void find_relation_one_partial(mpz_t**matrix_factorization,long size,int cardinality_factor_base,const mpz_t n,int k,char*factorized);
void concatenate_matrix_mpz(mpz_t***pmatrix1,int *num_row1,int num_col1,mpz_t**matrix2,int num_row2,int num_col2);
void concatenate_matrix_long(long***res,long**matrix1,int num_row1,int num_col1,long**matrix2,int num_row2,int num_col2);
mpz_t***alloc_array_matrix_mpz(int length_array_matrix);
void concatenate_all_matrix_mpz_same_dimension(mpz_t***array_matrix_mpz,int length_array_matrix,int *row_result,int num_row_matrix,int num_col_matrix);
void adjust_array_of_prime(long*array_of_prime_chosen_for_a,int length_array,long*best_q,int length_best_q);
int count_number_potential_B_smooth_matrix_sorted(mpz_t**single_matrix_factorization,int num_row,int num_col);
void copy_matrix(mpz_t**matrix1,mpz_t**matrix2,int num_row,int num_col);
void sort_square(mpz_t**matrix_factorization,long size,int cardinality_factor_base,long M);
int count_number_B_smooth_matrix_unsorted(mpz_t**single_matrix_factorization,int num_row,int num_col);
void copy_residuos_from_matrix(mpz_t**potential_matrix_B_smooth,int potential_num_B_smooth,mpz_t**array_matrix_mpz,int num_row_matrix,int num_col_matrix);
char check_if_matrix_is_reduce_mod_n(int**matrix,int num_row,int num_col,int n);
mpz_t** create_matrix_B_smooth(int num_B_smooth,mpz_t**array_matrix_mpz,int num_row_matrix,int num_col_matrix);
char check_if_array_var_is_valid(int **matrix_linear_system,int num_row,int num_col,char*array_var);
char verify_solution(int **matrix_linear_system,int num_row,int num_col,int*solution);
char check_solution_base_matrix(int**linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base);
char check_if_array_is_reduce_mod_n(int*array,int length,int n);
void reduce_array_mod_n(int*array,int length,int n);
struct matrix_factorization*create_matrix_factorization_f(int M,int card_f_base,const mpz_t a,
const mpz_t b,const mpz_t n);
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s);
struct matrix_factorization**alloc_array_matrix_factorization(int length_array_matrix);
void concatenate_all_matrix_factorization_same_dimension(struct matrix_factorization**array_matrix_factorization,int length_array_matrix,int *row_result,int num_row_matrix);
int count_number_B_smooth_matrix_unsorted_f(struct matrix_factorization *single_matrix_factorization,int num_row);
struct matrix_factorization*create_matrix_B_smooth_f(int num_B_smooth,struct matrix_factorization mat,int num_row_matrix,int num_col_matrix);
double*alloc_array_double(int length);
void free_memory_matrix_factorization(struct matrix_factorization *mat);
float*alloc_array_float(int length);
void concatenate_all_matrix_B_smooth(struct matrix_factorization**array_matrix_factorization,int length_array_matrix,int *row_result);
struct matrix_factorization* alloc_matrix_factorization(int num_row);
char is_in_array_int(int*array,long length,long p_i);
char*create_linear_system_f(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth);
int**calculate_base_linear_system_char(char*matrix_linear_system,int num_row,int num_col,int*dim_sol);
long get_index(int index_row,int index_col,int num_col);
void reduce_echelon_form_char(char*matrix,int num_row,int num_col);
int count_rows_not_null_char(char*matrix,int num_row,int num_col);
char check_solution_base_matrix_char(char*linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base);
char verify_solution_char(char*matrix_linear_system,int num_row,int num_col,const int*solution);
void divide_vector_multiple_of_2_by_2(int*vector,int length);
char**alloc_matrix_char(int num_row,int num_col);
void copy_matrix_with_array(char**linear_system2,char*linear_system,int num_row,int num_col);
void reduce_echelon_form_matrix_char(char**matrix,int num_row,int num_col);
char value_is_in_sorted_array(int index_of_prime,struct a_struct*array_a_struct,int length);
#endif




