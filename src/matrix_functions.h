#ifndef MATRIX_H
#define MATRIX_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "print.h"
#include "math_functions.h"
#include "basic.h"
#include "timing.h"

long*alloc_array_long(int length);
char array_is_fill_of_value(char*combination,int length,char value);
void free_memory_matrix_unsigned_long(unsigned long **matrix,int num_row,int num_col);
unsigned long**alloc_matrix_unsigned_long(int num_row,int num_col);
void free_memory_matrix_int(int **matrix,int num_row,int num_col);
int* sum_vector(int*vector1,int*vector2,int lenght1,int lenght2);
long* sum_elem_multiple_of_2(long*vector1,long*vector2,int length1,int length2);
void prod_vector_and_scalar(int*vector,int scalar,int lenght);
int* prod_vector_and_scalar_v2(int*vector,int scalar,int lenght);
void reduce_matrix_mod_n(int**linear_system,int cardinality_factor_base,int num_of_B_smooth,int n);
int** alloc_matrix_int_and_materialize(int num_row,int num_col);
char is_in_array_long(long*array,long length,long p_i);
int** alloc_linear_system(int cardinality_f_base,int num_B_smooth);//alloca matrice con la fattorizzazione dei numeri b_smmoth
long* alloc_array_long(int lenght);
int* alloc_array_int(int lenght);
int** alloc_array_pointer_to_int(int lenght);
int** alloc_matrix_int(int num_row,int num_col);
int** alloc_matrix_int(int num_row,int num_col);
long** alloc_matrix_long(int num_row,int num_col);
char* alloc_array_char(long lenght);
int* get_coli(int **matrix,int num_row,int num_col,int index_col);
int scan_array_to_find_element_not_null(int*array,int start,int lenght_array);
void swap_row(int**matrix,int num_row,int num_col,int row1,int row2);
void reduce_echelon_form(int**matrix,int num_row,int num_col);
int count_rows_not_null(int **matrix,int num_row,int num_col);
char* find_free_var(int**matrix,int num_row,int num_col);
int calculate_dim_sol(int**matrix,int num_row,int num_col);
int* alloc_vector_base(char *array_var,int lenght,int num_vector_base);
int** calculate_base_linear_system(int**matrix,int num_row,int num_col,int*dim_sol);
void calculate_vector_base(int **matrix,int num_row,int num_col,char*array_var,int*v,long thresold);
void concatenate_matrix_long(long***res,long**matrix1,int num_row1,int num_col1,long**matrix2,int num_row2,int num_col2);
char check_if_matrix_is_reduce_mod_n(int**matrix,int num_row,int num_col,int n);
char check_if_array_var_is_valid(int **matrix_linear_system,int num_row,int num_col,char*array_var);
char verify_solution(int **matrix_linear_system,int num_row,int num_col,int*solution);
char check_solution_base_matrix(int**linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base);
char check_if_array_is_reduce_mod_n(int*array,int length,int n);
void reduce_array_mod_n(int*array,int length,int n);
double* alloc_array_double(int length);
float* alloc_array_float(int length);
char is_in_array_int(int*array,long length,long p_i);
int** calculate_base_linear_system_char(char*matrix_linear_system,int num_row,int num_col,int*dim_sol,int max_dim_sol);
long get_index(int index_row,int index_col,int num_col);
void reduce_echelon_form_char(char*matrix,int num_row,int num_col);
int count_rows_not_null_char(char*matrix,int num_row,int num_col);
char check_solution_base_matrix_char(char*linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base);
char verify_solution_char(char*matrix_linear_system,int num_row,int num_col,const int*solution);
void divide_vector_multiple_of_2_by_2(int*vector,int length);
char** alloc_matrix_char(int num_row,int num_col);
void copy_matrix_with_array(char**linear_system2,char*linear_system,int num_row,int num_col);
void reduce_echelon_form_matrix_char(char**matrix,int num_row,int num_col);
void swap_row_unsigned_long(unsigned long**matrix,int num_row,int ind_row1,int ind_row2);
void reduce_echelon_form_binary_matrix(unsigned long**binary_matrix,int num_row,int num_col);
char* from_matrix_binary_to_matrix_char(unsigned long**binary_linear_system,int num_row,int num_col_binary_matrix,int*num_col_linear_system);
int calculate_dim_solution(int num_row_not_null,int num_col_not_null);
pthread_t *alloc_array_tid(int num_thread);
#endif




