#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"
#include "basic.h"

char not_print_list(int length,long thresold);
char not_print_array(long length,long thresold);
char not_print_matrix(long num_row,long num_col,long thresold_row,long thresold_col);
void print_array_char(char*array,int lenght,long thresold);
void print_array_long(long*array,int lenght,long thresold);
void print_array_int(int*array,int lenght,long thresold);
void print_vector_char(char*array,int lenght,long thresold);
void print_vector_long(long*array,int lenght,long thresold);
void print_vector_int(int*array,int lenght,long thresold);
void print_matrix_int(int**matrix,int num_row,int num_col,long thresold_row,long thresold_col);
void print_matrix_char(char**matrix,int num_row,int num_col,long thresold_row,long thresold_col);
void print_array_double(double*array,int lenght,long thresold);
void print_matrix_double(double**matrix,int num_row,int num_col,long thresold_row,long thresold_col);
void print_matrix_long(long**matrix,int num_row,int num_col,long thresold_row,long thresold_col);
void print_int(int a,char*string);
void print_char(char a,char*string);
void print_double(double a,char*string);
void print_long(long a,char*string);
void print_array_float(float*array,int length,long thresold);
void print_binary_matrix(unsigned long**binary_matrix,int num_row,int num_col,long thresold_row,long thresold_col);
void print_linear_system(char*linear_system,int num_rows,int num_cols,long thresold_col);
void print_base_linear_system(int**linear_system,int num_rows,int num_cols,long thresold_row,long thresold_col);

#endif
