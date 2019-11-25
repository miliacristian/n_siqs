#ifndef LINE_H
#define LINE_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>
#include <gmp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "properties.h"
#include "timing.h"
#include "parameters.h"
#include "square_relations_functions.h"
#include "list_factorization.h"
#include "array_numbers.h"
#include "a_b_c_BK_functions.h"
#include "thread_jobs.h"

char get_and_check_n(mpz_t n,FILE*file_number);
int* create_threads(pthread_t*array_tid,int num_thread);
void join_all_threads(pthread_t*array_tid,int length_array);
void test();
struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M);
void clear_struct_thread_data(struct thread_data t_data,int M);
void free_array_thread_data(struct thread_data*thread_data,int length_array_thread_data);
void free_memory_list_square_relation(struct node_square_relation*head);
struct factor_base_data*alloc_array_factor_base_data(int length);
int* create_factor_base_threads(pthread_t*array_tid,int num_thread);
struct factorization_thread_data* create_factorization_threads(pthread_t*array_tid,struct thread_data thread_data,const mpz_t a,int num_thread);
void free_list_factorization(struct node_factorization*head_factorization);
int calculate_start_factor_base(int id_thread);
int compare_a_struct( const void* a, const void* b);
void print_estimated_time(int cardinality_factor_base,int num_B_smooth);
#endif
