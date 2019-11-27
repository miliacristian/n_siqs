#include "square_relations_functions.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "array_numbers.h"
#include "basicv2.h"
#include <string.h>
#include "matrix_function_v2.h"
#include "properties.h"
#include "factor_base_functions.h"
#include "a_b_c_BK_functions.h"
#ifndef CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
#define CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H

extern int s;//numero di primi della factor base distinti che compongono a
extern struct factor_base_data*thread_factor_base_data;
extern struct a_struct*array_a_struct;
struct thread_data {
    float log_thresold;//valore della soglia
    mpz_t b;//coefficiente b
    int num_potential_B_smooth;
    int num_B_smooth;
    int num_semi_B_smooth;
    struct number*numbers;//array di struct di numberi
    struct node_square_relation*head_square;//ogni thread dopo il suo compito ha una lista di relazioni quadratiche
    struct node_square_relation*tail_square;//ogni thread dopo il suo compito ha una lista di relazioni quadratiche
    struct node_square_relation*head_residuos;
    struct node_square_relation*tail_residuos;
    long*j1_mod_p;//j1=-b+r)*a^-1 mod p
    long*j2_mod_p;//j2=(-b+r2)*a^-1 mod p

    //ogni volta che un thread analizza un dato appende la lista delle relazioni quadratiche a quella
    // precedentemente calcolata
};

struct factorization_thread_data{
    int id_thread;
    int start;
    int end;
    char is_a_default;
    pthread_mutex_t *mtx;
    void*pointer;
    struct thread_data thread_data;
};
void clear_struct_thread_data(struct thread_data t_data,int M);
void free_array_thread_data(struct thread_data*thread_data,int length_array_thread_data);
void join_all_threads(pthread_t*array_tid,int length_array);
int calculate_start_factor_base(int id_thread);
struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M);
int* create_factor_base_threads(pthread_t*array_tid,int num_thread);
#endif //CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
