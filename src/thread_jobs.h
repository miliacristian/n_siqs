#ifndef THREAD_JOBS_H
#define THREAD_JOBS_H

#include "square_relations_functions.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "array_numbers.h"
#include "basic.h"
#include <string.h>
#include <errno.h>
#include "matrix_functions.h"
#include "factor_base_functions.h"
#include "a_b_c_BK_functions.h"
#include "parameters.h"
#include "M_and_B_functions.h"
#include "list_factorization.h"
#include "list_square_relations.h"

extern int s;//numero di primi della factor base distinti che compongono a
extern struct factor_base_data*thread_factor_base_data;
extern struct a_struct*array_a_struct;
struct thread_data  {
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
} __attribute__ ((aligned(64)));

struct factorization_thread_data{
    int id_thread;
    int start;
    int end;
    char is_a_default;
    pthread_mutex_t *mtx;
    void*pointer;
    struct thread_data thread_data;
};
void clear_struct_thread_data(struct thread_data *t_data);
void free_array_thread_data(struct thread_data*thread_data,int length_array_thread_data);
void join_all_threads(pthread_t*array_tid,int length_array);
int calculate_start_factor_base(int id_thread);
struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M);
int* create_factor_base_threads(pthread_t*array_tid,int num_thread);
int thread_job_criv_quad(int id_thread);
struct node_factorization*factorize_num(const mpz_t num,int j_of_num,int first_index_f_base,int last_index_f_base,
                                           char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,const struct a_struct*array_a_struct,int s,struct thread_data *pthread_data);
int* create_threads(pthread_t*array_tid,int num_thread);
void factor_matrix(const mpz_t n,long M,struct thread_data *thread_data,const int cardinality_factor_base,const mpz_t a,
                     const struct a_struct*array_a_struct,int s);
void find_list_square_relation(struct thread_data *thread_data, int *num_B_smooth,int*num_semi_B_smooth, int *num_potential_B_smooth, long M,
                               struct node_square_relation **head_square, struct node_square_relation **tail_square,struct node_square_relation **head_residuos, struct node_square_relation **tail_residuos,
                               const mpz_t n,const mpz_t a,const struct a_struct*array_a_struct,int s);
#endif //CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
