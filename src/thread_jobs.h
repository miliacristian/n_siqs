#include "square_relations_functions.h"
#include <pthread.h>
#ifndef CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
#define CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
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
#endif //CRIVELLO_QUADRATICO__COPIA_1_THREAD_JOBS_H
