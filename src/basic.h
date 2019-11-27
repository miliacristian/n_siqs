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

#include "thread_jobs.h"
#include "main.h"

int* create_threads(pthread_t*array_tid,int num_thread);
struct factorization_thread_data* create_factorization_threads(pthread_t*array_tid,struct thread_data thread_data,const mpz_t a,int num_thread);
#endif
