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

int* create_threads(pthread_t*array_tid,int num_thread);
struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M);
int* create_factor_base_threads(pthread_t*array_tid,int num_thread);
struct factorization_thread_data* create_factorization_threads(pthread_t*array_tid,struct thread_data thread_data,const mpz_t a,int num_thread);
#endif
