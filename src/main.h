#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "dynamic_list.h"
#include "miller_rabin.h"
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <mpfr.h>
#include "list_factorization.h"
#include "list_square_relations.h"
#include <pthread.h>
#include "parameters.h"
#include "mpz_functions.h"
#include "matrix_functions.h"
#include "basic.h"
#include "thread_jobs.h"
#include "factor_base_functions.h"
#include "M_and_B_functions.h"
#include "a_b_c_BK_functions.h"
#include "square_relations_functions.h"
#include "n.h"
#include "timing.h"
#include "memory_limit.h"
#if LINEAR_PINNING!=1
#include "numa.h"
#endif
int thread_job_criv_quad(int i);
int thread_job_to_create_factor_base(int id_thread);
#endif
