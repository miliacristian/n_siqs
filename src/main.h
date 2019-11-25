#ifndef MAIN_H
#define MAIN_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "basic.h"
#include "math_function.h"
#include "dynamic_list.h"
#include "matrix_function.h"
#include "print.h"
#include "miller_rabin.h"
#include "criv_quad.h"
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <mpfr.h>
#include "list_factorization.h"
#include "list_square_relation.h"
#include <pthread.h>
#include "parameters.h"
int thread_job_criv_quad(int i);
int thread_job_to_create_factor_base(int id_thread);
#endif
