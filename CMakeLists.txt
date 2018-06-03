cmake_minimum_required(VERSION 3.10)
project(crivello_quadratico__copia_1)

set(CMAKE_CXX_STANDARD 11)

include_directories("crivello_quadratico (copia)")

add_executable(crivello_quadratico__copia_1
        "crivello_quadratico (copia)/basic.c"
        "crivello_quadratico (copia)/basic.h"
        "crivello_quadratico (copia)/biprime.txt"
        "crivello_quadratico (copia)/criv_quad.c"
        "crivello_quadratico (copia)/criv_quad.h"
        "crivello_quadratico (copia)/domini.txt"
        "crivello_quadratico (copia)/dyn_list_base.c"
        "crivello_quadratico (copia)/dyn_list_base.h"
        "crivello_quadratico (copia)/dynamic_list.c"
        "crivello_quadratico (copia)/dynamic_list.h"
        "crivello_quadratico (copia)/factorization_biprime.txt"
        "crivello_quadratico (copia)/file_log.txt"
        "crivello_quadratico (copia)/main.c"
        "crivello_quadratico (copia)/main.h"
        "crivello_quadratico (copia)/makefile"
        "crivello_quadratico (copia)/math_function.c"
        "crivello_quadratico (copia)/math_function.h"
        "crivello_quadratico (copia)/matrix_function.c"
        "crivello_quadratico (copia)/matrix_function.h"
        "crivello_quadratico (copia)/miller_rabin.c"
        "crivello_quadratico (copia)/miller_rabin.h"
        "crivello_quadratico (copia)/my_pari_functions.txt"
        "crivello_quadratico (copia)/number.txt"
        "crivello_quadratico (copia)/pollard_rho.c"
        "crivello_quadratico (copia)/print.c"
        "crivello_quadratico (copia)/print.h"
        "crivello_quadratico (copia)/README.txt"
        "crivello_quadratico (copia)/run"
        "crivello_quadratico (copia)/run.c"
        "crivello_quadratico (copia)/runtest"
        "crivello_quadratico (copia)/runtest.c"
        "crivello_quadratico (copia)/valgrind-out.txt")
