#ifndef MATH_H
#define MATH_H

#include "dynamic_list.h"
#include <gmp.h>
#include "basic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "gmp.h"
#include "matrix_function.h"
#include "list_factorization.h"
#include "list_square_relation.h"
#include <unistd.h>
#include <mpfr.h>
#include "parameters.h"
void square_root_mod_p_to_k(mpz_t rootpk,const mpz_t x,long p,const mpz_t n,int k);
long power_mod_n(long base,long exponent,long n);// base^exponent mod n

void select_square_mod_p(struct node**p,struct node** tail,const mpz_t n);

char is_prime(const mpz_t num,struct node*head);

void create_num(mpz_t num,const mpz_t a,const mpz_t b,const mpz_t n,long j);
char divide_all_by_p_to_k(const mpz_t r,long p,int index_of_prime,long k,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,const mpz_t n,const mpz_t a,const mpz_t b);//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione
	
char divide_all_by_min1(mpz_t*array_of_number,long M,mpz_t**matrix_factorization);//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione	
void factor_matrix(const mpz_t n,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,struct node*head_f_base,int cardinality_factor_base,const mpz_t a,const mpz_t b);

int count_number_B_smooth(const mpz_t*array_number,long M,struct node**head_index_B_smooth,struct node**tail_index_B_smooth,const mpz_t a,const mpz_t b,struct node**head_square_B_smooth,struct node**tail_square_B_smooth);
long**save_element_B_smooth(const mpz_t*array_number,mpz_t**matrix_factorization,long M,int cardinality_f_base,int* num_of_B_smooth,const mpz_t a,const mpz_t b,struct node**head_square_B_smooth,struct node**tail_square_B_smooth);//alloca matrice con la fattorizzazione dei numeri b_smmoth
void factor_matrix_f(const mpz_t n,long M,struct thread_data thread_data,int cardinality_factor_base,const mpz_t a,
                     struct a_struct*array_a_struct,int s);
int quadratic_residue(mpz_t x,const mpz_t q,const mpz_t n);
void add_exponent_of_a(long**matrix_B_smooth,int num_B_smooth,int s,const mpz_t a,long* array_of_prime_chosen_for_a);
void adjust_n(mpz_t n,int *k);
void reduce_array_mpz_mod_n(mpz_t*array,int length,const mpz_t a);
void set_to_odd_array_mpz_mod_n(mpz_t *array,int length,const mpz_t a);
mpz_t**save_element_B_smooth_in_matrix(mpz_t**matrix_factorization,long size,int cardinality_f_base,int* num_of_B_smooth);
float calculate_log_thresold(const mpz_t n,long M);
void find_list_square_relation(struct thread_data thread_data, int *num_B_smooth,int*num_semi_B_smooth, int *num_potential_B_smooth, long M,
                               struct node_square_relation **head_square, struct node_square_relation **tail_square,struct node_square_relation **head_residuos, struct node_square_relation **tail_residuos,
                               const mpz_t n,const mpz_t a,struct a_struct*array_a_struct,int s);
void calculate_index_min_max_a(int*number_prime_a,int*index_prime_a,int length,int*min_a,int*max_a);
struct a_struct*create_array_a_struct(int*number_prime_a,int*index_number_a,int length);
void*thread_factorization_job(void*arg);
char divide_all_by_p_to_k_with_thread(int rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b,struct a_struct*array_a_struct);
struct node_factorization*factorize_num_v2(const mpz_t num,int j_of_num,int first_index_f_base,int last_index_f_base,
                                           char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,struct a_struct*array_a_struct,int s,struct thread_data thread_data);
struct node_factorization*factorize_num_v1(const mpz_t num,int first_index_f_base,int last_index_f_base,
                                           char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,struct a_struct*array_a_struct,int s);
char calculate_num_from_factorization(mpz_t num_temp,struct node_factorization*head_factor);
char verify_square_relation(struct square_relation square_relation,const mpz_t n);
void calculate_thresold_large_prime(mpz_t thresold_large_prime,int max_prime);
char divide_all_by_2_log(long M,struct thread_data thread_data);
#endif











