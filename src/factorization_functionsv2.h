#ifndef X_H
#define X_H

#include "factor_base_functions.h"
#include "a_b_c_BK_functions.h"
struct node_factorization {
    int number;
    int exp_of_number;
    int index;
    struct node_factorization*next;
    struct node_factorization*prev;
};
struct row_factorization{
    int length;
	int*prime;//contiene primo factor base e
	int*log_prime;//valore del suo logaritmo
    int*root_n_mod_p;//prima radice di n mod p
    int*root2_n_mod_p;//seconda radice di n mod p
    int*inverse_a_mod_p;//a^-1 mod p
};
void free_list_factorization(struct node_factorization*head_factorization);
void print_factorization(const mpz_t num,struct node_factorization*head_factor);
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s);
char verify_factorization(const mpz_t num,mpz_t residuos,struct node_factorization*head_factor,const mpz_t a);
char calculate_num_from_factorization(mpz_t num_temp,struct node_factorization*head_factor);


#endif
