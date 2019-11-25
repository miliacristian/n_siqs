#ifndef X_H
#define X_H

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
#endif
