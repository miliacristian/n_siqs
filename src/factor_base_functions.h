
#ifndef FAC_BASE_FUNC_H
#define FAC_BASE_FUNC_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpz_functions.h"
#include <gmp.h>
#include "parameters.h"
#include "math_functions.h"

struct node_factor_base  {//struttura di un nodo della lista dinamica ordianta,la lista tiene traccia delle trasmissioni e delle ritrasmissioni
    int prime;
    struct node_factor_base* next;//puntatore al prossimo
    struct node_factor_base* prev;//puntatore al precedente
    int root_n_mod_prime;
};
struct factor_base_data {
    struct node_factor_base*head;
    struct node_factor_base*tail;
    int cardinality_factor_base;
    int last_prime_factor_base;
};
void print_list_factor_base(struct node_factor_base*head,int length);
void print_factor_base(struct node*head_f_base);

int count_element_linked_list_f(struct node_factor_base*head);
void get_element_linked_list_f(int *elem,struct node_factor_base*head,int index);
void remove_after_node_f(struct node_factor_base**ppos,struct node_factor_base**tail);
int delete_head_f(struct node_factor_base** head);
void insert_first_f(struct node_factor_base *new_node, struct node_factor_base **head, struct node_factor_base **tail);
void insert_at_tail_f(struct node_factor_base *new_node,struct node_factor_base**head,struct node_factor_base** tail);


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factor_base* get_new_node_f(int num,const mpz_t n);


void insert_at_head_f(struct node_factor_base* new_node,struct node_factor_base** head,struct node_factor_base** tail);

char verify_factor_base(struct node_factor_base*head,int cardinality_factor_base,int last_prime_factor_base);
char first_is_smaller_f(struct node_factor_base node1, struct node_factor_base node2);

void insert_ordered_f(int num,const mpz_t n, struct node_factor_base** head, struct node_factor_base** tail);
void free_memory_list_f(struct node_factor_base*head);
void union_list_factor_base(struct node_factor_base**head1,struct node_factor_base**tail1,int*cardinality_factor_base1,int*last_prime_factor_base1,
                            struct node_factor_base*head2,struct node_factor_base*tail2,int cardinality_factor_base2,int last_prime_factor_base2);
char verify_cardinality_list_factor_base(struct node_factor_base*head,int length);
struct factor_base_data* alloc_array_factor_base_data(int length);
void create_factor_base_f(int*cardinality_factor_base,long B,struct node_factor_base**head,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base);
#endif
