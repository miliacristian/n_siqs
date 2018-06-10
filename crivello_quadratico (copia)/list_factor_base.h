#include <gmp.h>
#ifndef DIN_F_H
#define DIN_F_H
struct node_factor_base  {//struttura di un nodo della lista dinamica ordianta,la lista tiene traccia delle trasmissioni e delle ritrasmissioni
    int prime;
    struct node_factor_base* next;//puntatore al prossimo
    struct node_factor_base* prev;//puntatore al precedente
    mpz_t root_n_mod_prime;
};

#endif
int count_element_linked_list_f(struct node_factor_base*head);
void get_element_linked_list_f(int *elem,struct node_factor_base*head,int index);
void remove_after_node_f(struct node_factor_base**ppos,struct node_factor_base**tail);
int delete_head_f(struct node_factor_base** head);
void insert_first_f(struct node_factor_base *new_node, struct node_factor_base **head, struct node_factor_base **tail);
void insert_at_tail_f(struct node_factor_base *new_node,struct node_factor_base**head,struct node_factor_base** tail);


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factor_base* get_new_node_f(int num,const mpz_t n);


void insert_at_head_f(struct node_factor_base* new_node,struct node_factor_base** head,struct node_factor_base** tail);


char first_is_smaller_f(struct node_factor_base node1, struct node_factor_base node2);

void insert_ordered_f(int num,const mpz_t n, struct node_factor_base** head, struct node_factor_base** tail);
void free_memory_list_f(struct node_factor_base*head);
