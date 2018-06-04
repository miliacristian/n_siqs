#include <gmp.h>
#ifndef CRIVELLO_QUADRATICO_COPIA_1_LIST_FACTORIZATION_H
#define CRIVELLO_QUADRATICO_COPIA_1_LIST_FACTORIZATION_H
struct node_factorization {
    int number;
    int exp_of_number;
    int index;
    struct node_factorization*next;
    struct node_factorization*prev;
};
#endif //CRIVELLO_QUADRATICO_COPIA_1_LIST_FACTORIZATION_H
int count_element_linked_list_factor(struct node_factorization*head);
void remove_after_node_factor(struct node_factorization**ppos,struct node_factorization**tail);

int delete_head_factor(struct node_factorization** head);

void insert_first_factor(struct node_factorization *new_node, struct node_factorization **head, struct node_factorization **tail);
void insert_at_tail_factor(struct node_factorization *new_node,struct node_factorization**head,struct node_factorization** tail);

//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factorization* get_new_node_factor(int number,int exp_of_number,int index);


void insert_at_head_factor(struct node_factorization* new_node,struct node_factorization** head,struct node_factorization** tail);


char first_is_smaller_factor(struct node_factorization node1, struct node_factorization node2);

void insert_ordered_factor(int number,int exp_of_number,int index, struct node_factorization** head, struct node_factorization** tail);
void free_memory_list_factor(struct node_factorization*head);
