#include <gmp.h>
#include "properties.h"
#include "factorization_functionsv2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "square_relations_functions.h"
#include <gmp.h>
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
char verify_cardinality_list_factorization(struct node_factorization*head,int length);

