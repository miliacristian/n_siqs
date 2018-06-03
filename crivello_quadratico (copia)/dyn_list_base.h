#include <gmp.h>
#ifndef DIN_F_H
#define DIN_F_H
struct node_f  {//struttura di un nodo della lista dinamica ordianta,la lista tiene traccia delle trasmissioni e delle ritrasmissioni
    int prime;
    struct node_f* next;//puntatore al prossimo
    struct node_f* prev;//puntatore al precedente
};
#endif
int count_element_linked_list_f(struct node_f*head);
void get_element_linked_list_f(int *elem,struct node_f*head,int index);
void remove_after_node_f(struct node_f**ppos,struct node_f**tail);
int delete_head_f(struct node_f** head);
void insert_first_f(struct node_f *new_node, struct node_f **head, struct node_f **tail);
void insert_at_tail_f(struct node_f *new_node,struct node_f**head,struct node_f** tail);


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_f* get_new_node_f(int num);


void insert_at_head_f(struct node_f* new_node,struct node_f** head,struct node_f** tail);


char first_is_smaller_f(struct node_f node1, struct node_f node2);

void insert_ordered_f(int num, struct node_f** head, struct node_f** tail);
void free_memory_list_f(struct node_f*head);
