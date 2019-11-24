#include "properties.h"
#include <gmp.h>
#ifndef DIN_H
#define DIN_H
struct node  {//struttura di un nodo della lista dinamica ordianta,la lista tiene traccia delle trasmissioni e delle ritrasmissioni
    mpz_t prime;
    struct node* next;//puntatore al prossimo
    struct node* prev;//puntatore al precedente
};
#endif
void insert_first(struct node *new_node, struct node **head, struct node **tail);//inserisce il primo nodo

void insert_at_tail(struct node *new_node,struct node**head,struct node** tail);//inserisce un nodo in coda

//alloca e inizializza un nodo della lista dinamica ordinata
struct node* get_new_node(mpz_t num);

void remove_after_node(struct node**ppos,struct node**tail);
void insert_at_head(struct node* new_node,struct node** head,struct node** tail);//inserisce un nodo in testa alla lista


char first_is_smaller(struct node node1, struct node node2);//verifica se il primo nodo contiene tempi più piccoli del secondo nodo

void insert_ordered(mpz_t num, struct node** head, struct node** tail);
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali

int count_element_linked_list(struct node*head);
void get_element_linked_list(mpz_t elem,struct node*head,int index);//index start at 0
int delete_head(struct node** head);//non è importante il valore iniziale di oldhead
    //initializza oldhead con il primo nodo della lista e distrugge il primo nodo della lista
void free_memory_list(struct node*head);
