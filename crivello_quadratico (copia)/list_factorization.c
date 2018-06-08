#include "list_factorization.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include <gmp.h>

int count_element_linked_list_factor(struct node_factorization*head){
    int count=0;
    if (head==NULL){
        handle_error_with_exit("head is NULL\n");
    }
    struct node_factorization *p=head;
    while(p!=NULL){
        count++;
        p=p->next;
    }
    return count;
}
/* da scrivere
  void get_element_linked_list_factor(mpz_t elem,struct node*head,int index){//index start at 0
    if(index<0 || head==NULL || elem==NULL){
        handle_error_with_exit("index must be zero or positive\n");
    }
    struct node*p=head;
    int count=0;
    while(p!=NULL){
        if(count==index){
            mpz_set(elem,p->prime);
            return;
        }
        count++;
        p=p->next;
    }
    return;
}*/
void remove_after_node_factor(struct node_factorization**ppos,struct node_factorization**tail){
    if(ppos==NULL || tail==NULL ){
        handle_error_with_exit("error in parameter remove_after_node\n");
    }
    if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
    }
    struct node_factorization *r = *ppos;
    struct node_factorization*q=r->next;
    if(q==NULL){//fine della lista,bisogna aggiornare la coda
        *ppos=r->next;
        *tail=r->prev;
        free(r);
        r=NULL;
    }
    else {
        q->prev = r->prev;
        *ppos = r->next;
        free(r);
        r=NULL;
    }
    return;
}

int delete_head_factor(struct node_factorization** head){//non è importante il valore iniziale di oldhead
    //initializza oldhead con il primo nodo della lista e distrugge il primo nodo della lista
    if(head==NULL){
        handle_error_with_exit("error in delete head\n");
    }
    if(*head == NULL){//nessuna testa da eliminare
        return -1;
    }
    if ((*head)-> next == NULL){//c'è un solo nodo in lista
        free(*head);
        *head = NULL;
    }else{
        struct node_factorization*temp=(*head)->next;
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first_factor(struct node_factorization *new_node, struct node_factorization **head, struct node_factorization **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail_factor(struct node_factorization *new_node,struct node_factorization**head,struct node_factorization** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_factor(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factorization* get_new_node_factor(int number,int exp_of_number,int index) {
    if( (number<=0 && number!=-1) || exp_of_number<0 || index<0){
        handle_error_with_exit("error in parameter get_new_node\n");
    }
    struct node_factorization* new_node = (struct node_factorization*)malloc(sizeof(struct node_factorization));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    new_node->exp_of_number=exp_of_number;
    new_node->index=index;
    new_node->number=number;
    new_node->prev = NULL;
    new_node->next = NULL;
    return new_node;
}


void insert_at_head_factor(struct node_factorization* new_node,struct node_factorization** head,struct node_factorization** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_factor(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller_factor(struct node_factorization node1, struct node_factorization node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(node1.index>node2.index){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}

void insert_ordered_factor(int number,int exp_of_number,int index, struct node_factorization** head, struct node_factorization** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    if( (number<=0 && number!=-1) || exp_of_number<0 || index<0){
        handle_error_with_exit("error in parameter get_new_node\n");
    }
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    struct node_factorization* temp = *tail;
    struct node_factorization* next_node = NULL;
    struct node_factorization* new_node = get_new_node_factor(number,exp_of_number,index);
    if(*head == NULL){
        insert_first_factor(new_node, head, tail);
        return;
    }

    if(first_is_smaller_factor((**tail),*new_node)){
        insert_at_tail_factor(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_factor(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_factor(new_node, head, tail);
                return;
            }
        }
        next_node = temp->next;
        new_node->prev = temp;
        new_node->next = next_node;
        temp->next = new_node;
        next_node->prev = new_node;
    }
    return;
}
void free_memory_list_factor(struct node_factorization*head){
    if(head==NULL){
        return;
    }
    struct node_factorization*p=head;
    struct node_factorization*q;
    while(p!=NULL){
        q=p->next;
        free(p);
        p=q;
    }
    return;
}

