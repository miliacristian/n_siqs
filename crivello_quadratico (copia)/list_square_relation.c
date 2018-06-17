#include "list_square_relation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "list_factorization.h"
#include <gmp.h>
#include <unistd.h>

int count_element_linked_list_square_rel(struct node_square_relation*head){
    int count=0;
    if (head==NULL){
        handle_error_with_exit("head is NULL\n");
    }
    struct node_square_relation *p=head;
    while(p!=NULL){
        count++;
        p=p->next;
    }
    return count;
}
/*void get_element_linked_list_square_rel(mpz_t elem,struct node_square_relation*head,int index){//index start at 0
    if(index<0 || head==NULL || elem==NULL){
        handle_error_with_exit("index must be zero or positive\n");
    }
    struct node_square_relation*p=head;
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
void remove_after_node_square_rel(struct node_square_relation**ppos,struct node_square_relation**tail){
    if(ppos==NULL || tail==NULL ){
        handle_error_with_exit("error in parameter remove_after_node\n");
    }
    if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
    }
    struct node_square_relation *r = *ppos;
    struct node_square_relation*q=r->next;
    if(q==NULL){//fine della lista,bisogna aggiornare la coda
        *ppos=r->next;
        *tail=r->prev;
        mpz_clear(r->square_relation.num);
        mpz_clear(r->square_relation.square);
        free_memory_list_factor((r)->square_relation.head_factorization);
        free(r);
        r=NULL;
    }
    else {
        q->prev = r->prev;
        *ppos = r->next;
        mpz_clear(r->square_relation.num);
        mpz_clear(r->square_relation.square);
        free_memory_list_factor((r)->square_relation.head_factorization);
        free(r);
        r=NULL;
    }
    return;
}

int delete_head_square_rel(struct node_square_relation** head){//non è importante il valore iniziale di oldhead
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
        struct node_square_relation*temp=(*head)->next;
        free_memory_list_factor((*head)->square_relation.head_factorization);
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first_square_rel(struct node_square_relation *new_node, struct node_square_relation **head, struct node_square_relation **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail_square_rel(struct node_square_relation *new_node,struct node_square_relation**head,struct node_square_relation** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_square_relation* get_new_node_square_rel(struct square_relation square_relation) {
    if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in parameter get_new_node\n");
    }
    struct node_square_relation* new_node = (struct node_square_relation*)malloc(sizeof(struct node_square_relation));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    new_node->square_relation=square_relation;
    new_node->prev = NULL;
    new_node->next = NULL;
    return new_node;
}


void insert_at_head_square_rel(struct node_square_relation* new_node,struct node_square_relation** head,struct node_square_relation** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.num,node2.square_relation.num)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}
/*void union_list_square(struct node_square_relation**head,struct node_square_relation**tail,struct node_square_relation*head_square,struct node_square_relation*tail_square){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    if(head==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    if(head_square==NULL){
        return;
    }
    if(*tail==NULL){//lista vuota,prendi l'altra lista
        *head=head_square;
        *tail=tail_square;
        return;
    }
    head_square->prev=(*tail);//il primo nodo della seconda lista punta al nodo ultimo della lista
    (*tail)->next=head_square;//l'ultimo nodo punta al primo nodo dell'altra lista
    *tail=tail_square;//la coda punta alla coda dell'altra lista
    return;
}*/
void remove_same_num(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth){
    if(head==NULL || tail==NULL || num_B_smooth==NULL){
        handle_error_with_exit("error in remove_same_num\n");
    }
    struct node_square_relation*l=*head;
    while(l!=NULL){
        while(l->next!=NULL && (mpz_cmp(l->square_relation.num,(l->next)->square_relation.num)==0) ){
            (*num_B_smooth)--;
            remove_after_node_square_rel(&(l->next),tail);
        }
        if(l->next==NULL){
            return;
        }
        l=l->next;
    }
    return;
}
void union_list_square(struct node_square_relation**head,struct node_square_relation**tail,struct node_square_relation*head_square,struct node_square_relation*tail_square){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    if(head==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    if(head_square==NULL){
        return;
    }
    if(*tail==NULL){//lista vuota,prendi l'altra lista
        *head=head_square;
        *tail=tail_square;
        return;
    }
    struct node_square_relation*p=head_square;
    while(p!=NULL){
        insert_ordered_square_rel(p->square_relation,head,tail);
        p=p->next;
    }
    return;
}
void insert_ordered_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
     if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in insert_ordered\n");
     }
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    struct node_square_relation* temp = *tail;
    struct node_square_relation* next_node = NULL;
    struct node_square_relation* new_node = get_new_node_square_rel(square_relation);

    if(*head == NULL){
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    if(first_is_smaller_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_square_rel(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_square_rel(new_node, head, tail);
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
void free_memory_list_square_rel(struct node_square_relation*head){
    if(head==NULL){
        return;
    }
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL){
        q=p->next;
        free_memory_list_factor(p->square_relation.head_factorization);
        free(p);
        p=q;
    }
    return;
}