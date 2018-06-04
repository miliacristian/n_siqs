#include "list_factor_base.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include <gmp.h>

int count_element_linked_list_f(struct node_factor_base*head){
	int count=0;
	if (head==NULL){
		handle_error_with_exit("head is NULL\n");
	}
	struct node_factor_base *p=head;
	while(p!=NULL){
		count++;
		p=p->next;
	}
	return count;
}
void get_element_linked_list_f(int *elem,struct node_factor_base*head,int index){//index start at 0
	if(index<0 || head==NULL || elem==NULL){
		handle_error_with_exit("index must be zero or positive\n");
	}
	struct node_factor_base*p=head;
	int count=0;
	while(p!=NULL){
		if(count==index){
			*elem=p->prime;
			return;
		}
		count++;
		p=p->next;
	}
	return;
}
void remove_after_node_f(struct node_factor_base**ppos,struct node_factor_base**tail){
	if(ppos==NULL || tail==NULL ){
		handle_error_with_exit("error in parameter remove_after_node\n");
	}
	 if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
  }
	struct node_factor_base *r = *ppos;
	struct node_factor_base*q=r->next;
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

int delete_head_f(struct node_factor_base** head){//non è importante il valore iniziale di oldhead
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
	struct node_factor_base*temp=(*head)->next;
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first_f(struct node_factor_base *new_node, struct node_factor_base **head, struct node_factor_base **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail_f(struct node_factor_base *new_node,struct node_factor_base**head,struct node_factor_base** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_f(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factor_base* get_new_node_f(int num) {
	if(num<=0 && num!=-1){
		handle_error_with_exit("error in get_new_node\n");
	}
    struct node_factor_base* new_node = (struct node_factor_base*)malloc(sizeof(struct node_factor_base));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    new_node->prime=num;
    new_node->prev = NULL;
    new_node->next = NULL;
    return new_node;
}


void insert_at_head_f(struct node_factor_base* new_node,struct node_factor_base** head,struct node_factor_base** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_f(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller_f(struct node_factor_base node1, struct node_factor_base node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
	if(node1.prime>=node2.prime){
		return 0;
	}
    return 1;//node1 è più piccolo di node 2
}

void insert_ordered_f(int num, struct node_factor_base** head, struct node_factor_base** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    struct node_factor_base* temp = *tail;
    struct node_factor_base* next_node = NULL;
    struct node_factor_base* new_node = get_new_node_f(num);
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    if(*head == NULL){
        insert_first_f(new_node, head, tail);
        return;
    }
    if(first_is_smaller_f((**tail),*new_node)){
        insert_at_tail_f(new_node,head, tail);
    }
		else{
        while(!first_is_smaller_f(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_f(new_node, head, tail);
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
void free_memory_list_f(struct node_factor_base*head){
	if(head==NULL){
		return;
	}
	struct node_factor_base*p=head;
	struct node_factor_base*q;
	while(p!=NULL){
		q=p->next;
		free(p);
		p=q;
	}
	return;
}
