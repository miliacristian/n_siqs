#include "list_square_relation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "list_factorization.h"
#include "print.h"
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
char first_is_smaller_sort_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.square,node2.square_relation.square)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}
char first_is_smaller_residuos_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.residuos,node2.square_relation.residuos)>=0){//node1 è più grande di node2
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
void sort_relation_by_residuos(struct node_square_relation*head,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos){
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL) {
        q = p->next;
        insert_ordered_residuos_square_rel(p->square_relation, head_sort_residuos, tail_sort_residuos);
        free(p);
        p = q;
    }
    return;
}
struct square_relation create_relation_large_prime(struct square_relation rel1,struct square_relation rel2,mpz_t n){
    struct square_relation new_relation;
    int number,exp_of_number,index;
    struct node_factorization*tail=NULL;

    mpz_init(new_relation.square);
    mpz_init(new_relation.residuos);
    mpz_init(new_relation.num);
    mpz_set_si(new_relation.residuos,1);//residuo=1
    mpz_set_si(new_relation.num,0);//num=0
    new_relation.head_factorization=NULL;
    //square=square1*square2*(residuos)^-1
    mpz_invert(new_relation.square,rel1.residuos,n);//square=residuos^-1 mod n
    mpz_mul(new_relation.square,new_relation.square,rel1.square);
    mpz_mod(new_relation.square,new_relation.square,n);//square=residuos^-1*square1 mod n
    mpz_mul(new_relation.square,new_relation.square,rel2.square);
    mpz_mod(new_relation.square,new_relation.square,n);//square=residuos^-1*square1*square2 mod n

    print_struct_square_relation(rel1);
    print_struct_square_relation(rel2);
    struct node_factorization*p1=rel1.head_factorization;
    struct node_factorization*p2=rel2.head_factorization;
    while(p1!=NULL || p2!=NULL){//finisci quando sei arrivato alla fine di entrambe le liste
        if(p1!=NULL && p2!=NULL && p1->index<p2->index){
            //se le liste non sono finite e indice minore nella lista 1 aggiungi il nodo della lista 1
            index=p1->index;
            exp_of_number=p1->exp_of_number;
            number=p1->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p1=p1->next;
        }
        else if(p1!=NULL && p2!=NULL && p1->index==p2->index){
            //se le liste non sono finite e indici uguali somma esponenti
            if(p1->number!=-1) {
                index = p1->index;
                exp_of_number = p1->exp_of_number + p2->exp_of_number;
                number = p1->number;
                insert_ordered_factor(number, exp_of_number, index, &(new_relation.head_factorization), &tail);
            }
            p1=p1->next;
            p2=p2->next;
        }
        else if(p1!=NULL && p2!=NULL && p1->index>p2->index){
            //se le liste non sono finite e indice della lista 2 è minore,aggiungi nodo della lista 2
            index=p2->index;
            exp_of_number=p2->exp_of_number;
            number=p2->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p2=p2->next;
        }
        else if(p1==NULL && p2!=NULL){//se la prima lista è finita aggiungi nodo della seconda lista
            index=p2->index;
            exp_of_number=p2->exp_of_number;
            number=p2->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p2=p2->next;
        }
        else if(p1!=NULL && p2==NULL){//se la seconda lista è finita aggiungi nodo della prima lista
            index=p1->index;
            exp_of_number=p1->exp_of_number;
            number=p1->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p1=p1->next;
        }
        else{
            handle_error_with_exit("error in create relation large prime,caso non gestito\n");
        }
    }
    printf("new relation\n");
    print_struct_square_relation(new_relation);
    return new_relation;
}
void insert_ordered_sort_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
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
    if(first_is_smaller_sort_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_sort_square_rel(*temp,*new_node)){
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
void combine_relation_B_smooth_and_semi_B_smooth(struct node_square_relation*head,struct node_square_relation*tail,struct node_square_relation**head_final_list_relation,struct node_square_relation**tail_final_list_relation,mpz_t n,int*num_B_smooth){
    if(head==NULL || tail==NULL || head_final_list_relation==NULL || tail_final_list_relation==NULL || num_B_smooth==NULL){
        handle_error_with_exit("error in combine_relation_B_smooth and semi_B_smooth\n");
    }
    struct node_square_relation*head_sort_residuos=NULL;
    struct node_square_relation*tail_sort_residuos=NULL;
    sort_relation_by_residuos(head,&head_sort_residuos,&tail_sort_residuos);
    char not_remove_residuos_one=1;
    struct node_square_relation*p=head_sort_residuos;//p=nodo
    struct node_square_relation*q;
    while(p!=NULL) {
        while(p!=NULL && not_remove_residuos_one && mpz_cmp_si(p->square_relation.residuos,1)==0){//se il residuo è uguale a 1
            gmp_printf("residuo uguale a 1 square=%Zd\n",p->square_relation.square);
            insert_ordered_sort_square_rel(p->square_relation,head_final_list_relation,tail_final_list_relation);
            q=p->next;
            free(p);
            p=q;
        }
        not_remove_residuos_one=0;
        q=p->next;
        while (p!=NULL && q!=NULL && mpz_cmp(p->square_relation.residuos, q->square_relation.residuos) == 0) {//residui uguali
            gmp_printf("residui uguali a %Zd\n", p->square_relation.residuos);
            (*num_B_smooth)++;
            struct square_relation new_square_relation=create_relation_large_prime(p->square_relation,q->square_relation,n);
            insert_ordered_sort_square_rel(new_square_relation,head_final_list_relation,tail_final_list_relation);
            q = q->next;
        }
        q=p->next;
        mpz_clear(p->square_relation.square);
        mpz_clear(p->square_relation.num);
        mpz_clear(p->square_relation.residuos);
        free_list_factorization(p->square_relation.head_factorization);
        free(p);
        p=q;
    }
    return;
}

void remove_same_num(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth){
    if(head==NULL || tail==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL){
        handle_error_with_exit("error in remove_same_num\n");
    }
    struct node_square_relation*l=*head;
    while(l!=NULL){
        while(l->next!=NULL && (mpz_cmp(l->square_relation.num,(l->next)->square_relation.num)==0) && (mpz_cmp_si(l->square_relation.num,0)!=0 )) {
            if (mpz_cmp_si(l->square_relation.residuos, 1) == 0){//trovati due numeri uguali con residuo uguale a 1
                (*num_B_smooth)--;
            }
            else{//trovati due numeri uguali con residuo diverso da 1
                (*num_semi_B_smooth)--;
            }
            remove_after_node_square_rel(&(l->next),tail);
        }
        if(l->next==NULL){
            return;
        }
        l=l->next;
    }
    return;
}
void remove_same_square(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth){
    if(head==NULL || tail==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL){
        handle_error_with_exit("error in remove_same_num\n");
    }
    struct node_square_relation*l=*head;
    while(l!=NULL){
        while(l->next!=NULL && (mpz_cmp(l->square_relation.square,(l->next)->square_relation.square)==0) ) {
            printf("remove\n");
            if (mpz_cmp_si(l->square_relation.residuos, 1) == 0){//trovati due numeri uguali con residuo uguale a 1
                (*num_B_smooth)--;
            }
            else{//trovati due numeri uguali con residuo diverso da 1
                (*num_semi_B_smooth)--;
            }
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
        struct node_square_relation*q=p->next;
        insert_ordered_square_rel(p->square_relation,head,tail);
        free(p);
        p=q;
    }
    return;
}
void insert_ordered_residuos_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
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
    if(first_is_smaller_residuos_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_residuos_square_rel(*temp,*new_node)){
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