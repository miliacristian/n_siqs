#include <stdio.h>
#include <stdlib.h>
#include "factorization_functionsv2.h"

void free_list_factorization(struct node_factorization*head_factorization){
    if(head_factorization==NULL){
        return;
    }
    struct node_factorization*p=head_factorization;
    struct node_factorization*q;
    while(p!=NULL){
        q=p->next;
        free(p);
        p=q;
    }
    return;
}