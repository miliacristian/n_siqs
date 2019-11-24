#include "factorization_functions.h"
void print_factorization(const mpz_t num,struct node_factorization*head_factor){
    if(head_factor==NULL){
        printf("no simple factorization found\n");
        return;
    }
	gmp_printf("%Zd=",num);
    struct node_factorization*p=head_factor;
    while(p!=NULL){
        printf("%d^%d,",p->number,p->exp_of_number);
        p=p->next;
    }
    printf("\n");
    return;
}
void print_list_square_relation(struct node_square_relation*head,int length){
	if (head==NULL){
        printf("impossible print list head is NULL\n");
        return;
    }
    if(length<0){
		handle_error_with_exit("error in print_list_square_relation\n");
	}
    if(length==0){
        printf("list is empty\n");
    }
    if(not_print_list(length)==1){
        return;
    }
    struct node_square_relation*p=head;
    while(p!=NULL){
        gmp_printf("square=%Zd,residuos=%Zd,",p->square_relation.square,p->square_relation.residuos);
        print_factorization(p->square_relation.num,p->square_relation.head_factorization);
        p=p->next;
    }
    printf("\n");
    return;
}
