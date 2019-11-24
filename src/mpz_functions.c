#include "mpz_functions.h"
//print mpz
void print_array_mpz(mpz_t*array,int length){
	if(length<=0 || array==NULL){
		handle_error_with_exit("error in print_array_mpz\n");
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		gmp_printf("%Zd ",array[i]);
	}
	printf("\n");
	return;
}

void print_list_mpz(struct node*head){
	if (head==NULL){
		printf("impossible print list head is NULL\n");
		return;
	}
	struct node *p=head;
	while(p!=NULL){
		gmp_printf("%Zd ",p->prime);
		p=p->next;
	}
	printf("\n");
	return;
}
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col){
	if(matrix==NULL){
		handle_error_with_exit("pointer null\n");
	}
	if(*matrix==NULL){
		handle_error_with_exit("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_mpz\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			print_array_mpz(matrix[i],num_col);
	}
	return;
}