#include "Bsmooth_functions.h"

void print_index_B_smooth(struct node*head){
	printf("list index_of number B_smooth:\n");
	print_list_mpz(head);
	return;
}
void print_matrix_B_smooth(mpz_t**matrix_B_smooth,int num_of_B_smooth,int cardinality_factor_base){//parametri matrice,righe,colonne
	if(matrix_B_smooth==NULL || *matrix_B_smooth==NULL){
		printf("matrix B_smooth is empty\n");
	}
	if(num_of_B_smooth <=0 || cardinality_factor_base <=0){
		handle_error_with_exit("error in print_matrix_B_smooth\n");
	}
	printf("matrix_B_smooth:\n");
	print_matrix_mpz(matrix_B_smooth,num_of_B_smooth,cardinality_factor_base*2+1);
	return;
}