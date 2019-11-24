#include "print.h"

//print siqs
void print_array_a_struct(struct a_struct*array_a_struct,int length){
	if(array_a_struct==NULL || length<=0){
		handle_error_with_exit("error in print_aray_a_struct\n");
	}
	for(int i=0;i<length;i++){
		printf("number=%d,index=%d\n",array_a_struct[i].number_prime_a,array_a_struct[i].index_prime_a);
	}
	return;
}
void print_array_Bk(mpz_t*array_Bk,long s){
	if(s<0){
		handle_error_with_exit("error in print_array_Bk\n");
	}
	if(array_Bk==NULL){
		printf("array_Bk is empty\n");
		return;
	}
	if(not_print_array(s)==1){
		return;
	}
	printf("array Bk:\n");
	print_array_mpz(array_Bk,s);
	return;
}
void print_array_bi(mpz_t* array_bi,long s){
	if(s<0){
		handle_error_with_exit("error in print_array_bi\n");
	}
	if(array_bi==NULL){
		printf("array_bi is empty\n");
		return;
	}
	long length=(long)pow(2,s-1);
	if(not_print_array(length)==1){
		return;
	}
	printf("array_bi:\n");
	print_array_mpz(array_bi,length);
	return;
	
}
void print_array_chosen_for_a(long*array_of_prime_chosen_for_a,int card_factor_base,int s){
	if(card_factor_base<=0){
		handle_error_with_exit("error in card factor base print_array_chosen_for_a\n");
	}
	if(array_of_prime_chosen_for_a==NULL){
		if(s==0){
			return;
		}
		else{
			handle_error_with_exit("error in print_array_chosen_for_a\n");
		}
	}
	printf("array of prime chosen for a:\n");
	print_array_long(array_of_prime_chosen_for_a,card_factor_base*2);
	return;
}
//print square relation
void print_struct_square_relation(struct square_relation square_relation){
	gmp_printf("square=%Zd,residuos=%Zd\n",square_relation.square,square_relation.residuos);
}






void print_matrix_factorization(mpz_t**matrix_factorization,int M,int cardinality_factor_base){
	if(matrix_factorization==NULL || *matrix_factorization==NULL){
		printf("matrix factorization is empty\n");
	}
	if(M <=0 || cardinality_factor_base <=0){
		handle_error_with_exit("error in print_matrix_factorization\n");
	}
	printf("matrix_factorization:\n");
	print_matrix_mpz(matrix_factorization,2*M+1,cardinality_factor_base*2+2);
	return;
}
void print_thread_data(struct thread_data thread_data,long M,int cardinality_factor_base){
	if(M<=0 || cardinality_factor_base<=0){
		handle_error_with_exit("error in print_thread_data\n");
	}
	if(M>THRESOLD_PRINT_ARRAY/2){
		return;
	}
	for(int i=0;i<2*M+1;i++){
		gmp_printf("b=%Zd,",thread_data.b);
		printf("first_index=%d,last_index=%d,sum_log=%d,j=%d\n",thread_data.numbers[i].first_index_f_base,thread_data.numbers[i].last_index_f_base,thread_data.numbers[i].sum_log,thread_data.numbers[i].j);
	}
	print_array_long(thread_data.j1_mod_p,cardinality_factor_base);
	print_array_long(thread_data.j2_mod_p,cardinality_factor_base);
}

void print_array_matrix_same_dimension(mpz_t***array_matrix_mpz,int length_array,int num_row,int num_col){
	if(array_matrix_mpz==NULL || *array_matrix_mpz==NULL || **array_matrix_mpz==NULL || length_array<=0
	   || num_row<=0 ||  num_col<=0){
		handle_error_with_exit("error in print_array_mpz_same_dimension\n");
	}
	for(int i=0;i<length_array;i++){
		print_matrix_mpz(array_matrix_mpz[i],num_row,num_col);
	}
	return;
}

void print_array_number(mpz_t*array,int M){
	if(M<=0){
		handle_error_with_exit("error in print_array_number\n");
	}
	if(array==NULL){
		printf("array_number is empty\n");
		return;
	}
	if(not_print_array(2*M+1)==1){
		return;
	}
	printf("array number:\n");
	print_array_mpz(array,2*M+1);
	return;
}






