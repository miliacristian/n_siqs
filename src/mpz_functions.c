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
//funzioni per allocare array o matrici 
void copy_array_mpz(mpz_t*array1,mpz_t*array2,int length){
	if(array1==NULL || array2==NULL || length<=0){
		handle_error_with_exit("error in copy_array_mpz\n");
	}
	for(int i=0;i<length;i++){
		mpz_set(array1[i],array2[i]);
	}
	return;
}
void copy_matrix_mpz(mpz_t**matrix1,mpz_t**matrix2,int num_row,int num_col){
	if(matrix1==NULL || *matrix1==NULL || matrix2==NULL || *matrix2==NULL
	|| num_row<=0 || num_col<=0){
		handle_error_with_exit("error in copy_array_mpz\n");
	}
	for(int i=0;i<num_row;i++){
		copy_array_mpz(matrix1[i],matrix2[i],num_col);//copia riga
	}
	return;
}

mpz_t*alloc_array_mpz(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_mpz\n");
	}
	mpz_t *array=malloc(sizeof(mpz_t)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc\n");
	}
	for(int i=0;i<length;i++){
		mpz_init(array[i]);
	}
	return array;
}

mpz_t***alloc_array_matrix_mpz(int length_array_matrix){
	mpz_t***array_matrix=NULL;
	if(length_array_matrix<=0){
		handle_error_with_exit("error in parameter alloc aray_matrix_mpz\n");
	}
	array_matrix=malloc(sizeof(mpz_t**)*length_array_matrix);
	if(array_matrix==NULL){
		handle_error_with_exit("error in malloc alloc array_matrix_mpz\n");
	}
	memset(array_matrix,0,sizeof(mpz_t**)*length_array_matrix);
	return array_matrix;
}

mpz_t**realloc_array_pointer_to_mpz(mpz_t**t,int new_length){
	if(new_length<=0 || t==NULL || *t==NULL){
		handle_error_with_exit("error in parameteralloc array pointer to mpz\n");
	}
	mpz_t **array=realloc(t,sizeof(mpz_t*)*(new_length));
	if(array==NULL){
		handle_error_with_exit("error in realloc alloc array pointer to long\n");
	}
	return array;
}

mpz_t**alloc_array_pointer_to_mpz(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to mpz\n");
	}
	mpz_t **array=malloc(sizeof(mpz_t*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to long\n");
	}
	memset(array,0,sizeof(mpz_t*)*(length));
	return array;
}

mpz_t **alloc_matrix_mpz(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_long\n");
	}
	mpz_t**matrix;
	matrix=alloc_array_pointer_to_mpz(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_mpz(num_col);
	}
	return matrix;
}
void copy_column_matrix_mpz(mpz_t *v,mpz_t**matrix_factorization,long num_row,int index_col){
	if(v==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || num_row<=0 || index_col<0){
		handle_error_with_exit("error in copy_column matrix mpz\n");
	}
	for(long i=0;i<num_row;i++){
		mpz_set(v[i],matrix_factorization[i][index_col]);
	}
	return;
}

void sum_elem_multiple_of_2_mpz(mpz_t*vector1,mpz_t*vector2,int length1,int length2){//vettore risultato=vector1=copia elementi di posto dispari(indice pari) e somma elementi di posto pari(indice dispari)
//vector1 e vector2 devono avere gli elementi di posto dispari uguali e la loro lunghezza deve essere pari
	if(vector1==NULL || vector2==NULL || length1<=0 || length2<=0 || length1!=length2 || length1%2!=0 ){
		handle_error_with_exit("invalid parameter sum_elem_multiple_of_2\n");
	}
	mpz_t t;
	mpz_init(t);
	for(int i=0;i<length1;i=i+2){
		mpz_add(t,vector1[i+1],vector2[i+1]);
		mpz_set(vector1[i+1],t);
	}
	mpz_clear(t);
	return;
}

void find_max_array_mpz(mpz_t max,mpz_t*array,long length){
	if(array==NULL || max==NULL || length<=0){
		handle_error_with_exit("error in sort_array_number\n");
	}
	mpz_set(max,array[0]);
	if(length==1){
		return;
	}
	for(long i=1;i<length;i++){
		if(mpz_cmp(max,array[i])<0){//se il valore nell'array è maggiore prendi quello
			mpz_set(max,array[i]);
		}
	}
	return;
}


void divide_elem_multiple_of_2_by_x(mpz_t*vector,int length,double x){//vector = num1 exp1 num2 exp2 ecc,divide elementi di posto pari per 2
	//il vettore deve avere lunghezza pari
	if(vector==NULL || length<=0 || length%2!=0 || x==0){
		handle_error_with_exit("divide_elem_multiple_of_2_by_2\n");
	}
	if(x!=2){
		handle_error_with_exit("error in parameter\n");
	}
	for(int i=1;i<length;i=i+2){
		if(mpz_divisible_ui_p(vector[i],(unsigned long)x)==0){
			handle_error_with_exit("error in mpz_divisible,divide_elem_multiple_of_2_by_2\n");
		}
		mpz_divexact_ui(vector[i],vector[i],x);
	}
	return;
}
