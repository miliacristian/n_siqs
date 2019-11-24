#include "print.h"

char not_print_list(int length){
	if(length>=THRESOLD_PRINT_LIST){
		return 1;
	}
	return 0;
}
char not_print_array(long length){//ritorna 1 se non bisogna stampare 0 altrimenti
	if(length<=THRESOLD_PRINT_ARRAY){
		return 0;
	}
	return 1;
}
char not_print_matrix(long num_row,long num_col){
	if(num_row<=THRESOLD_PRINT_MATRIX && num_col<=THRESOLD_PRINT_MATRIX){
		return 0;
	}
	return 1;
}

void print_binary_array_unsigned_long(unsigned long*binary_array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(binary_array==NULL){
		printf("array is empty in print_binary_array\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		unsigned long num=binary_array[i];
		unsigned long temp=1;
		temp=temp<<63;
		unsigned long j=0;
		while (j<BIT_OF_UNSIGNED_LONG) {
			if (num & temp)
				printf("1");
			else
				printf("0");
			temp >>= 1;
			j++;
		}
	}
	printf("\n");
	return;
}

void print_binary_matrix(unsigned long**binary_matrix,int num_row,int num_col){
	if(binary_matrix==NULL || *binary_matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
		printf("riga %d:",i);
		print_binary_array_unsigned_long(binary_matrix[i],num_col);//stampa 1 riga
	}
	return;
}
void print_int(int a,char*string){
	if(string==NULL){
		handle_error_with_exit("error in print_int\n");
		return;
	}
	printf("%s=%d\n",string,a);
	return;
}
void print_char(char a,char*string){
	if(string==NULL){
		handle_error_with_exit("error in print_int\n");
		return;
	}
	printf("%s=%c\n",string,a);
	return;
}
void print_double(double a,char*string){
	if(string==NULL){
		handle_error_with_exit("error in print_int\n");
		return;
	}
	printf("%s=%lf\n",string,a);
	return;
}
void print_long(long a,char*string){
	if(string==NULL){
		handle_error_with_exit("error in print_int\n");
		return;
	}
	printf("%s=%ld\n",string,a);
	return;
}

void print_array_char(char*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_char\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	printf("\n");
	return;
}

void print_array_float(float*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_char\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%f ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_double(double*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_double\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_double\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%lf ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_long(long*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_long\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_long\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%ld ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_int(int*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_array_long\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_int\n");
		return;
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	printf("\n");
	return;
}
void print_vector_char(char*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_vector_char\n");
	}
	printf("start print\n");
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("  %d\n",array[i]);
	}
	printf("end print\n");
	return;
}
void print_vector_long(long*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_vector_long\n");
	}
	if(not_print_array(length)==1){
		return;
	}
	printf("start print\n");
	for(int i=0;i<length;i++){
		printf("  %ld\n",array[i]);
	}
	printf("end print\n");
	return;
}
void print_vector_int(int*array,int length){
	if(length<=0){
		handle_error_with_exit("error in print_vector_int\n");
	}
	if(not_print_array(length)==1){
		return;
	}
	printf("start print\n");
	for(int i=0;i<length;i++){
		printf("  %d\n",array[i]);
	}
	printf("end print\n");
	return;
}

void print_matrix_int(int**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			printf("riga %d:",i);
			print_array_int(matrix[i],num_col);//stampa 1 riga
	}
	return;
}
void print_matrix_long(long**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			printf("riga %d:",i);
			print_array_long(matrix[i],num_col);//stampa 1 riga
	}
	return;
}
void print_matrix_double(double**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_double\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			print_array_double(matrix[i],num_col);
	}
	return;
}
void print_matrix_char(char**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_double\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
		print_array_char(matrix[i],num_col);
	}
	return;
}


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



void print_linear_system(char*linear_system,int cardinality_factor_base,int num_of_B_smooth){//parametri matrice,righe,colonne
	if(linear_system==NULL){
		printf("linear system is empty\n");
	}
	if(num_of_B_smooth <=0 || cardinality_factor_base <=0){
		handle_error_with_exit("error in print_linear_system\n");
	}
	printf("linear_system\n");
	char*pointer=linear_system;
	for(int i=0;i<cardinality_factor_base;i++) {
		print_array_char(pointer, num_of_B_smooth);
		pointer = pointer + num_of_B_smooth;
	}
	return;
}

void print_base_linear_system(int**linear_system,int num_B_smooth,int dim_sol){
	if(linear_system==NULL || *linear_system==NULL){
		printf("linear system is empty\n");
	}
	if(num_B_smooth <=0 || dim_sol <=0){
		handle_error_with_exit("error in print_base_linear_system\n");
	}
	printf("base of linear system\n");
	print_matrix_int(linear_system,num_B_smooth, dim_sol);
	return;
}
void print_all_solution(int**matrix_solution,int dim_sol,int num_B_smooth){
	if(matrix_solution==NULL || *matrix_solution==NULL){
		printf("linear system is empty\n");
	}
	if(num_B_smooth <=0 || dim_sol <=0){
		handle_error_with_exit("error in print_all_solution\n");
	}
	printf("solution of system:\n");
	print_matrix_int(matrix_solution,(int)pow(2,dim_sol)-1,num_B_smooth);
	return;
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






