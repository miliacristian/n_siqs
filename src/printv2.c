#include "printv2.h"

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







