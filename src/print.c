#include "print.h"

char not_print_list(int length,long thresold){
	if(length<=thresold){
		return 0;
	}
	return 1;
}
char not_print_array(long length,long thresold){//ritorna 1 se non bisogna stampare 0 altrimenti
	if(length<=thresold){
		return 0;
	}
	return 1;
}
char not_print_matrix(long num_row,long num_col,long thresold_row,long thresold_col){
	if(num_row<=thresold_row && num_col<= thresold_col){
		return 0;
	}
	return 1;
}

void print_binary_array_unsigned_long(unsigned long*binary_array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(binary_array==NULL){
		printf("array is empty in print_binary_array\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
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

void print_binary_matrix(unsigned long**binary_matrix,int num_row,int num_col,long thresold_row,long thresold_col){
	if(binary_matrix==NULL || *binary_matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col,thresold_row,thresold_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
		printf("riga %d:",i);
		print_binary_array_unsigned_long(binary_matrix[i],num_col,thresold_col);//stampa 1 riga
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

void print_array_char(char*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_char\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	printf("\n");
	return;
}

void print_array_float(float*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_char\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_char\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%f ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_double(double*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_double\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_double\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%lf ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_long(long*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_long\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_long\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%ld ",array[i]);
	}
	printf("\n");
	return;
}
void print_array_int(int*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_array_long\n");
	}
	if(array==NULL){
		printf("array is empty in print_array_int\n");
		return;
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	printf("\n");
	return;
}
void print_vector_char(char*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_vector_char\n");
	}
	printf("start print\n");
	if(not_print_array(length,thresold)==1){
		return;
	}
	for(int i=0;i<length;i++){
		printf("  %d\n",array[i]);
	}
	printf("end print\n");
	return;
}
void print_vector_long(long*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_vector_long\n");
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	printf("start print\n");
	for(int i=0;i<length;i++){
		printf("  %ld\n",array[i]);
	}
	printf("end print\n");
	return;
}
void print_vector_int(int*array,int length,long thresold){
	if(length<=0){
		handle_error_with_exit("error in print_vector_int\n");
	}
	if(not_print_array(length,thresold)==1){
		return;
	}
	printf("start print\n");
	for(int i=0;i<length;i++){
		printf("  %d\n",array[i]);
	}
	printf("end print\n");
	return;
}

void print_matrix_int(int**matrix,int num_row,int num_col,long thresold_row,long thresold_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col,thresold_row,thresold_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			printf("riga %d:",i);
			print_array_int(matrix[i],num_col,thresold_col);//stampa 1 riga
	}
	return;
}
void print_matrix_long(long**matrix,int num_row,int num_col,long thresold_row,long thresold_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_int\n");
	}
	if(not_print_matrix(num_row,num_col,thresold_row,thresold_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			printf("riga %d:",i);
			print_array_long(matrix[i],num_col,thresold_col);//stampa 1 riga
	}
	return;
}
void print_matrix_double(double**matrix,int num_row,int num_col,long thresold_row,long thresold_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_double\n");
	}
	if(not_print_matrix(num_row,num_col,thresold_row,thresold_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			print_array_double(matrix[i],num_col,thresold_col);
	}
	return;
}
void print_matrix_char(char**matrix,int num_row,int num_col,long thresold_row,long thresold_col){
	if(matrix==NULL || *matrix==NULL){
		printf("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_double\n");
	}
	if(not_print_matrix(num_row,num_col,thresold_row,thresold_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
		print_array_char(matrix[i],num_col,thresold_col);
	}
	return;
}


void print_linear_system(char*linear_system,int num_rows,int num_cols,long thresold_col){//parametri matrice,righe,colonne
	if(linear_system==NULL){
		printf("linear system is empty\n");
	}
	if(num_cols <=0 || num_rows <=0){
		handle_error_with_exit("error in print_linear_system\n");
	}
	printf("linear_system\n");
	char*pointer=linear_system;
	for(int i=0;i<num_rows;i++) {
		print_array_char(pointer, num_cols,thresold_col);
		pointer = pointer + num_cols;
	}
	return;
}

void print_base_linear_system(int**linear_system,int num_rows,int num_cols,long thresold_row,long thresold_col){
	if(linear_system==NULL || *linear_system==NULL){
		printf("linear system is empty\n");
		return;
	}
	if(num_rows <=0 || num_cols <=0){
		handle_error_with_exit("error in print_base_linear_system\n");
	}
	printf("base of linear system\n");
	print_matrix_int(linear_system,num_rows, num_cols,thresold_row,thresold_col);
	return;
}







