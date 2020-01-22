
#include "matrix_functions.h"


char* from_matrix_binary_to_matrix_char(unsigned long**binary_linear_system,int num_row,int num_col_binary_matrix,int*num_col_linear_system){
	if(binary_linear_system==NULL || *binary_linear_system==NULL || num_row<=0 || num_col_binary_matrix<=0 || num_col_linear_system==NULL){
		handle_error_with_exit("error in from_matrix_binary_to_matrix_char\n");
	}
	char*linear_system;
	int index_col_linear_system;
	int index;
	char bit;
	*num_col_linear_system=num_col_binary_matrix*BIT_OF_UNSIGNED_LONG;
	linear_system=alloc_array_char(num_row*(*num_col_linear_system));
	for(int i=0;i<num_row;i++){
		for(int j=0;j<num_col_binary_matrix;j++){
			for(int b=0;(unsigned long)b<BIT_OF_UNSIGNED_LONG;b++){
				bit=(binary_linear_system[i][j]>> (BIT_OF_UNSIGNED_LONG-b-1)) & 1U;
				index_col_linear_system=BIT_OF_UNSIGNED_LONG*j+b;
				index=get_index(i,index_col_linear_system,*num_col_linear_system);
				linear_system[index]=bit;
			}
		}
	}
	return linear_system;
}
void copy_matrix_with_array(char**linear_system2,char*linear_system,int num_row,int num_col){
    if(linear_system2==NULL || *linear_system2==NULL || linear_system==NULL || num_row<=0 || num_col<=0){
        handle_error_with_exit("error in copy_matrix_with_array\n");
    }
    long index;
    for(int i=0;i<num_row;i++){
        for(int j=0;j<num_col;j++){
            index=get_index(i,j,num_col);
            linear_system2[i][j]=linear_system[index];
        }
    }
    return;
}

pthread_t *alloc_array_tid(int num_thread){
	if(num_thread<0){
		handle_error_with_exit("error in create_thread\n");
	}
	if(num_thread==0){
		return NULL;
	}
	pthread_t*array_tid=malloc(sizeof(pthread_t)*num_thread);
	if(array_tid==NULL){
		handle_error_with_exit("error in malloc alloc array tid\n");
	}
	return array_tid;
}
pthread_t *alloc_array_tid_and_materialize(int num_thread){
	if(num_thread<0){
		handle_error_with_exit("error in create_thread\n");
	}
	if(num_thread==0){
		return NULL;
	}
	pthread_t*array_tid=malloc(sizeof(pthread_t)*num_thread);
	if(array_tid==NULL){
		handle_error_with_exit("error in malloc alloc array tid\n");
	}
	memset(array_tid,0,sizeof(pthread_t)*num_thread);
	return array_tid;
}

char* alloc_array_char(long length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_char\n");
	}
	char *array=malloc(sizeof(char)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc\n");
	}
	return array;
}
char* alloc_array_char_and_materialize(long length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_char\n");
	}
	char *array=malloc(sizeof(char)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc\n");
	}
	memset(array,0,sizeof(char)*length);
	return array;
}
long* alloc_array_long(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc aray_long\n");
	}
	long *array=malloc(sizeof(long)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array long\n");
	}
	return array;
}
long* alloc_array_long_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc aray_long\n");
	}
	long *array=malloc(sizeof(long)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array long\n");
	}
	memset(array,0,sizeof(long)*length);
	return array;
}
unsigned long* alloc_array_unsigned_long(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc aray_long\n");
	}
	unsigned long *array=malloc(sizeof(unsigned long)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array long\n");
	}
	return array;
}

unsigned long* alloc_array_unsigned_long_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc aray_long\n");
	}
	unsigned long *array=malloc(sizeof(unsigned long)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array long\n");
	}
	memset(array,0,sizeof(unsigned long)*length);
	return array;
}
int* alloc_array_int(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array int\n");
	}
	int *array=malloc(sizeof(int)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array int\n");
	}
	memset(array,0,sizeof(int)*length);
	return array;
}
int* alloc_array_int_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array int\n");
	}
	int *array=malloc(sizeof(int)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array int\n");
	}
	memset(array,0,sizeof(int)*length);
	return array;
}

double* alloc_array_double(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array double\n");
	}
	double *array=malloc(sizeof(double)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array double\n");
	}
	return array;
}
double* alloc_array_double_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array double\n");
	}
	double *array=malloc(sizeof(double)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array double\n");
	}
	memset(array,0,sizeof(double)*length);
	return array;
}
float* alloc_array_float(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array double\n");
	}
	float *array=malloc(sizeof(float)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array double\n");
	}
	return array;
}
float* alloc_array_float_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc array double\n");
	}
	float *array=malloc(sizeof(float)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array double\n");
	}
	memset(array,0,sizeof(float)*length);
	return array;
}
long** alloc_array_pointer_to_long(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to long\n");
	}
	long **array=malloc(sizeof(long*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to long\n");
	}
	return array;
}
long** alloc_array_pointer_to_long_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to long\n");
	}
	long **array=malloc(sizeof(long*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to long\n");
	}
	memset(array,0,sizeof(long*)*(length));
	return array;
}
unsigned long** alloc_array_pointer_to_unsigned_long(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to unsigned long\n");
	}
	unsigned long **array=malloc(sizeof(unsigned long*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to unsigned long\n");
	}
	return array;
}
unsigned long** alloc_array_pointer_to_unsigned_long_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to unsigned long\n");
	}
	unsigned long **array=malloc(sizeof(unsigned long*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to unsigned long\n");
	}
	memset(array,0,sizeof(unsigned long*)*(length));
	return array;
}

int** alloc_array_pointer_to_int(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc pointer to int\n");
	}
	int **array=malloc(sizeof(int*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to int\n");
	}
	memset(array,0,sizeof(int*)*(length));
	return array;
}
int** alloc_array_pointer_to_int_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc pointer to int\n");
	}
	int **array=malloc(sizeof(int*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to int\n");
	}
	memset(array,0,sizeof(int*)*(length));
	return array;
}


double** alloc_array_pointer_to_double(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_pointer_to_double\n");
	}
	double **array=malloc(sizeof(double*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to double\n");
	}
	return array;
}
double** alloc_array_pointer_to_double_and_materialize(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_pointer_to_double\n");
	}
	double **array=malloc(sizeof(double*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to double\n");
	}
	memset(array,0,sizeof(double*)*(length));
	return array;
}
char** alloc_array_pointer_to_char(int length){
    if(length<=0){
        handle_error_with_exit("error in parameter alloc_array_pointer_to_char\n");
    }
    char **array=malloc(sizeof(char*)*(length));
    if(array==NULL){
        handle_error_with_exit("error in malloc alloc array pointer to char\n");
    }
    return array;
}
char** alloc_array_pointer_to_char_and_materialize(int length){
    if(length<=0){
        handle_error_with_exit("error in parameter alloc_array_pointer_to_char\n");
    }
    char **array=malloc(sizeof(char*)*(length));
    if(array==NULL){
        handle_error_with_exit("error in malloc alloc array pointer to char\n");
    }
    memset(array,0,sizeof(char*)*(length));
    return array;
}
int** alloc_matrix_int(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_int\n");
	}
	int**matrix;
	matrix=alloc_array_pointer_to_int(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_int(num_col);
	}
	return matrix;
}
int** alloc_matrix_int_and_materialize(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_int\n");
	}
	int**matrix;
	matrix=alloc_array_pointer_to_int_and_materialize(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_int_and_materialize(num_col);
	}
	return matrix;
}


double** alloc_matrix_double(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_double\n");
	}
	double**matrix;
	matrix=alloc_array_pointer_to_double(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_double(num_col);
	}
	return matrix;
}
double** alloc_matrix_double_and_materialize(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_double\n");
	}
	double**matrix;
	matrix=alloc_array_pointer_to_double_and_materialize(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_double_and_materialize(num_col);
	}
	return matrix;
}
char** alloc_matrix_char(int num_row,int num_col){
    if(num_row<=0 || num_col <=0){
        handle_error_with_exit("error in alloc_matrix_char\n");
    }
    char**matrix;
    matrix=alloc_array_pointer_to_char(num_row);
    for(int i=0;i<num_row;i++){
        matrix[i]=alloc_array_char(num_col);
    }
    return matrix;
}
char** alloc_matrix_char_and_materialize(int num_row,int num_col){
    if(num_row<=0 || num_col <=0){
        handle_error_with_exit("error in alloc_matrix_char\n");
    }
    char**matrix;
    matrix=alloc_array_pointer_to_char_and_materialize(num_row);
    for(int i=0;i<num_row;i++){
        matrix[i]=alloc_array_char_and_materialize(num_col);
    }
    return matrix;
}


long** alloc_matrix_long(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_long\n");
	}
	long**matrix;
	matrix=alloc_array_pointer_to_long(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_long(num_col);
	}
	return matrix;
}
long** alloc_matrix_long_and_materialize(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_long\n");
	}
	long**matrix;
	matrix=alloc_array_pointer_to_long_and_materialize(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_long_and_materialize(num_col);
	}
	return matrix;
}
unsigned long** alloc_matrix_unsigned_long(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_unsigned_long\n");
	}
	unsigned long**matrix;
	matrix=alloc_array_pointer_to_unsigned_long(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_unsigned_long(num_col);
	}
	return matrix;
}
unsigned long** alloc_matrix_unsigned_long_and_materialize(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_unsigned_long\n");
	}
	unsigned long**matrix;
	matrix=alloc_array_pointer_to_unsigned_long_and_materialize(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_unsigned_long_and_materialize(num_col);
	}
	return matrix;
}


char is_in_array_long(long*array,long length,long p_i){//ritorna vero se p_i è presente nell'array,0 altrimenti
	if(array==NULL || length<=0){
		handle_error_with_exit("error in is_in_array_long\n");
	}
	for(long i=0;i<length;i++){
		if(array[i]==p_i){
			return 1;
		}
	}
	return 0;
}
char is_in_array_int(int*array,long length,long p_i){//ritorna vero se p_i è presente nell'array,0 altrimenti
	if(array==NULL || length<=0){
		handle_error_with_exit("error in is_in_array_long\n");
	}
	for(long i=0;i<length;i++){
		if(array[i]==p_i){
			return 1;
		}
	}
	return 0;
}


void divide_vector_multiple_of_2_by_2(int*vector,int length){
    //ogni elemento deve essere multiplo di 2
    if(vector==NULL || length<=0){
        handle_error_with_exit("divide_elem_by_2\n");
    }
    for(int i=0;i<length;i++){
        if(vector[i]==0){
            continue;
        }
        if((vector[i] & 1)!=0){
            handle_error_with_exit("error in divide_elem_by_2\n");
        }
        vector[i]=vector[i]/2;
    }
    return;
}
int** alloc_linear_system(int num_rows,int num_cols){//alloca sistema lineare
	if(num_rows<=0 || num_cols<=0){
		handle_error_with_exit("invalid parameter alloc linear system\n");
	}
	int**linear_system;
	linear_system=alloc_matrix_int(num_rows,num_cols);
	return linear_system;
}
int* get_coli(int **matrix,int num_row,int num_col,int index_col){//indice parte da 0,ottiene la colonna iesima della matrice
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || index_col<0){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	if(index_col>=num_col){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	int*col=alloc_array_int(num_row);
	for(int i=0;i<num_row;i++){
		col[i]=matrix[i][index_col];//index_col è fissato perché la colonna rimane la stessa durante il ciclo
	}
	return col;
}

#if DEBUG==1
char check_if_array_var_is_valid(int **matrix_linear_system,int num_row,int num_col,char*array_var){
	//verifica che l'array delle varaibili è stato calcolato correttamente
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0 || array_var==NULL){
		handle_error_with_exit("invalid parameter check_if_array_var_is_valid\n");
	}
	int count=0;
	char*array_var_temp=alloc_array_char(num_col);
	memcpy(array_var_temp,array_var,sizeof(char)*num_col);
	int num_row_not_null=count_rows_not_null(matrix_linear_system,num_row,num_col);
	for(int i=0;i<num_row_not_null;i++){
		count=0;
		for(int j=0;j<num_col;j++){
			if(matrix_linear_system[num_row_not_null-1-i][num_col-1-j]==0){
				continue;
			}
			else if(array_var_temp[num_col-1-j]==2){//se variabile ricavabile trovata marcala come calcolata,cosi alla prossima riga 					non viene ricontata,infatti vogliamo calcolare le variabili ricavabili per ogni colonna
				array_var_temp[num_col-1-j]=3;
				count++;
			}
			if(count>1){
				return 0;
			}
		}
		if(count!=1){
			return 0;
		}
	}
	free(array_var_temp);
	return 1;
}
char check_if_array_var_is_valid_char(char*matrix_linear_system,int num_row,int num_col,char*array_var){
	//verifica che l'array delle varaibili è stato calcolato correttamente
	if(matrix_linear_system==NULL || num_row<=0 || num_col<=0 || array_var==NULL){
		handle_error_with_exit("invalid parameter check_if_array_var_is_valid\n");
	}
	int count=0;
	char*array_var_temp=alloc_array_char(num_col);
	memcpy(array_var_temp,array_var,sizeof(char)*num_col);
	int num_row_not_null=count_rows_not_null_char(matrix_linear_system,num_row,num_col);
	for(int i=0;i<num_row_not_null;i++){
		count=0;
		for(int j=0;j<num_col;j++){
			long index=get_index(num_row_not_null-1-i,num_col-1-j,num_col);
			if(matrix_linear_system[index]==0){
				continue;
			}
			else if(array_var_temp[num_col-1-j]==2){//se variabile ricavabile trovata marcala come calcolata,cosi alla prossima riga 					non viene ricontata,infatti vogliamo calcolare le variabili ricavabili per ogni colonna
				array_var_temp[num_col-1-j]=3;
				count++;
			}
			if(count>1){
				return 0;
			}
		}
		if(count!=1){
			return 0;
		}
	}
	free(array_var_temp);
	return 1;
}
char verify_solution(int **matrix_linear_system,int num_row,int num_col,int*solution){
	//verifica che la soluzione del sistema lineare è corretta
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0 || solution==NULL){
		handle_error_with_exit("invalid parameter verify_solution\n");
	}
	int count=0;
	int num_row_not_null=count_rows_not_null(matrix_linear_system,num_row,num_col);
	for(int i=0;i<num_row_not_null;i++){
		count=0;
		for(int j=0;j<num_col;j++){
			if(matrix_linear_system[i][j]==0){
				continue;
			}
			else if(matrix_linear_system[i][j]==1){
				if(solution[j]==1){
					count++;
				}
				else if(solution[j]==0){
				}
				else{
					handle_error_with_exit("invalid solution\n");
				}
			}
			else{
				handle_error_with_exit("invalid matrix_linear_system\n");
			}
		}
		if((count & 1)!=0){
			return 0;
		}
	}
	return 1;
}
char verify_solution_char(char*matrix_linear_system,int num_row,int num_col,const int*solution){
    //verifica che la soluzione del sistema lineare è corretta
    if(matrix_linear_system==NULL || num_row<=0 || num_col<=0 || solution==NULL){
        handle_error_with_exit("invalid parameter verify_solution\n");
    }
    int count=0;
    int num_row_not_null=count_rows_not_null_char(matrix_linear_system,num_row,num_col);
    for(int i=0;i<num_row_not_null;i++){
        count=0;
        for(int j=0;j<num_col;j++){
            long index=get_index(i,j,num_col);
            if(matrix_linear_system[index]==0){
                continue;
            }
            else if(matrix_linear_system[index]==1){
                if(solution[j]==1){
                    count++;
                }
                else if(solution[j]==0){
                }
                else{
                    handle_error_with_exit("invalid solution\n");
                }
            }
            else{
                handle_error_with_exit("invalid matrix_linear_system\n");
            }
        }
        if((count & 1)!=0){
            return 0;
        }
    }
    return 1;
}
char check_solution_base_matrix(int**linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base){
	//verifica che le soluzioni provenienti dai vettori di base sono soluzioni ammissibili
	if(linear_system==NULL || *linear_system==NULL || num_row_system<=0 || num_col_system<=0 || base_matrix==NULL || *base_matrix==NULL ||
	num_row_base<=0 || num_col_base<=0){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	char test=0;
	for(int i=0;i<num_col_base;i++){
		int *column=get_coli(base_matrix,num_row_base,num_col_base,i);
		test=verify_solution(linear_system,num_row_system,num_col_system,column);
		if(test==0){
			return 0;
		}
		free(column);
		column=NULL;
	}
	return 1;
}
char check_solution_base_matrix_char(char*linear_system,int num_row_system,int num_col_system,int **base_matrix,int num_row_base,int num_col_base){
    //verifica che le soluzioni provenienti dai vettori di base sono soluzioni ammissibili
    if(linear_system==NULL || num_row_system<=0 || num_col_system<=0 || base_matrix==NULL || *base_matrix==NULL ||
       num_row_base<=0 || num_col_base<=0){
        handle_error_with_exit("error in parameter reduce echelon form\n");
    }
    char test=0;
    for(int i=0;i<num_col_base;i++){
        int *column=get_coli(base_matrix,num_row_base,num_col_base,i);
        test=verify_solution_char(linear_system,num_row_system,num_col_system,column);
        if(test==0){
            return 0;
        }
        free(column);
        column=NULL;
    }
    return 1;
}
char check_if_matrix_is_echelon_reduce(int**matrix_linear_system,int num_row,int num_col){
	//verifica che la matrice è ridotta a scala,non verifica che la matrice è rref
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in check_if_matrix_is_echelon_reduce\n");
	}
	int i=0;
	int j=0;
	while(i<num_row && j<num_col){
		if(matrix_linear_system[i][j]!=0){//un pivot è stato trovato verifica che tutta la colonna ha elementi nulli
			for(int temp=i+1;temp<num_row;temp++){//scansiona tutta la colonna partendo dall'elemento successivo
				if(matrix_linear_system[temp][j]!=0){//se l'elemento è diverso da 0 non è ridotta a scala
					return 0;
				}
			}
			i+=1;
			j+=1;
		}
		else{//matrix[i][j]==0,verifica che tutta la colonna è zero poi aumenta di uno solo j
			for(int temp=i+1;temp<num_row;temp++){//scansiona tutta la colonna partendo dall'elemento successivo
				if(matrix_linear_system[temp][j]!=0){//se l'elemento è diverso da 0 non è ridotta a scala
					return 0;
				}
			}
			j+=1;
		}
	}
	return 1;
}
char check_if_matrix_char_is_echelon_reduce(char*matrix_linear_system,int num_row,int num_col){
	//verifica che la matrice è ridotta a scala,non verifica che la matrice è rref
	if(matrix_linear_system==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in check_if_matrix_is_echelon_reduce\n");
	}
	int i=0;
	int j=0;
	while(i<num_row && j<num_col){
		long index=get_index(i,j,num_col);
		if(matrix_linear_system[index]!=0){//un pivot è stato trovato verifica che tutta la colonna ha elementi nulli
			for(int temp=i+1;temp<num_row;temp++){//scansiona tutta la colonna partendo dall'elemento successivo
				long index=get_index(temp,j,num_col);
				if(matrix_linear_system[index]!=0){//se l'elemento è diverso da 0 non è ridotta a scala
					return 0;
				}
			}
			i+=1;
			j+=1;
		}
		else{//matrix[i][j]==0,verifica che tutta la colonna è zero poi aumenta di uno solo j
			for(int temp=i+1;temp<num_row;temp++){//scansiona tutta la colonna partendo dall'elemento successivo
				long index=get_index(temp,j,num_col);
				if(matrix_linear_system[index]!=0){//se l'elemento è diverso da 0 non è ridotta a scala
					return 0;
				}
			}
			j+=1;
		}
	}
	return 1;
}
char check_if_array_var_is_correct(char*array_var,int length,int dim_sol){
	//il numero delle variabili libre deve essere uguale alla dimensione della soluzione
	if(array_var==NULL || length<=0 || dim_sol<=0){//la dimensione non può essere 0 altrimenti ho una sola soluzione
		handle_error_with_exit("error in check if array_var is correct\n");
	}
	int count_free_var=0;//conta le variabili libere
	for(int i=0;i<length;i++){
		if(array_var[i]!=1 && array_var[i]!=2 && array_var[i]!=3){
			handle_error_with_exit("error in check if array_var is correct cycle\n");
		}
		if(array_var[i]==1){
			count_free_var++;//variabile libera trovata
		}
	}
	if(count_free_var!=dim_sol){
		return 0;//se sono diversi ritorna 0
	}
	return 1;//tutto ok
}
char check_if_matrix_is_reduce_mod_n(int**matrix,int num_row,int num_col,int n){
//verifica che la matrice è ridotta modulo n
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || n<=0){
		printf("num_row=%d num_col=%d n=%d\n",num_row,num_col,n);
		handle_error_with_exit("error in check if matrix is reduce mod n\n");
	}
	for(int i=0;i<num_row;i++){
		for(int j=0;j<num_col;j++){
			if(matrix[i][j]>=n || matrix[i][j]<0){
				return 0;
			}
		}
	}
	return 1;
}
char check_if_matrix_char_is_reduce_mod_n(char*matrix,int num_row,int num_col,int n){
//verifica che la matrice è ridotta modulo n
    if(matrix==NULL || num_row<=0 || num_col<=0 || n<=0){
        printf("num_row=%d num_col=%d n=%d\n",num_row,num_col,n);
        handle_error_with_exit("error in check if matrix char is reduce mod n\n");
    }
    long index;
    for(int i=0;i<num_row;i++){
        for(int j=0;j<num_col;j++){
            index=get_index(i,j,num_col);
            if(matrix[index]>=n || matrix[index]<0){
                return 0;
            }
        }
    }
    return 1;
}
char check_if_array_is_reduce_mod_n(int*array,int length,int n){
//verifica che l'array è ridotto modulo n
	if(array==NULL || length<=0 || n<=0){
		handle_error_with_exit("error in check if array is reduce mod n\n");
	}
	for(int i=0;i<length;i++){
		if(array[i]>=n || array[i]<0){
			return 0;
		}
	}
	return 1;
}
#endif
void swap_row(int**matrix,int num_row,int num_col,int ind_row1,int ind_row2){//ind_row1,int ind_row2 start at 0
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || ind_row1==ind_row2 || ind_row1>=num_row || ind_row2>=num_row || ind_row1<0 || ind_row2<0){
		handle_error_with_exit("error in parameter swap_row\n");
	}
	int temp;
	for(int i=0;i<num_col;i++){//il numero di colonne corrisponde alla lunghezza di una riga
		temp=matrix[ind_row1][i];
		matrix[ind_row1][i]=matrix[ind_row2][i];
		matrix[ind_row2][i]=temp;
	}
	return;
}
void swap_row_char2(char**matrix,int num_row,int num_col,int ind_row1,int ind_row2){//ind_row1,int ind_row2 start at 0
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || ind_row1==ind_row2 || ind_row1>=num_row || ind_row2>=num_row || ind_row1<0 || ind_row2<0){
		handle_error_with_exit("error in parameter swap_row\n");
	}
	char*temp=matrix[ind_row1];
	matrix[ind_row1]=matrix[ind_row2];
	matrix[ind_row2]=temp;
	return;
}
void swap_row_unsigned_long(unsigned long**matrix,int num_row,int ind_row1,int ind_row2){//ind_row1,int ind_row2 start at 0
	if(matrix==NULL || *matrix==NULL || num_row<=0 || ind_row1==ind_row2 || ind_row1>=num_row || ind_row2>=num_row || ind_row1<0 || ind_row2<0){
		handle_error_with_exit("error in parameter swap_row\n");
	}
	unsigned long*temp=matrix[ind_row1];
	matrix[ind_row1]=matrix[ind_row2];
	matrix[ind_row2]=temp;
	return;
}
void swap_row_char(char*matrix,int num_row,int num_col,int ind_row1,int ind_row2){//ind_row1,int ind_row2 start at 0
    if(matrix==NULL || num_row<=0 || num_col<=0 || ind_row1==ind_row2 || ind_row1>=num_row || ind_row2>=num_row || ind_row1<0 || ind_row2<0){
        handle_error_with_exit("error in parameter swap_row\n");
    }
    int temp;
    long index1,index2;
    for(int i=0;i<num_col;i++){//il numero di colonne corrisponde alla lunghezza di una riga
        index1=get_index(ind_row1,i,num_col);//moltiplicazione
        index2=get_index(ind_row2,i,num_col);//moltiplicazione
        temp=matrix[index1];
        matrix[index1]=matrix[index2];
        matrix[index2]=temp;
    }
    return;
}
void reduce_matrix_mod_n(int**matrix,int num_row,int num_col,int n){//riduce gli elementi di una matrice a modulo n
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || n<=0){
		handle_error_with_exit("error in reduce matrix_mod_n\n");
	}
	for(int i=0;i<num_row;i++){
		for(int j=0;j<num_col;j++){
			matrix[i][j]=calculate_mod_n((long)matrix[i][j],(long)n);
		}
	}
	return;
}
void reduce_array_mod_n(int*array,int length,int n){//riduce gli elementi di un array modulo n
	if(array==NULL || length<=0 || n<=0){
		handle_error_with_exit("error in reduce matrix_mod_n\n");
	}
	for(int i=0;i<length;i++){
			array[i]=calculate_mod_n((long)array[i],(long)n);
	}
	return;
}
int count_rows_not_null(int **matrix,int num_row,int num_col){//la matrice deve essere ridotta a scala,conta righe non nulle sulla matrice ridotta a scala
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	int count=0;//conta righe non nulle
	for(int i=num_row-1;i>=0;i--){//ciclo inverso per le righe,scansiona prima le ultime righe poi le prime
		for(int j=0;j<num_col;j++){//ciclo colonne
			if(matrix[i][j]!=0){//la riga i-esima non è nulla aumenta il contatore e passa alla prossima riga
				count++;
				break;
			}
		}
	}
	return count;
}
int count_rows_not_null_char(char*matrix,int num_row,int num_col){//la matrice deve essere ridotta a scala,conta righe non nulle sulla matrice ridotta a scala
    if(matrix==NULL || num_row<=0 || num_col<=0){
        handle_error_with_exit("error in parameter get_coli\n");
    }
    int count=0;//conta righe non nulle
    for(int i=num_row-1;i>=0;i--){//ciclo inverso per le righe,scansiona prima le ultime righe poi le prime
        for(int j=0;j<num_col;j++){//ciclo colonne
            long index=get_index(i,j,num_col);
            if(matrix[index]!=0){//la riga i-esima non è nulla aumenta il contatore e passa alla prossima riga
                count++;
                break;
            }
        }
    }
    return count;
}
int count_cols_not_null(int **matrix,int num_row,int num_col){//la matrice deve essere ridotta a scala,conta righe non nulle sulla matrice ridotta a scala
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	int count=0;//conta colonne non nulle
	for(int j=num_col-1;j>=0;j--){//ciclo inverso per le colonne,scansiona prima le ultime colonne poi le prime
		for(int i=0;i<num_row;i++){//ciclo righe
			if(matrix[i][j]!=0){//la colonna i-esima non è nulla aumenta il contatore e passa alla prossima colonna
				count++;
				break;
			}
		}
	}
	return count;
}
int count_cols_not_null_char(char*matrix,int num_row,int num_col){//la matrice deve essere ridotta a scala,conta righe non nulle sulla matrice ridotta a scala
	if(matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	int count=0;//conta colonne non nulle
	for(int j=num_col-1;j>=0;j--){//ciclo inverso per le colonne,scansiona prima le ultime colonne poi le prime
		for(int i=0;i<num_row;i++){//ciclo righe
			long index=get_index(i,j,num_col);
			if(matrix[index]!=0){//la colonna i-esima non è nulla aumenta il contatore e passa alla prossima colonna
				count++;
				break;
			}
		}
	}
	return count;
}
int calculate_dim_sol(int**matrix,int num_row,int num_col){//matrice deve essere ridotta a scala
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_col\n");
	}
	int num_row_not_null=count_rows_not_null(matrix,num_row,num_col);
	int num_col_not_null=count_cols_not_null(matrix,num_row,num_col);
	return num_col_not_null-num_row_not_null;//numero incognite-rango matrice
}
int calculate_dim_solution(int num_row_not_null,int num_col_not_null){
	return num_col_not_null-num_row_not_null;
}
int calculate_dim_sol_char(char*matrix,int num_row,int num_col,int*num_row_not_null,int*num_col_not_null){//matrice deve essere ridotta a scala
    if(matrix==NULL || num_row<=0 || num_col<=0 || num_row_not_null==NULL || num_col_not_null==NULL){
        handle_error_with_exit("error in parameter get_col\n");
    }
    *num_row_not_null=count_rows_not_null_char(matrix,num_row,num_col);
    *num_col_not_null=count_cols_not_null_char(matrix,num_row,num_col);
    return *num_col_not_null-*num_row_not_null;//numero incognite-rango matrice
}
int* sum_vector(int*vector1,int*vector2,int length1,int length2){//somma 2 vettori
	if(vector1==NULL || vector2==NULL || length1<=0 || length2<=0 || length1!=length2){
		handle_error_with_exit("invalid parameter sum_vector\n");
	}
	int*vector_result;
	vector_result=alloc_array_int(length1);
	for(int i=0;i<length1;i++){
		vector_result[i]=vector1[i]+vector2[i];
	}
	return vector_result;
}
long* sum_elem_multiple_of_2(long*vector1,long*vector2,int length1,int length2){//vettore risultato=copia elementi di posto dispari(indice pari) e somma elementi di posto pari(indice dispari)
//vector1 e vector2 devono avere gli elementi di posto dispari uguali e la loro lunghezza deve essere pari
	if(vector1==NULL || vector2==NULL || length1<=0 || length2<=0 || length1!=length2 || length1%2!=0 ){
		handle_error_with_exit("invalid parameter sum_elem_multiple_of_2\n");
	}
	long*vector_result=alloc_array_long(length1);
	for(int i=0;i<length1;i=i+2){//i è sempre pari al massimo vale length1-2
		vector_result[i]=vector1[i];//non c'entra nulla con sum_elem_mult 2
		vector_result[i+1]=vector1[i+1]+vector2[i+1];
	}
	return vector_result;
}


void prod_vector_and_scalar(int*vector,int scalar,int length){//prodotto vettore e scalare
	if(vector==NULL || length<=0){
		handle_error_with_exit("invalid parameter prod vector and scalar\n");
	}
	for(int i=0;i<length;i++){
		vector[i]=vector[i]*scalar;
	}
	return;
}

int* prod_vector_and_scalar_v2(int*vector,int scalar,int length){//prodotto vettore*scalare,ritorna vettore risultato
	if(vector==NULL || length<=0){
		handle_error_with_exit("invalid parameter prod vector and scalar v2\n");
	}
	int*result=alloc_array_int(length);
	for(int i=0;i<length;i++){
		result[i]=vector[i]*scalar;
	}
	return result;
}




void calculate_vector_base(int **matrix_linear_system,int num_row,int num_col,char*array_var,int*v,long thresold){//calcola un vettore di base soluzione del sistema lineare,matrix linear system è ridotta a scala e binaria,calcola vettori di base binari
	//si concentra sulla sottomatrice con righe non nulle!
	//memorizza la soluzione in v
	//per ogni riga le variabili libere hanno valore 0 o 1 (già assegnato)->ogni riga ha esattamente un valore ricavabile da calcolare
	//per calcolare il valore ricavabile:si vede qual'è il valore ricavabile della riga e si sommano i valori noti tra di loro e si ottiene:
	// coeff*var_ricavabile+somma_valori_noti=0 <-> coeff*var_ricavabile=-somma_valori_noti <-> var_ricavabile=-somma_valori_noti/coeff
	//una volta calcolato il valore della variabile viene impostata la variabile nell'array_var a calcolata -> array_var[i]=3
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0 || array_var==NULL || v==NULL){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	int index_to_calculate=-1;//variabile da calcolare nella riga
	int num_row_not_null=count_rows_not_null(matrix_linear_system,num_row,num_col);
	int rows_null=num_row-num_row_not_null;
	int count=0;//conta in ogni riga il numero di variabili libere o assegnate che hanno valore 1
	for(int i=0;i<num_row-rows_null;i++){//ciclo per tutta la sottomatrice con righe non nulle
		count=0;
		index_to_calculate=-1;
		for(int j=0;j<num_col;j++){//per ogni riga trova il valore da calcolare
			if(matrix_linear_system[num_row_not_null-1-i][num_col-1-j]==0){//se l'elemento della riga non nulla è 0 passa al prossimo 					elemento della riga
				continue;
			}
			if(array_var[num_col-1-j]==1 || array_var[num_col-1-j]==3){//se liberi o già calcolati(entrambi sono già stati assegnati)
				if(v[num_col-1-j]==1){//il valore della variabile libera o assgnata è 1
					count++;
				}
				else if(v[num_col-1-j]==0){
					//nothing
				}
				else{
					handle_error_with_exit("error in calculate vector base v[num_col-1-j]\n");
				}
			}
			else if(array_var[num_col-1-j]==2){//array_var[num_col-1-j]!=1 || array_var[num_col-1-j]!=3 ->array_var==2==ricavabile
				index_to_calculate=num_col-1-j;//memorizza la colonna in cui c'è la variabile da calcolare,
			}
			else{
				handle_error_with_exit("error in calculate vector base\n");
			}
		}
		if(index_to_calculate<0){//riga in cui non è stato trovato nessun elemento da calcolare
			printf("index_to_calculate=%d\n",index_to_calculate);
			print_array_char(array_var,num_col,thresold);
			handle_error_with_exit("error in index_calculate index<0\n");
		}
		if(index_to_calculate>num_col){
			printf("num_col=%d,index_to_calculate=%d\n",num_col,index_to_calculate);
			handle_error_with_exit("error in index_calculate,index>num_col\n");
		}
		if(num_row_not_null-1-i>num_row || num_row_not_null-1-i<0){
			handle_error_with_exit("error in index calculate vector base\n");
		}
		if((count & 1)==0){//divisibile per 2
			v[index_to_calculate]=0;
		}
		else{
			v[index_to_calculate]=1;
		}
		array_var[index_to_calculate]=3;//imposta variabile come calcolata,piano piano le variabili diventano tutte calcolate
	}
	return;
}
void calculate_vector_base_char(char*matrix_linear_system,int num_row,int num_col,char*array_var,int*v,int num_row_not_null,long thresold){//calcola un vettore di base soluzione del sistema lineare,matrix linear system è ridotta a scala e binaria,calcola vettori di base binari
    //si concentra sulla sottomatrice con righe non nulle!
    //memorizza la soluzione in v
    //per ogni riga le variabili libere hanno valore 0 o 1 (già assegnato)->ogni riga ha esattamente un valore ricavabile da calcolare
    //per calcolare il valore ricavabile:si vede qual'è il valore ricavabile della riga e si sommano i valori noti tra di loro e si ottiene:
    // coeff*var_ricavabile+somma_valori_noti=0 <-> coeff*var_ricavabile=-somma_valori_noti <-> var_ricavabile=-somma_valori_noti/coeff
    //una volta calcolato il valore della variabile viene impostata la variabile nell'array_var a calcolata -> array_var[i]=3
    if(matrix_linear_system==NULL || num_row<=0 || num_col<=0 || array_var==NULL || v==NULL){
        handle_error_with_exit("error in parameter get_coli\n");
    }
    int index_to_calculate=-1;//variabile da calcolare nella riga
    int rows_null=num_row-num_row_not_null;
    int count=0;//conta in ogni riga il numero di variabili libere o assegnate che hanno valore 1
    for(int i=0;i<num_row-rows_null;i++){//ciclo per tutta la sottomatrice con righe non nulle
        count=0;
        index_to_calculate=-1;
        for(int j=0;j<num_col;j++){//per ogni riga trova il valore da calcolare
            long index=get_index(num_row_not_null-1-i,num_col-1-j,num_col);
            if(matrix_linear_system[index]==0){//se l'elemento della riga non nulla è 0 passa al prossimo 					elemento della riga
                continue;
            }
            if(array_var[num_col-1-j]==1 || array_var[num_col-1-j]==3){//se liberi o già calcolati(entrambi sono già stati assegnati)
                if(v[num_col-1-j]==1){//il valore della variabile libera o assgnata è 1
                    count++;
                }
                else if(v[num_col-1-j]==0){
                    //nothing
                }
                else{
                    handle_error_with_exit("error in calculate vector base v[num_col-1-j]\n");
                }
            }
            else if(array_var[num_col-1-j]==2){//array_var[num_col-1-j]!=1 || array_var[num_col-1-j]!=3 ->array_var==2==ricavabile
                index_to_calculate=num_col-1-j;//memorizza la colonna in cui c'è la variabile da calcolare,
            }
            else{
                handle_error_with_exit("error in calculate vector base\n");
            }
        }
        if(index_to_calculate<0){//riga in cui non è stato trovato nessun elemento da calcolare
            printf("index_to_calculate=%d\n",index_to_calculate);
            print_array_char(array_var,num_col,thresold);
            handle_error_with_exit("error in index_calculate index<0\n");
        }
        if(index_to_calculate>num_col){
            printf("num_col=%d,index_to_calculate=%d\n",num_col,index_to_calculate);
            handle_error_with_exit("error in index_calculate,index>num_col\n");
        }
        if(num_row_not_null-1-i>num_row || num_row_not_null-1-i<0){
            handle_error_with_exit("error in index calculate vector base\n");
        }
        if((count & 1)==0){//divisibile per 2
            v[index_to_calculate]=0;
        }
        else{
            v[index_to_calculate]=1;
        }
        array_var[index_to_calculate]=3;//imposta variabile come calcolata,piano piano le variabili diventano tutte calcolate
    }
    return;
}

int* alloc_vector_base(char*array_var,int length,int index_vector_base){//lenght=num_row sistema lineare,index parte da 0
	//crea un vettore della base:un 1 su una variabile libera e tutti 0 sulle altre variabili libere rimanenti
	//es array_var = [1 2 1 1 2] -> x1,x3,x4 sono libere
	//vector_base1=[1 x 0 0 y] trasposto
	//vector_base2=[0 x 1 0 y] trasposto
	//vecotr base3=[0 x 0 1 y] trasposto


	//in questo modo i vettori di base sono sicuramente indipendenti a*v1+b*v2+c*v3 è una soluzione del sistema lineare
	//i valore 1 e 0 non sono un caso,infatti quando si esegue l'algoritmo su carta e si trovano i 3 vettori in funzione di v1 v2 e v3 
        //mettendo a fattor comune per il primo vettore il secondo e il terzo rispettivamente v1,v2,v3 si vede che il primo vettore ha un 1 nella posizione v1 e 0 nelle posizioni v2 e v3,il vettore v2 ha un 1 nella posizione v2 e 0 nelle posizioni v1 e v3
	// e il vettore v3 ha un 1 nella posizione v3 e 0 nelle posizioni v1 e v2
	if(array_var==NULL || length<=0 || index_vector_base<0 || index_vector_base>=length){
		handle_error_with_exit("error in parameter alloc_vector_base\n");
	}
	int*vector_base=alloc_array_int(length);
	for(int i=0;i<length;i++){
		if(array_var[i]==3){
			vector_base[i]=0;
			continue;
		}
		else if(array_var[i]==2){//se non è una variabile libera vai avanti,non riempire il valore
			continue;
		}
		else if(array_var[i]!=1){
			handle_error_with_exit("error in alloc_vector_base\n");
		}
		//array_var[i]==1,variabilie iesima è libera
		if(i==index_vector_base){//se il contatore è arrivato all'indice imposta l'unico 1 da impostare per il vettore
			vector_base[i]=1;
		}
		else{//altrimenti metti 0 sulle altre variabili libere
			vector_base[i]=0;
		}
	}
	return vector_base;
}

char* find_free_var(int**matrix_linear_system,int num_row,int num_col){//matrice ridotta a scala,ritorna array con elementi 1 e 2 (ricordiamo che 1=variabile libera 2=variabile ricavabile 3=variabile calcolata)
//es array_var=[1 2 1 1 2 1] ci dice che x1 è libera x2 ricavabile x3 libera x4 libera x5 ricavabile x6 libera
//scansiona tutta la matrice
         
	//si concentra sulla sottomatrice che contiene righe non nulle(solo la sottomatrice non nulla ci da informazioni sulle variabili)
	//idea== se le variabili sono libere hanno già un valore fissato se le variabili sono ricavabili allo step successivo avranno un valore,(es x1+x2=3 partendo da destra x2 è ricavabile x1 è libera)
	//operazioni sulla riga:la prima variabile che si incontra non ancora assegnata diventa ricavabile e le altre non ancora assegnate tutte libere
	//la scansione avviene dall'ultima colonna dell'ultima  riga non nulla fino all'inizio della riga non nulla e poi si passa alla riga superiore sempre partendo dalla fine
	//nota bene:la matrice è ridotta a scala quindi contiene num_row_not_null pivot e quindi deve contenere num_row_not_null variabili 		ricavabili
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_col\n");
	}
	int num_row_not_null=count_rows_not_null(matrix_linear_system,num_row,num_col);
	int rows_null=num_row-num_row_not_null;
	int num_free_var_left=calculate_dim_sol(matrix_linear_system,num_row,num_col);//il numero delle variabili libere è lo stesso della dimensione della base
	int num_var=num_col;//il numero delle variabili è lo stesso del numero delle colonne
	char*array_var=alloc_array_char(num_var);
	for(int i=0;i<num_row-rows_null;i++){//scansiona solo righe non nulle
		char not_null=0;//numero di elementi non nulli incontrati
		for(int j=0;j<num_col;j++){//dopo il primo elemento di riga non nullo gli altri sono tutti variabili libere
			//n.b. num_var=num_col
			if(matrix_linear_system[num_row_not_null-1-i][num_col-1-j]!=0){
				if(array_var[num_var-1-j]==1 || array_var[num_var-1-j]==2){//variabile già assegnata come variabile libera o 						ricavabile
					continue;
				}
				else if(array_var[num_var-1-j]==0 && not_null==0){//se la variabile iesima ancora non è stata assegnata e non sono ancora state incontrate variabili da risolvere o libere 
					//allora assegna la variabile iesima come ricavabile e imposta not_null=1
					array_var[num_var-1-j]=2;//assegnala come ricavabile
					not_null=1;
				}
				else if(array_var[num_var-1-j]==0 && not_null==1){//se la variabile iesima ancora non è stata assegnata ed è stata già incontrata una variabile non nulla(ricavabile o libera)
					//allora assegna la variabile iesima come libera 
					array_var[num_var-1-j]=1;//libera
					num_free_var_left--;//aggiorna il numero di variabili libere ancora da assegnare
				}
				else{
					handle_error_with_exit("error in fin_free_var\n");
				}
			}
		}
	}
	for(int j=0;j<num_var;j++){//gli elementi che sono rimasti a zero diventano variabili già calcolate,necessario se dopo tutto il ciclo alcune variabili non sono state assegnate
		if(array_var[j]==0){
			array_var[j]=3;
		}
	}
	if(num_free_var_left>0){
		handle_error_with_exit("error in find_free_var exit2\n");
	}
	return array_var;
}
char* find_free_var_char(char*matrix_linear_system,int num_row,int num_col,int num_row_not_null,int num_col_not_null){//matrice ridotta a scala,ritorna array con elementi 1 e 2 (ricordiamo che 1=variabile libera 2=variabile ricavabile 3=variabile calcolata)
//es array_var=[1 2 1 1 2 1] ci dice che x1 è libera x2 ricavabile x3 libera x4 libera x5 ricavabile x6 libera
//scansiona tutta la matrice

	//si concentra sulla sottomatrice che contiene righe non nulle(solo la sottomatrice non nulla ci da informazioni sulle variabili)
	//idea== se le variabili sono libere hanno già un valore fissato se le variabili sono ricavabili allo step successivo avranno un valore,(es x1+x2=3 partendo da destra x2 è ricavabile x1 è libera)
	//operazioni sulla riga:la prima variabile che si incontra non ancora assegnata diventa ricavabile e le altre non ancora assegnate tutte libere
	//la scansione avviene dall'ultima colonna dell'ultima  riga non nulla fino all'inizio della riga non nulla e poi si passa alla riga superiore sempre partendo dalla fine
	//nota bene:la matrice è ridotta a scala quindi contiene num_row_not_null pivot e quindi deve contenere num_row_not_null variabili 		ricavabili
	if(matrix_linear_system==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter get_col\n");
	}
	int rows_null=num_row-num_row_not_null;
	int num_free_var_left=calculate_dim_solution(num_row_not_null,num_col_not_null);//il numero delle variabili libere è lo stesso della dimensione della base
	int num_var=num_col;//il numero delle variabili è lo stesso del numero delle colonne
	char*array_var=alloc_array_char(num_var);
	for(int i=0;i<num_row-rows_null;i++){//scansiona solo righe non nulle
		char not_null=0;//numero di elementi non nulli incontrati
		for(int j=0;j<num_col;j++){//dopo il primo elemento di riga non nullo gli altri sono tutti variabili libere
			//n.b. num_var=num_col
			long index=get_index(num_row_not_null-1-i,num_col-1-j,num_col);
			if(matrix_linear_system[index]!=0){
				if(array_var[num_var-1-j]==1 || array_var[num_var-1-j]==2){//variabile già assegnata come variabile libera o 						ricavabile
					continue;
				}
				else if(array_var[num_var-1-j]==0 && not_null==0){//se la variabile iesima ancora non è stata assegnata e non sono ancora state incontrate variabili da risolvere o libere
					//allora assegna la variabile iesima come ricavabile e imposta not_null=1
					array_var[num_var-1-j]=2;//assegnala come ricavabile
					not_null=1;
				}
				else if(array_var[num_var-1-j]==0 && not_null==1){//se la variabile iesima ancora non è stata assegnata ed è stata già incontrata una variabile non nulla(ricavabile o libera)
					//allora assegna la variabile iesima come libera
					array_var[num_var-1-j]=1;//libera
					num_free_var_left--;//aggiorna il numero di variabili libere ancora da assegnare
				}
				else{
					handle_error_with_exit("error in fin_free_var\n");
				}
			}
		}
	}
	for(int j=0;j<num_var;j++){//gli elementi che sono rimasti a zero diventano variabili già calcolate,necessario se dopo tutto il ciclo alcune variabili non sono state assegnate
		if(array_var[j]==0){
			array_var[j]=3;
		}
	}
	if(num_free_var_left>0){
		handle_error_with_exit("error in find_free_var exit2\n");
	}
	return array_var;
}

int** calculate_base_linear_system_char(char*matrix_linear_system,int num_row,int num_col,int*dim_sol,int max_dim_sol){//matrice ridotta modulo n,calcola una base del sistema lineare
    #if DEBUG==1
    if(matrix_linear_system==NULL || num_row<=0 || num_col<=0 || dim_sol==NULL){
        handle_error_with_exit("error in parameter get_coli\n");
    }
    if(check_if_matrix_char_is_reduce_mod_n(matrix_linear_system,num_row,num_col,2)==0){
        handle_error_with_exit("matrix is not reduce mod n");
    }
    if(check_if_matrix_char_is_echelon_reduce(matrix_linear_system,num_row,num_col)==0){
        handle_error_with_exit("error in calculate_base_linear_system\n");
    }
    #endif
    int num_row_not_null,num_col_not_null;
    *dim_sol=calculate_dim_sol_char(matrix_linear_system,num_row,num_col,&num_row_not_null,&num_col_not_null);//calcola la dimensione della base del sistema lineare
    print_time_elapsed("time to calculate_dim_sol");
    if(*dim_sol<0){
        return NULL;
    }
    if(*dim_sol==0){
        return NULL;
    }
    //n.b. int free_var=*dim_sol;
    if(*dim_sol>max_dim_sol){
        *dim_sol=max_dim_sol;
    }
    int**base_linear_system=alloc_matrix_int(num_col,*dim_sol);//la base del sistema lineare ha come righe il numero di colonne della matrice 		(numero delle variabili) e come colonne il numero di vettori linearmente indipendenti
    char*array_var=NULL;
    array_var=find_free_var_char(matrix_linear_system,num_row,num_col,num_row_not_null,num_col_not_null);//calcola array delle variabili libere 
    //che specifica tutte le variabili che sono state impostate come libere per tutto il sistema,è lungo num_col,è necessario 
    //allocarlo ogni volta perchè viene sporcato dalla funzione calculate_vector_base
    
    char*array_var_temp=alloc_array_char(num_col);
    int s=0;//s==start,indice della posizione delle variabili libere
    for(int i=0;i<*dim_sol;i++){//ripeti il procedimento per il numero di vettori linearmente indipendenti,e crea ad ogni ciclo un vettore
        //di base
        memcpy(array_var_temp,array_var,sizeof(char)*num_col);
        //s non va resettato a zero altrimenti troverà sempre l'indice della prima variabile libera
        while(array_var[s]!=1 && s<num_col){//finquando non si arriva ad una variabile libera vai avanti
            s++;
        }
        if(array_var[s]==1){//per ogni variabile libera va calcolato un vettore di base
            int*v=alloc_vector_base(array_var_temp,num_col,s);//alloca e inizializza vettore di base
            //risolvi il sistema lineare per sostituzione sapendo che ogni vettore di base fornisce una soluzione parziale
            calculate_vector_base_char(matrix_linear_system,num_row,num_col,array_var_temp,v,num_row_not_null,max_dim_sol);//memorizza in v la soluzione
            //errore in calculate vector base
            for(int j=0;j<num_col;j++){//copia l'array v nella colonna della matrice base_sistema_lineare
                base_linear_system[j][i]=v[j];//ok
                reduce_int_mod_n(&(base_linear_system[j][i]),2);//ridurre elemento modificato mod 2
            }
            free(v);
            v=NULL;
            s++;
        }
        else{//leggi handle_error
            handle_error_with_exit("la lunghezza dell'array_var è stata superata e non sono stati trovati dim_sol vettori di base\n");
        }
    }
    free(array_var);
    array_var=NULL;
    free(array_var_temp);
    array_var_temp=NULL;
    return base_linear_system;
}
int scan_array_to_find_element_not_null(int*array,int start,int lenght_array){//start=punto di partenza,inizia da 0 lenght=lunghezza vettore,
//scansiona gli elementi del vettore e ritorna l'indice del primo diverso da zero se lo trova,altrimenti -1
	//ritorna la prima occorrenza o -1 se nessuna occorrenza
	if(array==NULL || start<0 || lenght_array<=0 || start>=lenght_array){//start è un indice non può essere uguale a length_array
		handle_error_with_exit("error in parameter scan_array_to_find_element_not_null\n");
	}
	for(int i=start;i<lenght_array;i++){//start è minore di length_array,se l'array non è vuoto fa almeno un ciclo
		if(array[i]!=0){
			return i;//occorrenza trovata
		}
	}
	return -1;//occorrenza non trovata
}



void reduce_echelon_form(int**matrix,int num_row,int num_col){//versione rref,fare il test per verificare
		// che è effettivamente ridotta in modo rref
	int lead=0;
	int i;
	#if DEBUG==1
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	if(check_if_matrix_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
		handle_error_with_exit("matrix is not reduce mod n");
	}
	#endif
	for(int r=0;r<num_row;r++){
		if(num_col<=lead){
			return;
		}
		i=r;
		while(matrix[i][lead]==0){
			i=i+1;
			if(num_row==i){
				i=r;
				lead=lead+1;
				if(num_col==lead){
					return;
				}
			}
		}
		if(i!=r){
			swap_row(matrix,num_row,num_col,i,r);
		}
		#if TODO==1//this check is correct?program break in this point with debug enabled
		if(matrix[r][lead]==0){
			handle_error_with_exit("invalid pivot\n");
		}
		#endif
		for(int i=0;i<num_row;i++){//riduzione della matrice versione rref
			if(i!=r){
				if(matrix[i][lead]!=0){
					for(int f=lead;f<num_col;f++){//in realtà il ciclo si può fare da lead a num_col,la parte prima di lead 						sono sottrazioni per 0
						matrix[i][f]=matrix[i][f]-matrix[r][f];//si fa la differenza perchè il pivot è sempre 1
						matrix[i][f]=reduce_mod_2(matrix[i][f]);
					}
				}
			}			
		}
		lead=lead+1;
	}
	return;
}
void reduce_echelon_form_binary_matrix(unsigned long**binary_matrix,int num_row,int num_col){
	if(binary_matrix==NULL || *binary_matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter reduce echelon form binary matrix\n");
	}
	int remainder;
	int index_element;
	unsigned long bit;
	int lead=0;//colonna
	int i;//indice di riga dove è il pivot
	for(int r=0;r<num_row;r++){//r riga
		if(num_col*BIT_OF_UNSIGNED_LONG<=(unsigned long)lead){
			return;
		}
		i=r;
		remainder=lead%BIT_OF_UNSIGNED_LONG;
		index_element=lead/BIT_OF_UNSIGNED_LONG;
		bit = (binary_matrix[i][index_element] >> (BIT_OF_UNSIGNED_LONG-remainder-1)) & 1U;
		while(bit==0){
			i=i+1;
			if(num_row==i){//se sei arrivato alla fine della matrice ricomincia ma vai alla colonna successiva
				i=r;//rimettiti alla riga r-esima
				lead=lead+1;//colonna successiva
				if(num_col*BIT_OF_UNSIGNED_LONG==(unsigned long)lead){
					return;
				}
			}
			remainder=lead%BIT_OF_UNSIGNED_LONG;
			index_element=lead/BIT_OF_UNSIGNED_LONG;
			bit = (binary_matrix[i][index_element] >> (BIT_OF_UNSIGNED_LONG-remainder-1)) & 1U;
		}
		if(i!=r){//se gli indici di riga sono diversi swappa le righe
			swap_row_unsigned_long(binary_matrix,num_row,i,r);
		}
		#if TODO==1
		#if DEBUG==1 //this check is valid?program break in this point
		if(binary_matrix[r][lead]==0){
			handle_error_with_exit("invalid pivot\n");
		}
		#endif
		#endif

		//in r,lead c'è il pivot sottrai tutta la sottomatrice
		for(int i=r+1;i<num_row;i++){//versione ref,sottrai sottomatrice sotto pivot
		//for(int i=0;i<num_row;i++){//riduzione della matrice versione rref
			if(i!=r){
				remainder=lead%BIT_OF_UNSIGNED_LONG;
				index_element=lead/BIT_OF_UNSIGNED_LONG;
				bit = (binary_matrix[i][index_element] >> (BIT_OF_UNSIGNED_LONG-remainder-1)) & 1U;
				if(bit!=0){
					for(int j1=0;j1<num_col;j1++){//xor tra le 2 rihe
						binary_matrix[i][j1]=binary_matrix[i][j1]^binary_matrix[r][j1];
					}
				}
			}
		}
		lead=lead+1;//vai alla riga successiva
	}
	return;
}
void reduce_echelon_form_matrix_char(char**matrix,int num_row,int num_col){//versione rref,fare il test per verificare
	// che è effettivamente ridotta in modo rref
	int lead=0;
	int i;
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	/*if(check_if_matrix_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
		handle_error_with_exit("matrix is not reduce mod n");
	}*/
	for(int r=0;r<num_row;r++){
		printf("r=%d,num_row=%d\n",r,num_row);
		if(num_col<=lead){
			return;
		}
		i=r;
		while(matrix[i][lead]==0){
			i=i+1;
			if(num_row==i){
				i=r;
				lead=lead+1;
				if(num_col==lead){
					return;
				}
			}
		}
		if(i!=r){
			swap_row_char2(matrix,num_row,num_col,i,r);
		}
		if(matrix[r][lead]==0){
			handle_error_with_exit("invalid pivot\n");
		}
		//in r,lead c'è il pivot sottrai tutta la sottomatrice
        //for(int i=r+1;i<num_row;i++){//versione ref
		for(int i=0;i<num_row;i++){//riduzione della matrice versione rref
			if(i!=r){
				if(matrix[i][lead]!=0){
					for(int f=lead;f<num_col;f++){//in realtà il ciclo si può fare da lead a num_col,la parte prima di lead 						sono sottrazioni per 0
						//matrix[i][f]=matrix[i][f]-matrix[r][f];//si fa la differenza perchè il pivot è sempre 1,nessun moltiplicatore
						//matrix[i][f]=(char)reduce_mod_2((char)matrix[i][f]);
						if(matrix[i][f]!=matrix[r][f]){
							matrix[i][f]=1;
						}
						else{
							matrix[i][f]=0;
						}
					}
				}
			}
		}
		lead=lead+1;
	}
	return;
}
void reduce_echelon_form_char(char*matrix,int num_row,int num_col){//versione rref,fare il test per verificare
    // che è effettivamente ridotta in modo rref
    int lead=0;
    int i;
    long index,index1;
    #if DEBUG==1
    if(matrix==NULL || num_row<=0 || num_col<=0){
        handle_error_with_exit("error in parameter reduce echelon form\n");
    }
    if(check_if_matrix_char_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
        handle_error_with_exit("matrix is not reduce mod n");
    }
    #endif
    printf("inizio echelon form\n");
    for(int r=0;r<num_row;r++){
        printf("r=%d,num_row=%d\n",r,num_row);
        if(num_col<=lead){
            return;
        }
        i=r;
        index=get_index(i,lead,num_col);//moltiplicazione
        while(matrix[index]==0){
            i=i+1;
            if(num_row==i){
                i=r;
                lead=lead+1;
                if(num_col==lead){
                    return;
                }
            }
            index=get_index(i,lead,num_col);//moltiplicazione
        }
        if(i!=r){
            swap_row_char(matrix,num_row,num_col,i,r);
        }
        index=get_index(r,lead,num_col);
        if(matrix[index]==0){
            handle_error_with_exit("invalid pivot\n");
        }
        for(int i=r+1;i<num_row;i++){//versione ref
        //for(int i=0;i<num_row;i++){//riduzione della matrice versione rref
            if(i!=r){
                index=get_index(i,lead,num_col);
                if(matrix[index]!=0){
                    for(int f=lead;f<num_col;f++){//in realtà il ciclo si può fare da lead a num_col,la parte prima di lead	sono sottrazioni per 0
                        index=get_index(i,f,num_col);
                        index1=get_index(r,f,num_col);
                        matrix[index]=matrix[index]-matrix[index1];//si fa la differenza perchè il pivot è sempre 1
                        matrix[index]=reduce_mod_2((int)matrix[index]);
                    }
                }
            }
        }
        lead=lead+1;
    }
    return;
}



long get_index(int index_row,int index_col,int num_col){
	//ritorna l'indice della posizione numerica associata all'elemento index_row,index_col
	long index;
	if(index_col<0 || index_row<0 || num_col<=0){
		handle_error_with_exit("error in get_index\n");
	}
	index=(long)index_row*num_col;//elemento della riga
	index+=index_col;//shift della colonna
	return index;
}
char array_is_fill_of_value(char*combination,int length,char value){//verifica che ogni elemento dell'array è uguale a value
	if(length<=0 || combination==NULL){
		handle_error_with_exit("error in parameter\n");
	}
	for(int i=0;i<length;i++){
		if(combination[i]!=value){
			return 0;
		}
	}
	return 1;
}

void free_memory_matrix_long(long **matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix int\n");
	}
	for(int i=0;i<num_row;i++){
		free(matrix[i]);
	}
	free(matrix);
	matrix=NULL;
	return;
}
void free_memory_matrix_int(int **matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix int\n");
	}
	for(int i=0;i<num_row;i++){
		free(matrix[i]);
	}
	free(matrix);
	matrix=NULL;
	return;
}
void free_memory_matrix_unsigned_long(unsigned long **matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix int\n");
	}
	for(int i=0;i<num_row;i++){
		free(matrix[i]);
	}
	free(matrix);
	matrix=NULL;
	return;
}
void free_memory_matrix_char(char **matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix int\n");
	}
	for(int i=0;i<num_row;i++){
		free(matrix[i]);
	}
	free(matrix);
	matrix=NULL;
	return;
}









