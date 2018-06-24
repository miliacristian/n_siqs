#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix_function.h"
#include "basic.h"
#include "list_factor_base.h"
#include "list_square_relation.h"
#include "list_factorization.h"
#include <mpfr.h>
extern struct row_factorization r;
extern struct timespec timer;
extern struct timespec time_start;
extern FILE*file_log;
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
char value_is_in_sorted_array(int index_of_prime,struct a_struct*array_a_struct,int length){
    if(array_a_struct==NULL || length<0){
        handle_error_with_exit("error in value is in sorted_array\n");
    }
    //se è minore del primo o maggiore dell'ultimo allora non è nell'array
    if(index_of_prime<array_a_struct[0].index_prime_a || index_of_prime>array_a_struct[length-1].index_prime_a){
        return 0;
    }
    //se sono uguali ritorna 1 altrimenti 0
    for(int i=0;i<length;i++){
        if(index_of_prime==array_a_struct[i].index_prime_a){
            return 1;
        }
    }
    return 0;
}

/*void free_memory_matrix_factorization(struct matrix_factorization *mat){
	if(mat==NULL){
		handle_error_with_exit("error in free memory matrix factorization\n");
	}
	for(int i=0;i<mat->num_row;i++){
		mpz_clear((mat->row[i]).num);
		mpz_clear((mat->row[i]).square);
	}
	free(mat->row);
	mat->row=NULL;
	free(mat);
	mat=NULL;
	return;
}

struct row*alloc_array_row(int num_row){
	if(num_row<0){
		handle_error_with_exit("error in row_factorization\n");
	}
	struct row* row;
	row=malloc(sizeof(struct row)*num_row);
	if(row==NULL){
		handle_error_with_exit("error in row_factorization\n");
	}
	memset(row,0,sizeof(struct row)*num_row);
	for(int i=0;i<num_row;i++){
		row[i].index_last_prime=-1;
		row[i].index_first_prime=-1;
		mpz_init(row[i].num);
		mpz_init(row[i].square);
		row[i].sum_log=0;
		row[i].log=-1;
	}
	return row;
}

struct matrix_factorization* alloc_matrix_factorization(int num_row){
	if(num_row<0){
		handle_error_with_exit("error in matrix_factorization alloc_matrix_factorization\n");
	}
	struct matrix_factorization*m=malloc(sizeof(struct matrix_factorization)*1);//ritorna il puntatore all'area di memoria dove risiede
	//la matrice di fattorizzazione
	if(m==NULL){
		handle_error_with_exit("error in matrix_factorization alloc_matrix_factorization\n");
	}
	m->num_row=num_row;
	m->row=alloc_array_row(num_row);
	return m;
}

struct matrix_factorization**alloc_array_matrix_factorization(int length_array_matrix){
	if(length_array_matrix<=0){
		handle_error_with_exit("error in alloc array_matrix_factorization\n");
	}
	struct matrix_factorization**array_matrix_factorization;
	array_matrix_factorization=malloc(sizeof(struct matrix_factorization*)*length_array_matrix);
	if(array_matrix_factorization==NULL){
		handle_error_with_exit("error in malloc array_matrix_factorization\n");
	}
	return array_matrix_factorization;
}
*/

//scriverla in modo parallelizzato
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s){
	if(head_f_base_f==NULL || card_f_base<=0 || (array_a_struct==NULL && s>0) || s<0){
		handle_error_with_exit("error in create row factorization\n");
	}
	mpz_t temp;
	int index=0;
	mpz_init(temp);
	struct node_factor_base*p=head_f_base_f;
	for(int i=0;i<card_f_base;i++){
        r.prime[i]=p->prime;//metti il primo della factor base in posizione pari
		r.root_n_mod_p[i]=p->root_n_mod_prime;
		r.root2_n_mod_p[i]= r.prime[i]-r.root_n_mod_p[i];//rad2=p-rad1 mod p
		if((i==0 && s>0) || (i==1 && s>0) ||
		   (index<s && i==array_a_struct[index].index_prime_a) ) {//p divide a non esiste inverso modulo p
			r.inverse_a_mod_p[i]=-1;
            if(i!=0 && i!=1) {
				index++;
			}
        }
        else if(s==0){//a=1
			r.inverse_a_mod_p[i]=1;
		}
        else{//p non divide a ->esiste l'inverso
		    mpz_set_si(temp,r.prime[i]);//temp=p
            if(mpz_invert(temp,a,temp)==0){//temp=inverse of a mod p
				if(i==array_a_struct[index].index_prime_a){
					handle_error_with_exit("exit\n");
				}
            	handle_error_with_exit("error in mpz_invert create_row_factorization\n");
            }
			r.inverse_a_mod_p[i]=mpz_get_si(temp);
		}
		if(p->prime!=-1){
			r.log_prime[i]=(int)round(log2f((float)p->prime));
		}
		p=p->next;
	}
	if(index!=s){
	    handle_error_with_exit("error in initialize inverse mod p create_row_factorization\n");
	}
	mpz_clear(temp);
	return;
}
/*struct matrix_factorization*create_matrix_factorization_f(int M,int card_f_base,const mpz_t a,const mpz_t b,const mpz_t n){
	if(M<=0 || card_f_base<=0 || a==NULL || b==NULL || n==NULL){
		handle_error_with_exit("error in create_matrix_factorization_f\n");
	}
	struct matrix_factorization *m;
	m=alloc_matrix_factorization(2*M+1);
	char test=((mpz_get_si(a)!=0) && (mpz_get_si(a)!=1));//test==1 se a è contemporaneamente diverso da 1 e da 0
	long j=-M;
	long h=1+2*j;
	int index;
	mpz_t c;
	mpz_t v,v1;
	mpz_t double_a,double_b,t_j,a_temp;
	mpfr_t log_value;

	mpz_init(c);
	mpz_init(v);
	mpz_init(v1);
	mpz_init(double_a);
	mpz_init(double_b);
	mpz_init(t_j);
	mpz_init(a_temp);
	mpfr_init(log_value);

	//calcolo di c
	mpz_mul(v,b,b);//v=b^2
	mpz_sub(v,v,n);//v=b^2-n
	if(mpz_divisible_p(v,a)==0){
		handle_error_with_exit("error in mpz_divisible create array_of_number\n");
	}
	mpz_divexact(c,v,a);//c=(b^2-n)/a
	//calcolo di double_a
	mpz_set(double_a,a);//double_a=a
	mpz_mul_ui(double_a,double_a,2);//double_a=2*a
	//calcolo di double_b
	mpz_set(double_b,b);//double_b=b
	mpz_mul_ui(double_b,double_b,2);//double_b=2*b
	//calcolo di t_j=a*(1+2*j)
	mpz_mul_si(t_j,a,h);//t_j=a*(1+2*j)

	//quadrato da aggiungere nella matrice e numero
	//imposta primo valore dell'array
	mpz_mul_si(v,a,j*j);//v=a*j^2
	mpz_set_si(v1,2*j);//v1=2*j
	mpz_mul(v1,v1,b);//v1=2*j*b
	mpz_add(v,v,v1);//v=a*j^2+2*j*b
	mpz_add(v,v,c);//v=a*j^2+2*j*b+c
	//array[i]=a*j^2+2*j*b+c;//a(j)=a*(j^2)+2*bj+c se a=1 b=x0 a(j)=(x0+j)^2-n
	mpz_set((*m).row[0].num,v);//array[i]=a*j^2+2*j*b+c;
	mpz_mul_si(a_temp,a,j);//a_temp=a*j
	mpz_add(a_temp,a_temp,b);//a_temp=a*j+b
	mpz_set((*m).row[0].square,a_temp);//imposta alla prima riga (aj+b) dopo di che ai prossimi elementi aggiungi iterativamente a
	mpfr_set_z(log_value,(*m).row[0].num,MPFR_RNDN);
	if(mpz_cmp_si((*m).row[0].num,0)<0){//se è negativo
		mpfr_neg(log_value,log_value,MPFR_RNDN);
	}
	mpfr_log(log_value,log_value,MPFR_RNDN);
	(*m).row[0].log=mpfr_get_d(log_value,MPFR_RNDN);
	//imposta valori successivi

	for(int i=1;i<2*M+1;i++){
		mpz_add(v,v,t_j);//v=v+t_j
		mpz_add(v,v,double_b);//v=v+t_j+2b
		mpz_set((*m).row[i].num,v);//imposta il valore dell'array
		mpfr_set_z(log_value,(*m).row[i].num,MPFR_RNDN);
		if(mpz_cmp_si((*m).row[i].num,0)<0){//se è negativo
			mpfr_neg(log_value,log_value,MPFR_RNDN);
		}
		mpfr_log(log_value,log_value,MPFR_RNDN);
		(*m).row[i].log=mpfr_get_d(log_value,MPFR_RNDN);
		mpz_add(t_j,t_j,double_a);//calcola il nuovo t_j
		mpz_add((*m).row[i].square,(*m).row[i-1].square,a);
	}
	mpz_clear(c);
	mpz_clear(v);
	mpz_clear(v1);
	mpz_clear(double_a);
	mpz_clear(double_b);
	mpz_clear(t_j);
	mpz_clear(a_temp);
	mpfr_clear(log_value);
	mpfr_free_cache();
	return m;
}*/

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
void copy_matrix(mpz_t**matrix1,mpz_t**matrix2,int num_row,int num_col){
	if(matrix1==NULL || *matrix1==NULL || matrix2==NULL || *matrix2==NULL
	|| num_row<=0 || num_col<=0){
		handle_error_with_exit("error in copy_array_mpz\n");
	}
	for(int i=0;i<num_row;i++){
		copy_array_mpz(matrix1[i],matrix2[i],num_col);//copia riga
	}
	return;
}
char*alloc_array_char(long length){
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

long*alloc_array_long(int length){
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
unsigned long*alloc_array_unsigned_long(int length){
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
int*alloc_array_int(int length){
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
double*alloc_array_double(int length){
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
float*alloc_array_float(int length){
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
long**alloc_array_pointer_to_long(int length){
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
unsigned long**alloc_array_pointer_to_unsigned_long(int length){
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
int**alloc_array_pointer_to_int(int length){
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

/*struct row_factorization*realloc_array_row_factorization(struct row_factorization*row,int new_length){
	if(new_length<=0 || row==NULL){
		handle_error_with_exit("error in parameter realloc array_row_factorization\n");
	}
	struct row_factorization*r=realloc(row,sizeof(struct row_factorization)*new_length);
	if(r==NULL){
		handle_error_with_exit("error in realloc array row factorization\n");
	}
	return r;
}
struct row*realloc_array_row(struct row*row,int new_length){
	if(new_length<0 || row==NULL){
		handle_error_with_exit("error in parameter realloc_array_row\n");
	}
	struct row*r=realloc(row,sizeof(struct row)*new_length);
	if(r==NULL){
		handle_error_with_exit("error in realloc array row\n");
	}
	return r;
}*/
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

double**alloc_array_pointer_to_double(int length){
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
char**alloc_array_pointer_to_char(int length){
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
double**alloc_matrix_double(int num_row,int num_col){
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
char**alloc_matrix_char(int num_row,int num_col){
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
int**alloc_matrix_int(int num_row,int num_col){
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
long**alloc_matrix_long(int num_row,int num_col){
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
unsigned long**alloc_matrix_unsigned_long(int num_row,int num_col){
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
/*void concatenate_matrix_long(long***res,long**matrix1,int num_row1,int num_col1,long**matrix2,int num_row2,int num_col2){//concatena 2 matrici
//allungando le righe della matrice che contiene più righe
	if(matrix1==NULL || *matrix1==NULL || matrix2==NULL || *matrix2==NULL || num_row1<=0 
		|| num_col1<=0 || num_row2<=0 || num_col2<=0
		|| num_col1!=num_col2){
		handle_error_with_exit("invalid parameter concatenate_matrix_mpz\n");
	}
	*res=alloc_array_pointer_to_long(num_row1+num_row2);
	for(int i=0;i<num_row1;i++){
		(*res)[i]=matrix1[i];
	}
	print_matrix_long(*res,num_row1,num_col1);
	for(int j=num_row1;j<num_row1+num_row2;j++){
		(*res)[j]=matrix2[j-num_row1];
	}
	return;
}


void concatenate_matrix_mpz(mpz_t***pmatrix1,int *num_row1,int num_col1,mpz_t**matrix2,int num_row2,int num_col2){//concatena 2 matrici
//allungando le righe della prima matrice
	if(pmatrix1==NULL || *pmatrix1==NULL || **pmatrix1==NULL || matrix2==NULL || *matrix2==NULL || num_row1==NULL || *num_row1<=0 
		 || num_col1<=0  || num_row2<=0 || num_col2<=0
		|| num_col1!=num_col2){
		handle_error_with_exit("invalid parameter concatenate_matrix_mpz\n");
	}
	*pmatrix1=realloc_array_pointer_to_mpz(*pmatrix1,*num_row1+num_row2);//rialloca prima matrice

	for(int j=*num_row1;j<*num_row1+num_row2;j++){
		(*pmatrix1)[j]=matrix2[j-*num_row1];
	}
	*num_row1=*num_row1+num_row2;
	free(matrix2);//libera la memoria solamente dell'array pointer to mpz,quello che contiene i puntatori a riga della matrice
	matrix2=NULL;
	return;
}
void concatenate_matrix_factorization(struct matrix_factorization*pmatrix1,int *num_row1,struct matrix_factorization *matrix2,int num_row2){//concatena 2 matrici
//allungando le righe della prima matrice
	if(pmatrix1==NULL || num_row1==NULL || *num_row1<0  || num_row2<0 || matrix2==NULL || matrix2->num_row<0){
		handle_error_with_exit("invalid parameter concatenate_matrix_factorization\n");
	}
	struct matrix_factorization {//matrice che contiene tutte le fattorizazzioni
	struct row*row;//2*M+1
	int length;
	};
	if(num_row2==0){//libera memoria
		free(matrix2->row);//libera la memoria solamente dell'array pointer to mpz,quello che contiene i puntatori a riga della matrice
		matrix2->row=NULL;
		free(matrix2);
		matrix2=NULL;
		return;
	}
	pmatrix1->row=realloc_array_row(pmatrix1->row,*num_row1+num_row2);//cambia dimensione puntatore
	pmatrix1->num_row=*num_row1+num_row2;//cambia lunghezza matrice

	for(int j=*num_row1;j<*num_row1+num_row2;j++){//copia delle righe
		(*pmatrix1).row[j]=(*matrix2).row[j-*num_row1];
	}
	*num_row1=*num_row1+num_row2;
	free(matrix2->row);//libera la memoria solamente dell'array pointer to mpz,quello che contiene i puntatori a riga della matrice
	matrix2->row=NULL;
	free(matrix2);
	matrix2=NULL;
	return;
}
void concatenate_all_matrix_mpz_same_dimension(mpz_t***array_matrix_mpz,int length_array_matrix,int *row_result,int num_row_matrix,int num_col_matrix){
	if(array_matrix_mpz==NULL || *array_matrix_mpz==NULL || **array_matrix_mpz==NULL || length_array_matrix<=0 || row_result==NULL || 			num_row_matrix<=0 || num_col_matrix<=0){
		handle_error_with_exit("error in concatenate_all_matrix_mpz_same_dimension\n");
	}
	if(length_array_matrix==1){
		*row_result=num_row_matrix;
		return;
	}
	
	//lenght_array_matrix>1
	*row_result=num_row_matrix;//inizializza riga dei risultati
	for(int i=0;i<length_array_matrix-1;i++){
		concatenate_matrix_mpz(array_matrix_mpz,row_result,num_col_matrix,array_matrix_mpz[i+1],num_row_matrix,num_col_matrix);
	}
	return;
}
void concatenate_all_matrix_B_smooth(struct matrix_factorization**array_matrix_factorization,int length_array_matrix,int *row_result){
	//le matrici non possono essere anche ==NULL
	if(array_matrix_factorization==NULL || *array_matrix_factorization==NULL || length_array_matrix<=0 || row_result==NULL){
		handle_error_with_exit("error in concatenate_all_matrix_mpz_same_dimension\n");
	}
	if(length_array_matrix==1){//prendi solamente la prima matrice B_smooth
		*row_result=(*array_matrix_factorization)[0].num_row;
		return;
	}
	//lenght_array_matrix>1
	*row_result=(*(array_matrix_factorization[0])).num_row;//inizializza riga dei risultati
	for(int i=0;i<length_array_matrix-1;i++){
		concatenate_matrix_factorization(*array_matrix_factorization,row_result,array_matrix_factorization[i+1],
		(*(array_matrix_factorization[i+1])).num_row);
	}
	return;
}*/

/*void adjust_array_of_prime(long*array_of_prime_chosen_for_a,int length_array,long*best_q,int length_best_q){//mette degli 1 in posizione opportuna
	long index=-1;
	if(array_of_prime_chosen_for_a==NULL || length_array<=0 || best_q==NULL || length_best_q<=0){
		handle_error_with_exit("error in adjust_array_of_prime\n");
	}
	for(int i=0;i<length_best_q;i++){
		index=best_q[i];
		if(index<2){
			handle_error_with_exit("error in adjust_array_of_prime,primes chosen not good\n");
		}
		if(2*index+1>=length_array){
			handle_error_with_exit("error in index into array best_q\n");
		}
		array_of_prime_chosen_for_a[2*index+1]=1;
	}
	return;
}*/

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
/*void find_relation_one_partial(mpz_t**matrix_factorization,long size,int cardinality_factor_base,const mpz_t n,int k,char*factorized){
	if(matrix_factorization==NULL || *matrix_factorization==NULL || size<=0 || cardinality_factor_base<=0 || n==NULL || factorized==NULL){
		handle_error_with_exit("error in find_relation_one_partial\n");
	}
	mpz_t inverse,temp;
	long len_only_prime=cardinality_factor_base*2;//numero di colonne dove risiedono solo i primi della factor base
	int index_last_element=cardinality_factor_base*2+1;
	//matrix_factorization[i][index_last_element-1]==quadrato
	//matrix_factorization[i][index_last_element]==residuo della divisione
	mpz_init(inverse);
	mpz_init(temp);
	for(long i=0;i<size-1;i++){
		if(mpz_cmp_si(matrix_factorization[i][index_last_element],1)==0){//se il residuo è 1 continua,nessuna relazione parziale
			continue;
		}
		//residuo diverso da 1
		if(mpz_cmp(matrix_factorization[i][index_last_element],matrix_factorization[i+1][index_last_element])==0){
			//se i residui presenti in 2 righe diverse sono uguali->possibile relazione 1 parziale trovata
			if(mpz_cmp(matrix_factorization[i][index_last_element-1],matrix_factorization[i+1][index_last_element-1])==0){
				//se i quadrati sono uguali(è ovvio che i residui sono uguali) e non va bene come relazione 1 parziale
				continue;//relazione 1 parziale trovata ma doppiona(i 2 quarati sono uguali->anche i residui sono uguali)
			}
			else{//possibile relazione 1 parziale trovata
				mpz_gcd(inverse,matrix_factorization[i][index_last_element],n);//gcd(residuo,n)
				if(mpz_cmp_si(inverse,1)!=0 && mpz_cmp(inverse,n)!=0){//gcd !=0 e gcd!=n
					if(mpz_cmp_si(inverse,k)!=0){//gcd !=0 e gcd!=n gcd!k
						mpz_divexact(temp,n,inverse);//temp=n/inverse
						if(mpz_divisible_ui_p(inverse,k)!=0){//inverse è divisibile per k
							mpz_divexact_ui(inverse,inverse,k);
						}
						else{
							mpz_divexact_ui(temp,temp,k);
						}
						gmp_printf("factorization of n=%Zd*%Zd\n",inverse,temp);//n=(inverse)*(n/inverse)
						//tempo totale,imposta il tempo iniziale alla struct,tempo totale=get_time-tempo iniziale
						fprintf(file_log,"p=");
						mpz_out_str(file_log,10,inverse);
						fprintf(file_log," ");
						fprintf(file_log,"q=");
						mpz_out_str(file_log,10,temp);
						fprintf(file_log," ");
						*factorized=1;
						mpz_clear(inverse);
						mpz_clear(temp);
						return;
					}
					else{//gcd !=0 e gcd!=n gcd=k,fattore banale di n,(n ha il fattore k)
						continue;//non si può trovare la relazione 1-parziale
					}
				}
				//altrimenti esiste l'inverso e calcola relazione 1 parziale
				mpz_mul(matrix_factorization[i][index_last_element-1],matrix_factorization[i][index_last_element-1],
				matrix_factorization[i+1][index_last_element-1]);//moltiplicazione dei quadrati nelle 2 righe
				mpz_mod(matrix_factorization[i][index_last_element-1],matrix_factorization[i][index_last_element-1],n);//riduzione 					modulo n del quadrato
				if(mpz_invert(inverse,matrix_factorization[i][index_last_element],n)==0){//calcolo inverso del residuo mod n
					handle_error_with_exit("error in mpz_invert\n");
				}
				mpz_mul(matrix_factorization[i][index_last_element-1],
				matrix_factorization[i][index_last_element-1],inverse);//moltiplicazione quadrati e inverso del residuo
				mpz_mod(matrix_factorization[i][index_last_element-1],matrix_factorization[i][index_last_element-1],n);
				sum_elem_multiple_of_2_mpz(matrix_factorization[i],matrix_factorization[i+1],len_only_prime,len_only_prime);
				mpz_set_si(matrix_factorization[i][index_last_element],1);//imposta il residuo a 1
			}
		}
	}
	mpz_clear(inverse);
	mpz_clear(temp);
	return;
}*/

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
/*void sort_square(mpz_t**matrix_factorization,long size,int cardinality_factor_base,long M){//ordina i residui in modo crescente
	//ordina usando radixsort(solo per numeri interi>=0),per numeri interi <0 bisogna fare delle piccole modifiche
	//size==righe,cardinality_factor_base*2+2=num_col
	if(size<=0 || matrix_factorization==NULL || *matrix_factorization==NULL || cardinality_factor_base<=0 || M<=0){
		handle_error_with_exit("error in sort_array_number\n");
	}
	long *bucket=NULL;
	int i;
	mpz_t max,temp,significant_digit,q,base;
	mpz_t*semi_sorted=alloc_array_mpz(size);
	mpz_t*array=alloc_array_mpz(size);
	copy_column_matrix_mpz(array,matrix_factorization,size,cardinality_factor_base*2);//copia numeri che devono essere ordinati
	int num_row=size;
	int num_col=cardinality_factor_base*2+2;
	mpz_t**matrix_semi_sorted=alloc_matrix_mpz(num_row,num_col);
	
	mpz_init(max);
	mpz_init(temp);
	mpz_init(q);
	mpz_init(base);
	mpz_init(significant_digit);

	mpz_set_si(base,M);//rimuoverlo
	find_max_array_mpz(max,array,size);
	bucket=alloc_array_long(mpz_get_si(base));
	mpz_set_si(significant_digit,1);//significant_digit=1
	mpz_fdiv_q(q,max,significant_digit);//q*significant_digit+r=max
	while(mpz_cmp_si(q,0)>0){//finquando non si arriva alla prima cifra più significativa del massimo nell'array
		memset(bucket,0,sizeof(long)*mpz_get_si(base));//azzera array bucket
    	// Counts the number of "keys" or digits that will go into each bucket
    		for (i=0;i<size; i++){//conta occorrenze dei numeri che vanno in un certo bucket es a=[123 453 768] se stiamo scansionando l'ultima 		cifra restituisce bucket b=[0 0 0 2 0 0 0 0 1 0]
			mpz_fdiv_q(temp,array[i],significant_digit);//parte intera inferiore(temp)=array[i]/significant_digit
			mpz_mod(temp,temp,base);
			bucket[mpz_get_si(temp)]++;
		}
    		
     		 //Add the count of the previous buckets,
     		// Acquires the indexes after the end of each bucket location in the array
			//  Works similar to the count sort algorithm
		     
    		for (i = 1; i <mpz_get_si(base); i++){
      			bucket[i] += bucket[i - 1];
		}
    		// Use the bucket to fill a "semiSorted" array
    		for (i = size - 1; i >= 0; i--){
			mpz_fdiv_q(temp,array[i],significant_digit);//parte intera inferiore(temp)=array[i]/significant_digit 
			mpz_mod(temp,temp,base);
			bucket[mpz_get_si(temp)]--;
      			mpz_set(semi_sorted[bucket[mpz_get_si(temp)]],array[i]);
			copy_array_mpz(matrix_semi_sorted[bucket[mpz_get_si(temp)]],matrix_factorization[i],num_col);
		}
    		for (i = 0; i < size; i++){
      			mpz_set(array[i],semi_sorted[i]);//array[i]=semisorted[i]
			copy_array_mpz(matrix_factorization[i],matrix_semi_sorted[i],num_col);
		}
    		// Move to next significant digit
    		mpz_mul(significant_digit,significant_digit,base);//significant_digit=significant_digit*base,passa alla prossima cifra(più significativa rispetto a quella attuale)
		mpz_fdiv_q(q,max,significant_digit);//q*significant_digit+r=max
  	}
	copy_matrix(matrix_factorization,matrix_semi_sorted,num_row,num_col);
	free_memory_matrix_mpz(matrix_semi_sorted,num_row,num_col);
	free_memory_array_mpz(semi_sorted,size);
	free_memory_array_mpz(array,size);
	free(bucket);
	mpz_clear(q);
	mpz_clear(temp);
	mpz_clear(max);
	mpz_clear(base);
	mpz_clear(significant_digit);
	return;
}*/

/*void sort_residuos(long size,mpz_t**matrix_factorization,int cardinality_factor_base,long M){//ordina i residui in modo crescente
	//ordina usando radixsort(solo per numeri interi>=0),per numeri interi <0 bisogna fare delle piccole modifiche
	//size==righe,cardinality_factor_base*2+2=num_col
	if(size<=0 || matrix_factorization==NULL || *matrix_factorization==NULL || cardinality_factor_base<=0 || M<=0){
		handle_error_with_exit("error in sort_array_number\n");
	}
	long *bucket=NULL;
	int i;
	mpz_t max,temp,significant_digit,q,base;
	mpz_t*semi_sorted=alloc_array_mpz(size);
	mpz_t*array=alloc_array_mpz(size);
	copy_column_matrix_mpz(array,matrix_factorization,size,cardinality_factor_base*2+1);//copia numeri che devono essere ordinati
	int num_row=size;
	int num_col=cardinality_factor_base*2+2;
	mpz_t**matrix_semi_sorted=alloc_matrix_mpz(num_row,num_col);
	
	mpz_init(max);
	mpz_init(temp);
	mpz_init(q);
	mpz_init(base);
	mpz_init(significant_digit);

	mpz_set_si(base,M);
	find_max_array_mpz(max,array,size);
	bucket=alloc_array_long(mpz_get_si(base));
	mpz_set_si(significant_digit,1);//significant_digit=1
	mpz_fdiv_q(q,max,significant_digit);//q*significant_digit+r=max
	while(mpz_cmp_si(q,0)>0){//finquando non si arriva alla prima cifra più significativa del massimo nell'array
		memset(bucket,0,sizeof(long)*mpz_get_si(base));//azzera array bucket
    	// Counts the number of "keys" or digits that will go into each bucket
    		for (i=0;i<size; i++){//conta occorrenze dei numeri che vanno in un certo bucket es a=[123 453 768] se stiamo scansionando l'ultima 		cifra restituisce bucket b=[0 0 0 2 0 0 0 0 1 0]
			mpz_fdiv_q(temp,array[i],significant_digit);//parte intera inferiore(temp)=array[i]/significant_digit
			mpz_mod(temp,temp,base);
			bucket[mpz_get_si(temp)]++;
		}
    		
     		 //Add the count of the previous buckets,
     		// Acquires the indexes after the end of each bucket location in the array
			//  Works similar to the count sort algorithm
		     
    		for (i = 1; i <mpz_get_si(base); i++){
      			bucket[i] += bucket[i - 1];
		}
    		// Use the bucket to fill a "semiSorted" array
    		for (i = size - 1; i >= 0; i--){
			mpz_fdiv_q(temp,array[i],significant_digit);//parte intera inferiore(temp)=array[i]/significant_digit 
			mpz_mod(temp,temp,base);
			bucket[mpz_get_si(temp)]--;
      			mpz_set(semi_sorted[bucket[mpz_get_si(temp)]],array[i]);
			copy_array_mpz(matrix_semi_sorted[bucket[mpz_get_si(temp)]],matrix_factorization[i],num_col);
		}
    		for (i = 0; i < size; i++){
      			mpz_set(array[i],semi_sorted[i]);//array[i]=semisorted[i]
			copy_array_mpz(matrix_factorization[i],matrix_semi_sorted[i],num_col);
		}
    		// Move to next significant digit
    		mpz_mul(significant_digit,significant_digit,base);//significant_digit=significant_digit*base,passa alla prossima cifra(più significativa rispetto a quella attuale)
		mpz_fdiv_q(q,max,significant_digit);//q*significant_digit+r=max
  	}
	copy_matrix(matrix_factorization,matrix_semi_sorted,num_row,num_col);
	free_memory_matrix_mpz(matrix_semi_sorted,num_row,num_col);
	free_memory_array_mpz(semi_sorted,size);
	free_memory_array_mpz(array,size);
	free(bucket);
	mpz_clear(q);
	mpz_clear(temp);
	mpz_clear(max);
	mpz_clear(base);
	mpz_clear(significant_digit);
	return;
}*/
/*mpz_t**create_matrix_factorization(long M,int cardinality_factor_base,struct node*head_f_base,long*array_of_prime_chosen_for_a,const mpz_t a,const mpz_t b){
	if(M<=0 || cardinality_factor_base<=0 || (array_of_prime_chosen_for_a==NULL && mpz_get_si(a)!=1) || a==NULL || b==NULL){
		handle_error_with_exit("error in parameter create matrix factorization\n");
	}
	char test=((mpz_get_si(a)!=0) && (mpz_get_si(a)!=1));//test==1 se a è contemporaneamente diverso da 1 e da 0
	long length=cardinality_factor_base*2;
	long j=-M;
	mpz_t v,a_temp,j_temp;
	mpz_init(v);
	mpz_init(a_temp);
	mpz_init(j_temp);
	printf("inizio creazione matrix_factorization\n");
	mpz_t**matrix_factor=alloc_matrix_mpz(2*M+1,length+2);//gli ultimi 2 elementi servono per memorizzare il quadrato e il residuo(se residuo è 1 il num è b_smooth)
	for(long j=0;j<cardinality_factor_base;j++){
			get_element_linked_list(v,head_f_base,j);//riempie gli elementi di posto dispari(indice pari) delle righe con i primi della factor base 		
			if(test){//a è diverso da 1
				if(array_of_prime_chosen_for_a[2*j+1]==1){//array of prime posizione pari è 1
					mpz_set_si(matrix_factor[0][2*j+1],1);//aggiungi esponente di a alla fattorizzazione
				}
			}
			mpz_set_si(matrix_factor[0][2*j],mpz_get_si(v));//scrive il valore del primo della factor base in posizione dispari
			
	}
	print_time_elapsed("inizio creazione matrix_factorization");
	printf("mmm\n");
	//aggiunge alla matrice delle fattorizzazioni il fattore square (aj+b)
	mpz_mul_si(a_temp,a,j);//a_temp=a*j
	mpz_add(a_temp,a_temp,b);//a_temp=a*j+b
	mpz_set(matrix_factor[0][length],a_temp);//imposta alla prima riga (aj+b) dopo di che ai prossimi elementi aggiungi iterativamente a
	j++;
	print_time_elapsed("inizio creazione matrix_factorization");
	printf("qqqq\n");
	for(long i=1;i<2*M+1;i++){
		for(long j=0;j<length;j++){
			mpz_set(matrix_factor[i][j],matrix_factor[0][j]);
		}
		mpz_add(matrix_factor[i][length],matrix_factor[i-1][length],a);//aggiungi a all'elemento precedente
		j++;
	}
	print_time_elapsed("inizio creazione matrix_factorization");
	printf("tttt\n");
	for(long i=0;i<2*M+1;i++){
		if(mpz_cmp_si(matrix_factor[i][length],0)<0){//se gli elementi sono negativi
			mpz_neg(matrix_factor[i][length],matrix_factor[i][length]);//inverti il segno
		}
	}
	print_time_elapsed("inizio creazione matrix_factorization");
	printf("kkkk\n");
	mpz_clear(v);
	mpz_clear(a_temp);
	mpz_clear(j_temp);
	return matrix_factor;
}*/

/*mpz_t**create_matrix_factorization(long M,int cardinality_factor_base,struct node*head_f_base,long*array_of_prime_chosen_for_a,const mpz_t a,const mpz_t b){
	if(M<=0 || cardinality_factor_base<=0 || (array_of_prime_chosen_for_a==NULL && mpz_get_si(a)!=1) || a==NULL || b==NULL){
		handle_error_with_exit("error in parameter create matrix factorization\n");
	}
	char test=((mpz_get_si(a)!=0) && (mpz_get_si(a)!=1));//test==1 se a è contemporaneamente diverso da 1 e da 0
	long length=cardinality_factor_base*2;
	long j=-M;
	mpz_t v,a_temp,j_temp;
	mpz_init(v);
	mpz_init(a_temp);
	mpz_init(j_temp);
	printf("inizio creazione matrix_factorization\n");
	mpz_t**matrix_factor=alloc_matrix_mpz(2*M+1,length+2);//gli ultimi 2 elementi servono per memorizzare il quadrato e il residuo(se residuo 		è 1 il num è b_smooth)
	struct node*p=head_f_base;
	//calcola prima riga matrice_fattorizzazioni
	for(long j=0;j<cardinality_factor_base;j++){
			mpz_set(matrix_factor[0][2*j],p->prime);
			p=p->next;//riempie gli elementi di posto dispari(indice pari) delle righe con i primi della factor base 		
			if(test){//a è diverso da 1
				if(array_of_prime_chosen_for_a[2*j+1]==1){//array of prime posizione pari è 1
					mpz_set_si(matrix_factor[0][2*j+1],1);//aggiungi esponente di a alla fattorizzazione
				}
			}
			
	}
	print_time_elapsed("inizio");
	//aggiunge alla matrice delle fattorizzazioni il fattore square (aj+b)
	mpz_mul_si(a_temp,a,j);//a_temp=a*j
	mpz_add(a_temp,a_temp,b);//a_temp=a*j+b
	mpz_set(matrix_factor[0][length],a_temp);//imposta alla prima riga (aj+b) dopo di che ai prossimi elementi aggiungi iterativamente a
	j++;
	for(long i=1;i<2*M+1;i++){
		for(long j=0;j<length;j++){
			mpz_set(matrix_factor[i][j],matrix_factor[0][j]);
		}
		mpz_add(matrix_factor[i][length],matrix_factor[i-1][length],a);//aggiungi a all'elemento precedente
		j++;
	}
	print_time_elapsed("durante");
	for(long i=0;i<2*M+1;i++){
		if(mpz_cmp_si(matrix_factor[i][length],0)<0){//se gli elementi sono negativi
			mpz_neg(matrix_factor[i][length],matrix_factor[i][length]);//inverti il segno
		}
	}
	mpz_clear(v);
	mpz_clear(a_temp);
	mpz_clear(j_temp);
	print_time_elapsed("fine");
	return matrix_factor;
}*/

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
            handle_error_with_exit("error in divde_elem_by_2\n");
        }
        vector[i]=vector[i]/2;
    }
    return;
}
/*int count_number_potential_B_smooth_matrix_sorted(mpz_t**single_matrix_factorization,int num_row,int num_col){
	//assume che la matrice ha i residui ordinati
	int potential_number_B_smooth=0;
	if(single_matrix_factorization==NULL || *single_matrix_factorization==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in count_number_potential_B_smooth\n");
	}
	if(mpz_cmp_si(single_matrix_factorization[0][num_col-1],1)!=0){//se il primo residuo è diverso da 1 ritorna 0
		return 0;
	}
	int i=1;
	potential_number_B_smooth++;
	while(mpz_cmp_si(single_matrix_factorization[i][num_col-1],1)==0){
		potential_number_B_smooth++;//aumenta il numero di b_smooth potenziali
		i++;//passa alla prossima riga
	}
	return potential_number_B_smooth;
}*/

/*int count_number_B_smooth_matrix_unsorted(mpz_t**single_matrix_factorization,int num_row,int num_col){
	//non assume che la matrice ha i residui ordinati
	int number_B_smooth=0;
	if(single_matrix_factorization==NULL || *single_matrix_factorization==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in count_number_potential_B_smooth\n");
	}
	int i=0;
	while(i<num_row){
		if(mpz_cmp_si(single_matrix_factorization[i][num_col-1],1)==0){
			number_B_smooth++;//aumenta il numero di b_smooth potenziali
		}
		i++;//passa alla prossima riga
	}
	return number_B_smooth;
}*/

/*int count_number_B_smooth_matrix_unsorted_f(struct matrix_factorization *single_matrix_factorization,int num_row){
	//non assume che la matrice ha i residui ordinati
	if(num_row<=0 || single_matrix_factorization==NULL){
		handle_error_with_exit("error in count number_B_smooth_matrix_unsorted_f\n");
	}
	int number_B_smooth=0;
	if(num_row<=0){
		handle_error_with_exit("error in count_number_potential_B_smooth\n");
	}
	int i=0;
	double ratio;
	while(i<num_row){
		if((*single_matrix_factorization).row[i].log>(*single_matrix_factorization).row[i].sum_log){
			ratio=(*single_matrix_factorization).row[i].sum_log/(*single_matrix_factorization).row[i].log;
		}
		else{
			ratio=(*single_matrix_factorization).row[i].log/(*single_matrix_factorization).row[i].sum_log;
		}
		if(ratio>=PERCENT_B_SMOOTH){
			number_B_smooth++;//aumenta il numero di b_smooth potenziali
		}
		i++;//passa alla prossima riga
	}
	return number_B_smooth;
}*/
/*void copy_residuos_from_matrix(mpz_t**potential_matrix_B_smooth,int potential_num_B_smooth,mpz_t**array_matrix_mpz,int num_row_matrix,int num_col_matrix){
	//copia le righe B_smooth della matrice e le mette in potential_matrix_B_smooth
	if(potential_matrix_B_smooth==NULL || *potential_matrix_B_smooth==NULL || potential_num_B_smooth<=0 || 
		array_matrix_mpz==NULL || *array_matrix_mpz==NULL || num_row_matrix<=0 || num_col_matrix<=0){
		handle_error_with_exit("error in copy_residuos_from_matrix\n");
	}
	int count=0;
	for(int i=0;i<num_row_matrix;i++){
		if(mpz_cmp_si(array_matrix_mpz[i][num_col_matrix-1],1)==0){//se il residuo della riga i è 1
			copy_array_mpz(potential_matrix_B_smooth[count],array_matrix_mpz[i],num_col_matrix);
			count++;
		}
	}
	if(count!=potential_num_B_smooth){
		handle_error_with_exit("error in copy_residuos_from_matrix\n");
	}
	return;
}*/


/*void copy_row_matrix(struct row*row_dest,struct row row_source){
	if(row_dest==NULL){
		handle_error_with_exit("error in copy_row_matrix\n");
	}
	row_dest->index_first_prime=row_source.index_first_prime;
	row_dest->index_last_prime=row_source.index_last_prime;
	mpz_set(row_dest->num,row_source.num);
	mpz_set(row_dest->square,row_source.square);
	row_dest->sum_log=row_source.sum_log;
	row_dest->log=row_source.log;
	return;
}


mpz_t** create_matrix_B_smooth(int num_B_smooth,mpz_t**array_matrix_mpz,int num_row_matrix,int num_col_matrix){
	//copia le righe B_smooth della matrice e le mette in potential_matrix_B_smooth
	if(num_B_smooth<=0 || array_matrix_mpz==NULL || *array_matrix_mpz==NULL || num_row_matrix<=0 || num_col_matrix<=0){
		handle_error_with_exit("error in create_matrix_B_smooth\n");
	}
	mpz_t**matrix_B_smooth=alloc_matrix_mpz(num_B_smooth,num_col_matrix-1);
	int count=0;
	for(int i=0;i<num_row_matrix;i++){
		if(mpz_cmp_si(array_matrix_mpz[i][num_col_matrix-1],1)==0){//se il residuo della riga i è 1
			copy_array_mpz(matrix_B_smooth[count],array_matrix_mpz[i],num_col_matrix-1);//copia tutto tranne il residuo
			count++;
		}
	}
	if(count!=num_B_smooth){
		handle_error_with_exit("error in create_matrix_B_smooth\n");
	}
	return matrix_B_smooth;
}
struct matrix_factorization*create_matrix_B_smooth_f(int num_B_smooth,struct matrix_factorization mat,int num_row_matrix,int num_col_matrix){
	//copia le righe B_smooth della matrice e le mette in potential_matrix_B_smooth
	if(num_B_smooth<=0 || num_row_matrix<=0 || num_col_matrix<=0){
		handle_error_with_exit("error in create_matrix_B_smooth\n");
	}
	char B_smooth=0;
	double ratio;
	int count=0;
	struct matrix_factorization*matrix_B_smooth=alloc_matrix_factorization(num_B_smooth);//-1 non copia il residuo
	for(int i=0;i<num_row_matrix;i++){
		B_smooth=0;
		if(mat.row[i].log>mat.row[i].sum_log){
			ratio=mat.row[i].sum_log/mat.row[i].log;
		}
		else{
			ratio=mat.row[i].log/mat.row[i].sum_log;
		}
		if(ratio>=PERCENT_B_SMOOTH){
			B_smooth=1;//aumenta il numero di b_smooth potenziali
		}
		if(B_smooth==1){//se il residuo della riga i è 1
			copy_row_matrix(&(*matrix_B_smooth).row[count],mat.row[i]);
			count++;
		}
	}
	if(count!=num_B_smooth){
		handle_error_with_exit("error in create_matrix_B_smooth\n");
	}
	return matrix_B_smooth;
}*/

mpz_t*create_array_temp_factorization(int card_f_base,mpz_t*row_matrix_b_smooth){// crea array -1 0 2 0 30...primo factor base 0...di 2*cardinalità 		della factor base elementi
	if(card_f_base<=0 || row_matrix_b_smooth==NULL){
		handle_error_with_exit("error in parameter create_array_temp_factorization\n");
	}
	mpz_t*v;
	v=alloc_array_mpz(card_f_base*2);//alloca memoria
	for(int i=0;i<card_f_base*2;i=i+2){//ciclo a salti di 2
		mpz_set(v[i],row_matrix_b_smooth[i]);//sui posti dispari (indice pari) viene messo il primo della factor base
	}
	return v;
}
int**alloc_linear_system(int cardinality_f_base,int num_B_smooth){//alloca sistema lineare
	if(cardinality_f_base<=0 || num_B_smooth<=0 || num_B_smooth<cardinality_f_base){
		handle_error_with_exit("invalid parameter alloc linear system\n");
	}
	int**linear_system;
	linear_system=alloc_matrix_int(cardinality_f_base,num_B_smooth);
	return linear_system;
}
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
int**create_linear_system(int cardinality_f_base,mpz_t** matrix_B_smooth,int num_B_smooth){//alloca matrice con la fattorizzazione dei numeri b_smmoth
	if(cardinality_f_base<=0 || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0){
		handle_error_with_exit("error in create_linear_system\n");
	}
	int**linear_system=alloc_linear_system(cardinality_f_base,num_B_smooth);
	for(int i=0;i<cardinality_f_base;i++){
		for(int j=0;j<num_B_smooth;j++){
			if(mpz_divisible_2exp_p(matrix_B_smooth[j][2*i+1],1)==0){//non è divisibile per 2
				linear_system[i][j]=1;
			}
			else{//è divisibile per 2
				linear_system[i][j]=0;
			}
			//linear system[i][j]=esponente del'iesimo elemento della factor base preso dalla fattorizzazione j-esimo numero B-smooth
		}	
	}
	return linear_system;
}

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
int calculate_dim_sol_char(char*matrix,int num_row,int num_col){//matrice deve essere ridotta a scala
    if(matrix==NULL || num_row<=0 || num_col<=0){
        handle_error_with_exit("error in parameter get_col\n");
    }
    int num_row_not_null=count_rows_not_null_char(matrix,num_row,num_col);
    int num_col_not_null=count_cols_not_null_char(matrix,num_row,num_col);
    return num_col_not_null-num_row_not_null;//numero incognite-rango matrice
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

/*void calculate_vector_base(int **matrix_linear_system,int num_row,int num_col,char*array_var,int*v){//calcola un vettore di base soluzione del sistema lineare,matrix linear system è ridotta a scala
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
	double sum_negative=0;
	for(int i=0;i<num_row-rows_null;i++){//ciclo per tutta la sottomatrice con righe non nulle
		sum_negative=0;
		index_to_calculate=-1;
		for(int j=0;j<num_col;j++){//per ogni riga trova il valore da calcolare
			if(matrix_linear_system[num_row_not_null-1-i][num_col-1-j]==0){//se l'elemento della riga non nulla è 0 passa al prossimo 					elemento della riga
				continue;
			}
			if(array_var[num_col-1-j]==1 || array_var[num_col-1-j]==3){//se liberi o già calcolati(entrambi sono già stati assegnati)
				sum_negative=sum_negative-matrix_linear_system[num_row_not_null-1-i][num_col-1-j]*v[num_col-1-j];
			}
			else{//array_var[num_col-1-j]!=1 || array_var[num_col-1-j]!=3 ->array_var==2==ricavabile
				index_to_calculate=num_col-1-j;//memorizza la colonna in cui c'è la variabile da calcolare,
			}
		}
		if(index_to_calculate<0){//riga in cui non è stato trovato nessun elemento da calcolare
			printf("index_to_calculate=%d\n",index_to_calculate);
			print_array_char(array_var,num_col);
			handle_error_with_exit("error in index_calculate index<0\n");
		}
		if(index_to_calculate>num_col){
			printf("num_col=%d,index_to_calculate=%d\n",num_col,index_to_calculate);
			handle_error_with_exit("error in index_calculate,index>num_col\n");
		}
		if(num_row_not_null-1-i>num_row || num_row_not_null-1-i<0){
			handle_error_with_exit("error in index calculate vector base\n");
		}
	  	sum_negative=sum_negative/matrix_linear_system[num_row_not_null-1-i][index_to_calculate];//calcola il valore della variabile
		v[index_to_calculate]=sum_negative;//imposta il valore della variabile nel vettore di base
		array_var[index_to_calculate]=3;//imposta variabile come calcolata,piano piano le variabili diventano tutte calcolate
	}
	return;
}*/

void calculate_vector_base(int **matrix_linear_system,int num_row,int num_col,char*array_var,int*v){//calcola un vettore di base soluzione del sistema lineare,matrix linear system è ridotta a scala e binaria,calcola vettori di base binari
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
			print_array_char(array_var,num_col);
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
void calculate_vector_base_char(char*matrix_linear_system,int num_row,int num_col,char*array_var,int*v){//calcola un vettore di base soluzione del sistema lineare,matrix linear system è ridotta a scala e binaria,calcola vettori di base binari
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
    int num_row_not_null=count_rows_not_null_char(matrix_linear_system,num_row,num_col);
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
            print_array_char(array_var,num_col);
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
char* find_free_var_char(char*matrix_linear_system,int num_row,int num_col){//matrice ridotta a scala,ritorna array con elementi 1 e 2 (ricordiamo che 1=variabile libera 2=variabile ricavabile 3=variabile calcolata)
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
	int num_row_not_null=count_rows_not_null_char(matrix_linear_system,num_row,num_col);
	int rows_null=num_row-num_row_not_null;
	int num_free_var_left=calculate_dim_sol_char(matrix_linear_system,num_row,num_col);//il numero delle variabili libere è lo stesso della dimensione della base
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
/*int**calculate_base_linear_system(int**matrix_linear_system,int num_row,int num_col,int*dim_sol){//matrice ridotta modulo n,calcola una base del sistema lineare
	if(matrix_linear_system==NULL || *matrix_linear_system==NULL || num_row<=0 || num_col<=0 || dim_sol==NULL){
		handle_error_with_exit("error in parameter get_coli\n");
	}
	reduce_echelon_form(matrix_linear_system,num_row,num_col);//riduce la matrice a scala
	if(check_if_matrix_is_reduce_mod_n(matrix_linear_system,num_row,num_col,2)==0){
		handle_error_with_exit("matrix is not reduce mod n");
	}
	if(check_if_matrix_is_echelon_reduce(matrix_linear_system,num_row,num_col)==0){
		handle_error_with_exit("error in calculate_base_linear_system\n");
	}
	//printf("matrice_ridotta\n");
	//print_matrix_int(matrix_linear_system,num_row,num_col);
	*dim_sol=calculate_dim_sol(matrix_linear_system,num_row,num_col);//calcola la dimensione della base del sistema lineare
	fprintf(file_log,"dim_sol=%d ",*dim_sol);
	if(*dim_sol<0){
		handle_error_with_exit("error in calculate dim sol linear system\n");
	}
	if(*dim_sol==0){
		return NULL;
	}
	int free_var=*dim_sol;
	if(*dim_sol>MAX_DIM_SOL){
		*dim_sol=MAX_DIM_SOL;
	}
	int**base_linear_system=alloc_matrix_int(num_col,*dim_sol);//la base del sistema lineare ha come righe il numero di colonne della matrice 		(numero delle variabili) e come colonne il numero di vettori linearmente indipendenti
	char*array_var=NULL;
	array_var=find_free_var(matrix_linear_system,num_row,num_col);//calcola array delle variabili libere che specifica tutte le 				variabili che sono state impostate come libere per tutto il sistema,è lungo num_col,è necessario 				allocarlo ogni volta perchè viene sporcato
		//dalla funzione calculate_vector_base
	if(check_if_array_var_is_correct(array_var,num_col,free_var)==0){
			handle_error_with_exit("error in calculate array_var\n");
	}
	if(check_if_array_var_is_valid(matrix_linear_system,num_row,num_col,array_var)==0){
		printf("array :");
		print_array_char(array_var,num_col);
		handle_error_with_exit("array_var is not valid\n");
	}
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
			calculate_vector_base(matrix_linear_system,num_row,num_col,array_var_temp,v);//memorizza in v la soluzione
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
}*/
int**calculate_base_linear_system_char(char*matrix_linear_system,int num_row,int num_col,int*dim_sol){//matrice ridotta modulo n,calcola una base del sistema lineare
    if(matrix_linear_system==NULL || num_row<=0 || num_col<=0 || dim_sol==NULL){
        handle_error_with_exit("error in parameter get_coli\n");
    }
    reduce_echelon_form_char(matrix_linear_system,num_row,num_col);//riduce la matrice a scala
    if(check_if_matrix_char_is_reduce_mod_n(matrix_linear_system,num_row,num_col,2)==0){
        handle_error_with_exit("matrix is not reduce mod n");
    }
    if(check_if_matrix_char_is_echelon_reduce(matrix_linear_system,num_row,num_col)==0){
        handle_error_with_exit("error in calculate_base_linear_system\n");
    }
    printf("matrice_ridotta a scala\n");
    print_linear_system(matrix_linear_system,num_row,num_col);
    *dim_sol=calculate_dim_sol_char(matrix_linear_system,num_row,num_col);//calcola la dimensione della base del sistema lineare
    printf("dim_sol=%d\n",*dim_sol);
    fprintf(file_log,"dim_sol=%d ",*dim_sol);
    if(*dim_sol<0){
        handle_error_with_exit("error in calculate dim sol linear system\n");
    }
    if(*dim_sol==0){
        return NULL;
    }
    int free_var=*dim_sol;
    if(*dim_sol>MAX_DIM_SOL){
        *dim_sol=MAX_DIM_SOL;
    }
    int**base_linear_system=alloc_matrix_int(num_col,*dim_sol);//la base del sistema lineare ha come righe il numero di colonne della matrice 		(numero delle variabili) e come colonne il numero di vettori linearmente indipendenti
    char*array_var=NULL;
    array_var=find_free_var_char(matrix_linear_system,num_row,num_col);//calcola array delle variabili libere che specifica tutte le 				variabili che sono state impostate come libere per tutto il sistema,è lungo num_col,è necessario 				allocarlo ogni volta perchè viene sporcato
    //dalla funzione calculate_vector_base
    if(check_if_array_var_is_correct(array_var,num_col,free_var)==0){
        handle_error_with_exit("error in calculate array_var\n");
    }
    if(check_if_array_var_is_valid_char(matrix_linear_system,num_row,num_col,array_var)==0){
        printf("array :");
        print_array_char(array_var,num_col);
        handle_error_with_exit("array_var is not valid\n");
    }
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
            calculate_vector_base_char(matrix_linear_system,num_row,num_col,array_var_temp,v);//memorizza in v la soluzione
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

int*get_coli(int **matrix,int num_row,int num_col,int index_col){//indice parte da 0,ottiene la colonna iesima della matrice
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
/*char find_pivot_not_null_in_column(int **matrix,int num_row,int num_col,int pivot[2]){//ritorna l'indice del primo elemento non nullo in 		colonna,partendo da un indice assegnato(pivot[0][1]) 
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 || pivot[0]>=num_row || pivot[1]>=num_col){
		handle_error_with_exit("error in parameter find_pivot\n");
	}
	int*col=get_coli(matrix,num_row,num_col,pivot[1]);
	int start=pivot[0];//inizia dalla riga pivot[0] a scansionare
	int result=scan_array_to_find_element_not_null(col,start,num_row);//scansiona colonna
	if(result==-1){
		free(col);
		col=NULL;
		return result;//ritorna -1==not found
	}
	pivot[0]=result;//imposta pivot[0] a result
	free(col);
	col=NULL;
	return result;
}*/
void reduce_echelon_form(int**matrix,int num_row,int num_col){//versione rref,fare il test per verificare
		// che è effettivamente ridotta in modo rref
	int lead=0;
	int i;
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	if(check_if_matrix_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
		handle_error_with_exit("matrix is not reduce mod n");
	}
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
		if(matrix[r][lead]==0){
			handle_error_with_exit("invalid pivot\n");
		}
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
		//if(binary_matrix[r][lead]==0){
		//	handle_error_with_exit("invalid pivot\n");
		//}

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
    if(matrix==NULL || num_row<=0 || num_col<=0){
        handle_error_with_exit("error in parameter reduce echelon form\n");
    }
    if(check_if_matrix_char_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
        handle_error_with_exit("matrix is not reduce mod n");
    }
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

/*void reduce_echelon_form(int**matrix,int num_row,int num_col){//versione precedente della riduzione a scala,versione ref
	//conviene implementare rref?
	int lead=0;
	int i;
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	if(check_if_matrix_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
		handle_error_with_exit("matrix is not reduce mod n");
	}
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
		if(matrix[r][lead]==0){
			handle_error_with_exit("invalid pivot\n");
		}
		for(int i=r+1;i<num_row;i++){//riduzione della matrice
			if(matrix[i][lead]!=0){
				for(int f=0;f<num_col;f++){//in realtà il ciclo si può fare da lead a num_col,la parte prima di lead sono 						sottrazioni per 0
					matrix[i][f]=matrix[i][f]-matrix[r][f];//si fa la differenza perchè il pivot è sempre 1
					matrix[i][f]=reduce_mod_2(matrix[i][f]);
				}
			}
		}
		lead=lead+1;
	}
	return;
}*/

/*void reduce_echelon_form(int**matrix,int num_row,int num_col){//riduce la matrice ridotta modulo 2 a scala
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0 ){
		handle_error_with_exit("error in parameter reduce echelon form\n");
	}
	if(check_if_matrix_is_reduce_mod_n(matrix,num_row,num_col,2)==0){
		handle_error_with_exit("error in reduce echelon form\n");
	}
	int pivot[2],index_actual_row=0;//pivot[0]=riga in cui è il pivot pivot[1]=colonna in cui è il pivot
	char founded;
	//parti dal pivot in alto a sinistra della matrice
	pivot[0]=0;
	pivot[1]=0;
	while(pivot[0]<num_row && pivot[1]<num_col){//finquando il pivot non è sceso fino alla fine della matrice
		if(matrix[pivot[0]][pivot[1]]==0){//se pivot è zero cerca un altro pivot e scambia le righe
			founded=find_pivot_not_null_in_column(matrix,num_row,num_col,pivot);//pivot[0] viene impostato al valore non nullo della 				colonna pivot[1]
			
			if(founded==-1){//pivot not found in column
				pivot[1]=pivot[1]+1;//vai alla colonna successiva,mantieni invariata la riga senza cambiare pivot[0]
				continue;//continua a vedere se il pivot trovato è !=0
				//questa procedura garantisce che ci sia nella riga i-esima che si sta considerando un pivot,
				//anche dopo molti elementi uguali a zero
			}
			else{//pivot founded->scambia le righe e riassegna il pivot[0]
				swap_row(matrix,num_row,num_col,index_actual_row,pivot[0]);//scambia le righe con il pivot trovato
				pivot[0]=index_actual_row;//pivot[0] impostato a indice di riga attuale
				//questa procedura ritorna sicuramente un pivot preso da matrix[pivot[0]]matrix[pivot[1]]!=0
			}
		}
		//il pivot è stato assegnato,sottrai alla sottomatrice m volte la riga dove è il pivot per annullare tutti gli elementi nella 			colonna del pivot
		for(int i=pivot[0]+1;i<num_row;i++){//dalla riga successiva a fine matrice
			if(matrix[i][pivot[1]]==0){//se nella colonna del pivot l'elemento è nullo vai alla riga successiva,
				continue;
			}
			double mult=matrix[i][pivot[1]]/matrix[pivot[0]][pivot[1]];//moltiplicatore
			for(int j=pivot[1];j<num_col;j++){//dalla colonna del pivot fino a fine matrice
				matrix[i][j]=matrix[i][j]-mult*matrix[pivot[0]][j];
				reduce_int_mod_n(&(matrix[i][j]),2);//ridurre elemento modificato mod 2
			}
		}
		//una volta ridotti a zero tutti gli elementi nella colonna del pivot fai scendere il pivot di una riga e di una colonna e 			ri-assegna index_actual_row
		pivot[0]=pivot[0]+1;//aumenta una riga del pivot
		pivot[1]=pivot[1]+1;//aumenta una colonna del pivot
		index_actual_row=pivot[0];//indice di riga attuale diventa quella del pivot[0]
		//reduce_matrix_mod_n(matrix,num_row,num_col,2);//senza ridurre la matrice mod n si ottiene una riduzione a scala errata
		//ridurre la matrice modulo n viola la riduzione a scala del sistema?
	}
	return;
}*/


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
char*create_linear_system_f(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth){
	if(head==NULL || cardinality_factor_base<=0 || num_B_smooth<=0){
		handle_error_with_exit("error in create_linear_system\n");
	}
	char*linear_system=alloc_array_char(cardinality_factor_base*num_B_smooth);
	struct node_square_relation*p=head;
	int col_index=0;
	long index;
	while(p!=NULL){//cicla su tutte le relazioni quadratiche
		struct node_factorization*sq=p->square_relation.head_factorization;
		while(sq!=NULL){//cicla su tutti i fattori della singola relazione quadratica
			if((sq->exp_of_number & 1)!=0){//se l'esponente non è divisibile per 2
				//metti 1 in posizione indice nella colonna iesima
				index=get_index(sq->index,col_index,num_B_smooth);
				linear_system[index]=1;
			}
			sq=sq->next;
		}
		p=p->next;//passa alla prossima relazione quadratica
		col_index++;
	}
	if(col_index!=num_B_smooth){
		handle_error_with_exit("error initialize_linear_system\n");
	}

	return linear_system;
}
unsigned long**create_binary_linear_system(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth,int*num_col_binary_matrix){
    if(head==NULL || cardinality_factor_base<=0 || num_B_smooth<=0 || num_col_binary_matrix==NULL){
        handle_error_with_exit("error in create_linear_system\n");
    }
    int remainder=num_B_smooth%BIT_OF_UNSIGNED_LONG;
    if(remainder==0){
        *num_col_binary_matrix=num_B_smooth/BIT_OF_UNSIGNED_LONG;
    }
    else{
        *num_col_binary_matrix=((num_B_smooth-remainder)/BIT_OF_UNSIGNED_LONG)+1;
    }
    printf("bit unsigned long=%lu\n",BIT_OF_UNSIGNED_LONG);
    printf("num_col=%d\n",*num_col_binary_matrix);
    unsigned long**binary_linear_system=alloc_matrix_unsigned_long(cardinality_factor_base,*num_col_binary_matrix);
    struct node_square_relation*p=head;
    int col_index=0;
    while(p!=NULL){//cicla su tutte le relazioni quadratiche
        struct node_factorization*sq=p->square_relation.head_factorization;
        while(sq!=NULL){//cicla su tutti i fattori della singola relazione quadratica
            if((sq->exp_of_number & 1)!=0){//se l'esponente non è divisibile per 2
                //metti 1 in posizione indice nella colonna iesima
                int remainder=col_index%BIT_OF_UNSIGNED_LONG;//va da 0 a bit_unsigned_long->ci dice qual'è l'indice del bit da mettere a 1
                int index_element=col_index/BIT_OF_UNSIGNED_LONG;
                binary_linear_system[sq->index][index_element]^= 1UL << (BIT_OF_UNSIGNED_LONG-remainder-1);//toggle bit
            }
            sq=sq->next;
        }
        p=p->next;//passa alla prossima relazione quadratica
        col_index++;
    }
    if(col_index!=num_B_smooth){
        handle_error_with_exit("error initialize_linear_system\n");
    }

    return binary_linear_system;
}
/*char*create_linear_system_f(struct matrix_factorization *mat,int cardinality_factor_base){
	//mat contiene solamente i B_smooth;
	if(mat->num_row<=0 || cardinality_factor_base<=0 || mat==NULL){
		handle_error_with_exit("error in create linear_system_f\n");
	}
	int*array_factorization=NULL;
	int num_of_prime;
	char is_B_smooth=0;
	char*linear_system=alloc_array_char(mat->num_row*cardinality_factor_base);
	for(int i=0;i<mat->num_row;i++){
		num_of_prime=(mat->row[i].index_last_prime-mat->row[i].index_first_prime+1);
		array_factorization=alloc_array_int(num_of_prime*2);
		gmp_printf("num %Zd\n",mat->row[i].num);
		is_B_smooth=factorize_num_B_smooth(array_factorization,num_of_prime*2,&(mat->row[i]));
		print_array_int(array_factorization,num_of_prime*2);
		if(is_B_smooth==1 || is_B_smooth==2){
			if(is_B_smooth==1){
				//linear_system[0][i]=1;
			}
			//le prime colonne e le ultime rimangono riempite di zeri
			for(int j=0;j<num_of_prime;j++){//farlo in modo tale che non serve sapere gli indici dei primi che
				//compaiono nella fattorizzazione
				//riempire la colonna iesima con il primo esponente ridotto modulo 2 poi con il secondo esponente poi con il terzo 					fino alla fine
				//dell'array,gli altri valori saranno tutti 0
			}
		}
		free(array_factorization);
		array_factorization=NULL;
	}
	return linear_system;
}*/







