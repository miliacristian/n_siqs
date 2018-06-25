#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "miller_rabin.h"
#include "dynamic_list.h"
#include <pthread.h>
#include "matrix_function.h"
#include "main.h"
#include "list_square_relation.h"
#include "list_factorization.h"

struct timespec;
struct timespec timer;
struct timespec time_start;
FILE*file_log;
extern int M;
extern int cardinality_factor_base;
void destroy_mtx(pthread_mutex_t *mtx){
	if(mtx==NULL){
		handle_error_with_exit("error in destroy_mtx mtx is NULL\n");
	}
	if(pthread_mutex_destroy(mtx)!=0){
		handle_error_with_exit("error in pthread_mutex_destroy\n");
	}
	return;
}
void initialize_mtx(pthread_mutex_t *mtx){//inizializza mutex
	if(mtx==NULL){
		handle_error_with_exit("error in initialize_mtx mtx is NULL\n");
	}
	if(pthread_mutex_init(mtx,NULL)!=0){
		handle_error_with_exit("error in initialize mtx\n");
	}
	return;
}
void lock_mtx(pthread_mutex_t *mtx){//lock mutex
	if(mtx==NULL){
		handle_error_with_exit("error in lock_mtx mtx is NULL\n");
	}
	if(pthread_mutex_lock(mtx)!=0){
		handle_error_with_exit("error in pthread_mutex_lock\n");
	}
	return;
}
void unlock_mtx(pthread_mutex_t *mtx){//unlock mutex
	if(mtx==NULL){
		handle_error_with_exit("error in unlock_mtx mtx is NULL\n");
	}
	if(pthread_mutex_unlock(mtx)!=0){
		handle_error_with_exit("error in pthread_mutex_unlock\n");
	}
	return;
}
int compare_a_struct( const void* a, const void* b)
{
	struct a_struct s_a = * ( (struct a_struct*) a );
	struct a_struct s_b = * ( (struct a_struct*) b );
	if(s_a.index_prime_a==s_b.index_prime_a){
		return 0;
	}
	else if(s_a.index_prime_a>s_b.index_prime_a){
		return 1;
	}
	return -1;
}

void free_list_factorization(struct node_factorization*head_factorization){
    if(head_factorization==NULL){
    	return;
    }
    struct node_factorization*p=head_factorization;
    struct node_factorization*q;
    while(p!=NULL){
        q=p->next;
        free(p);
        p=q;
    }
    return;
}
void free_memory_list_square_relation(struct node_square_relation*head){
    if(head==NULL){
    	return;
    }
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL){
        q=p->next;
        mpz_clear(p->square_relation.square);
        mpz_clear(p->square_relation.num);
        mpz_clear(p->square_relation.residuos);
        free_list_factorization(p->square_relation.head_factorization);
        free(p);
        p=q;
    }


    return;
}

void free_array_thread_data(struct thread_data*thread_data,int length_array_thread_data){
	if(thread_data==NULL || length_array_thread_data<=0){
		handle_error_with_exit("error in free_array_thread_data\n");
	}
	for(int i=0;i<length_array_thread_data;i++){
		mpz_clear(thread_data[i].b);
		thread_data[i].head=NULL;
		thread_data[i].tail=NULL;
		free(thread_data[i].numbers);
	}
	free(thread_data);
	return;
}
void clear_struct_thread_data(struct thread_data t_data,int M) {
	for (int i = 0; i < 2 * M + 1; i++) {
		//t_data.numbers[i].j = i - M;
		t_data.numbers[i].first_index_f_base = -1;
		t_data.numbers[i].sum_log = 0;
		t_data.numbers[i].last_index_f_base = -1;
	}
	return;
}
struct factor_base_data*alloc_array_factor_base_data(int length){
	if(length<=0){
		handle_error_with_exit("error in alloc_array_factor_base_data\n");
	}
	struct factor_base_data*array_factor_base=malloc(sizeof(struct factor_base_data)*length);
	if(array_factor_base==NULL){
		handle_error_with_exit("error in malloc array_factor_base\n");
	}
	memset(array_factor_base,0, sizeof(struct factor_base_data)*length);
	return array_factor_base;
}

struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M){
	if(length_array_thread_data<=0 || M<=0){
		handle_error_with_exit("error in alloc_array_thread_data\n");
	}
	struct thread_data*t_data=malloc(sizeof(struct thread_data)*length_array_thread_data);
	if(t_data==NULL){
		handle_error_with_exit("error in malloc alloc array thread_data\n");
	}
	for(int i=0;i<length_array_thread_data;i++){
		mpz_init(t_data[i].b);
		t_data[i].log_thresold=0;
		t_data[i].num_potential_B_smooth=0;
		t_data[i].num_B_smooth=0;
		t_data[i].num_semi_B_smooth=0;
		t_data[i].head=NULL;
		t_data[i].tail=NULL;
		t_data[i].numbers=malloc(sizeof(struct number)*(2*M+1));
		if(t_data[i].numbers==NULL){
			handle_error_with_exit("error in malloc alloc numbers");
		}
	}
	for(int i=0;i<2*M+1;i++){
		t_data[0].numbers[i].j=i-M;
		t_data[0].numbers[i].first_index_f_base=-1;
		t_data[0].numbers[i].sum_log=0;
		t_data[0].numbers[i].last_index_f_base=-1;
	}
	for(int i=0;i<length_array_thread_data;i++){
		t_data[i].log_thresold=0;
		memcpy(t_data[i].numbers,t_data[0].numbers,sizeof(struct number)*(2*M+1));
	}
	return t_data;
}

FILE*open_file_log(){
	FILE*file=fopen("file_log.txt","a");//modalità append
	if(file==NULL){
		handle_error_with_exit("error in fopen,impossible open file\n");
	}
	return file;
}
void test(){
	//operazioni da fare in fase di test
	handle_error_with_exit("\n");
}

void*thread_job_polynomial(void*arg){
	int*id=arg;
	thread_job_criv_quad(*id);
	return NULL;
}
void*thread_job_factor_base(void*arg){
	int*id=arg;
	thread_job_to_create_factor_base(*id);
	pthread_exit(NULL);
	return NULL;
}

int* create_threads(pthread_t*array_tid,int num_thread){//ritorna il numero di thread creati con successo
	if(num_thread<0 || array_tid==NULL){
		handle_error_with_exit("error in create_thread\n");
	}
	if(num_thread==0){
		return NULL;
	}
	int*array_id=alloc_array_int(num_thread);
	for(int i=0;i<num_thread;i++){
		array_id[i]=i;
		if(pthread_create(&(array_tid[i]),NULL,thread_job_polynomial,&(array_id[i]))!=0){
			handle_error_with_exit("error in pthred_create create_thread\n");
		}
	}
	return array_id;
}
struct factorization_thread_data* create_factorization_threads(pthread_t*array_tid,struct thread_data thread_data,const mpz_t a,int num_thread){//ritorna il numero di thread creati con successo
    if(num_thread<0 || array_tid==NULL){
        handle_error_with_exit("error in create_thread\n");
    }
    if(num_thread==0){
        return NULL;
    }
    struct factorization_thread_data*factorization_data=malloc(sizeof(struct factorization_thread_data)*num_thread);
    if(factorization_data==NULL){
    	handle_error_with_exit("error in malloc factorization_data\n");
    }
    struct thread_data*t=alloc_array_polynomial_thread_data(num_thread,M);
    factorization_data[0].pointer=t;
	long remainder=reduce_int_mod_n_v2(cardinality_factor_base,num_thread);
	long length=(cardinality_factor_base-remainder)/(num_thread);
    for(int id=0;id<num_thread;id++){
        int start=id*length;
        int end=start+length-1;
        if(id==num_thread-1){//se è l'ultimo thread aggiungi anche il resto
        	end+=remainder;
        }
        if(start>end){
        	handle_error_with_exit("error in calculate start,end\n");
        }
        factorization_data[id].thread_data=t[id];
        mpz_set(factorization_data[id].thread_data.b,thread_data.b);
        factorization_data[id].id_thread=id;//assegna ad ogni thread l'indice
        factorization_data[id].start=start;
        factorization_data[id].end=end;
        if(mpz_cmp_si(a,1)==0) {
			factorization_data[id].is_a_default = 1;
		}
		else{
			factorization_data[id].is_a_default = 0;
        }
        if(pthread_create(&(array_tid[id]),NULL,thread_factorization_job,&(factorization_data[id]))!=0){
            handle_error_with_exit("error in pthread_create create_thread\n");
        }
    }
    return factorization_data;
}

int* create_factor_base_threads(pthread_t*array_tid,int num_thread){//ritorna il numero di thread creati con successo
	if(num_thread<0 || array_tid==NULL){
		handle_error_with_exit("error in create_thread\n");
	}
	if(num_thread==0){
		return NULL;
	}
	int*array_id=alloc_array_int(num_thread);
	for(int i=0;i<num_thread;i++){
		array_id[i]=i;
		if(pthread_create(&(array_tid[i]),NULL,thread_job_factor_base,&(array_id[i]))!=0){
			handle_error_with_exit("error in pthred_create create_thread\n");
		}
	}
	return array_id;
}
void join_all_threads(pthread_t*array_tid,int length_array){
	if(array_tid==NULL){
		return;//no thread to wait
	}
	if(length_array<=0){
		handle_error_with_exit("error in join all_thread\n");
	}
	for(int i=0;i<length_array;i++){
		if(pthread_join(array_tid[i],NULL)!=0){
			handle_error_with_exit("error in pthread_join join all_thread\n");
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
	memset(array_tid,0,sizeof(pthread_t)*num_thread);
	return array_tid;
}


long get_file_size(char*path){//ritorna la dimensione di un file dato un path
    if(path==NULL){
        handle_error_with_exit("error in get_file_size\n");
    }
    struct stat st;
    if(stat(path,&st)==-1){
        handle_error_with_exit("error get file_status\n");
    }
    return st.st_size;
}
FILE*open_file(char*path){
	if(path==NULL){
		handle_error_with_exit("error in open_file path is NULL\n");
	}
	FILE*file=fopen(path,"r");
	if(file==NULL){
		handle_error_with_exit("error in opne_file,fopen function\n");
	}
	return file;
}
char get_and_check_n(mpz_t n,FILE*file_number){
	if(file_number==NULL || n==NULL){
		handle_error_with_exit("error in get and check n,invalid filename\n");
	}
	int digit=mpz_inp_str(n,file_number,10);//leggi dal file n in base 10
	if(mpz_sgn(n)<0){
		mpz_neg(n,n);//n=-n;//inverte il segno di n
	}
	if(mpz_cmp_si(n,0)==0 || mpz_cmp_si(n,1)==0){//se n==0 o n==1 errore
		handle_error_with_exit("n must be different than 0,1 and -1\n");
	}
	test_n(n,NUM_TEST_MILLER_RABIN);//verifica che n non è un numero primo
	return digit;
}


void gettime(struct timespec*timer){
	if(timer==NULL){
		handle_error_with_exit("error in gettime\n");
	}
	if(clock_gettime(CLOCK_MONOTONIC,(struct timespec*)timer)!=0){
        	handle_error_with_exit("error in clock gettime\n");
    	}
}

struct timespec diff_timespec(struct timespec time_current,struct timespec timer){
	//timecurrent>timer
	struct timespec time_sub;
	if(time_current.tv_nsec>=timer.tv_nsec){
		time_sub.tv_nsec=time_current.tv_nsec-timer.tv_nsec;
		time_sub.tv_sec=time_current.tv_sec-timer.tv_sec;
	}
	else{//riporto nanosecondi di time_current minori di nanosecondi time
		time_current.tv_sec-=1;
		timer.tv_nsec-=time_current.tv_nsec;
		time_current.tv_nsec=0;//evita l'overflow nei nanosecondi
		time_current.tv_nsec+=1000000000;
		time_sub.tv_nsec=time_current.tv_nsec-timer.tv_nsec;
		time_sub.tv_sec=time_current.tv_sec-timer.tv_sec;
	}
	return time_sub;
}
void print_time_elapsed(char*string){
	//ns=1000000000=1 sec
	//ns=1000000 1 ms
	long ns,ms,sec,min,hour,temp;
	if(string==NULL){
		handle_error_with_exit("error in print time\n");
	}
	struct timespec time_current,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_current)!=0){
        	handle_error_with_exit("error in print_time\n");
    	}
	time_sub=diff_timespec(time_current,timer);
	if(clock_gettime(CLOCK_MONOTONIC,&timer)!=0){
        	handle_error_with_exit("error in clockgettime\n");
    	}
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	printf("%s:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",string,hour,min,sec,ms,ns);
	return;
}

void print_time_elapsed_local(char*string,struct timespec*timer_thread){
	//ns=1000000000=1 sec
	//ns=1000000 1 ms
	long ns,ms,sec,min,hour,temp;
	if(string==NULL || timer_thread==NULL){
		handle_error_with_exit("error in print time elapsed local\n");
	}
	struct timespec time_current,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_current)!=0){
		handle_error_with_exit("error in print_time\n");
	}
	time_sub=diff_timespec(time_current,*timer_thread);
	if(clock_gettime(CLOCK_MONOTONIC,timer_thread)!=0){
		handle_error_with_exit("error in clockgettime\n");
	}
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	printf("%s:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",string,hour,min,sec,ms,ns);
	return;
}
void print_time_elapsed_on_file_log(char*string){
	//ns=1000000000=1 sec
	//ns=1000000 1 ms
	long ns,ms,sec,min,hour,temp;
	if(string==NULL){
		handle_error_with_exit("error in print time\n");
	}
	struct timespec time_current,time_sub;
	if(clock_gettime(CLOCK_MONOTONIC,&time_current)!=0){
        	handle_error_with_exit("error in print_time\n");
    	}
	time_sub=diff_timespec(time_current,timer);
	ns=time_sub.tv_nsec%1000000;
	time_sub.tv_nsec-=ns;
	ms=time_sub.tv_nsec/1000000;
	sec=time_sub.tv_sec%60;//i secondi sono modulo 60
	min=(time_sub.tv_sec-sec)/60;
	temp=min%60;//temp min
	hour=(min-temp)/60;
	min=temp;
	fprintf(file_log,"%s:hour=%ld min=%ld sec:%ld ms=%ld ns:%ld",string,hour,min,sec,ms,ns);
	return;
}
void free_memory_matrix_mpz(mpz_t**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix mpz\n");
	}
	for(int i=0;i<num_row;i++){
		free_memory_array_mpz(matrix[i],num_col);
	}
	free(matrix);
	return;
}
void free_memory_array_mpz(mpz_t*array,long length){
	if(length<0){
		handle_error_with_exit("error in free memory array mpz\n");
	}
	if(length==0){
 		if(array==NULL){
			return;
		}
		else{
			handle_error_with_exit("error in parameter array free memory array mpz\n");
		}
	}
	for(long i=0;i<length;i++){
		mpz_clear(array[i]);
	}
	free(array);
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

void handle_error_with_exit(char*error_string){//uccide il processo dopo essersi accorto di un errore
    if(error_string==NULL){
	printf("error string is NULL\n");
	perror("");
        exit(EXIT_FAILURE);
    }
    printf("%s",error_string);
    //verifica che la stringa contenga il carattere newline
    fprintf(file_log,"%s",error_string);
    exit(EXIT_FAILURE);
}
