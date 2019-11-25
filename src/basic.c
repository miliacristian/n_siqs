
#include "basic.h"
#include "main.h"
#include "matrix_function.h"
#include "list_square_relation.h"
#include "list_factorization.h"

extern struct timespec timer;
struct timespec time_start;
struct timespec timer_test;

extern long M;
extern int cardinality_factor_base;
extern long B;
extern double thresold_relation;




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
		thread_data[i].head_square=NULL;
		thread_data[i].tail_square=NULL;
        thread_data[i].head_residuos=NULL;
        thread_data[i].tail_residuos=NULL;
		free(thread_data[i].numbers);
		thread_data[i].numbers=NULL;
		free(thread_data[i].j1_mod_p);
		thread_data[i].j1_mod_p=NULL;
		free(thread_data[i].j2_mod_p);
		thread_data[i].j2_mod_p=NULL;

	}
	free(thread_data);
	return;
}
void clear_struct_thread_data(struct thread_data t_data,int M) {
	for (int i = 0; i < 2 * M + 1; i++) {
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
		t_data[i].head_square=NULL;
		t_data[i].tail_square=NULL;
		t_data[i].head_residuos=NULL;
		t_data[i].tail_residuos=NULL;
		t_data[i].numbers=malloc(sizeof(struct number)*(2*M+1));
		if(t_data[i].numbers==NULL){
			handle_error_with_exit("error in malloc alloc numbers");
		}
		t_data[i].j1_mod_p=alloc_array_long(cardinality_factor_base);
		t_data[i].j2_mod_p=alloc_array_long(cardinality_factor_base);
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




void print_estimated_time(int cardinality_factor_base,int num_B_smooth){
    if(cardinality_factor_base<=0 || num_B_smooth<0){
        handle_error_with_exit("error in print_estimated_time");
    }
    long ns,ms,sec,min,hour,temp;
    struct timespec empty_struct,time_sub;
    empty_struct.tv_nsec=0;
    empty_struct.tv_sec=0;
    time_sub=diff_timespec(timer_test,empty_struct);
    ns=time_sub.tv_nsec%1000000;
    time_sub.tv_nsec-=ns;
    ms=time_sub.tv_nsec/1000000;
    sec=time_sub.tv_sec%60;//i secondi sono modulo 60
    min=(time_sub.tv_sec-sec)/60;
    temp=min%60;//temp min
    hour=(min-temp)/60;
    min=temp;

    long times=((double)cardinality_factor_base*thresold_relation)/(double)num_B_smooth;
    ms=ms*times;
    long remainder=ms%1000;
    int plus_sec=(ms-remainder)/1000;//secondi in più
    ms=remainder;
    sec=sec*times;
    remainder=sec%60;
    int plus_min=(sec-remainder)/60;
    sec=remainder;
    min+=plus_min;
    sec+=plus_sec;
	remainder=sec%60;
	plus_min=(sec-remainder)/60;
	sec=remainder;
	min+=plus_min;
    printf("hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",hour,min,sec,ms,ns);
    return;
}





