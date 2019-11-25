
#include "basic.h"
#include "main.h"
#include "matrix_function.h"
#include "list_square_relation.h"


extern struct timespec timer;
struct timespec time_start;
struct timespec timer_test;

extern long M;
extern int cardinality_factor_base;
extern long B;
extern double thresold_relation;

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





