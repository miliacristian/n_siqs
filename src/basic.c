#include "basic.h"

extern struct timespec timer;
struct timespec time_start;
struct timespec timer_test;

extern long M;
extern int cardinality_factor_base;
extern long B;
extern double thresold_relation;





void*thread_job_polynomial(void*arg){
	int*id=arg;
	thread_job_criv_quad(*id);
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




