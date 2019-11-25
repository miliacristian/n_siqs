#include "basicv2.h"
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
