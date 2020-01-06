#define _GNU_SOURCE
#include "basic.h"
int pin_thread_to_core(int core,cpu_set_t*oldset)
{
    int sched_cpu;
    cpu_set_t cpuset;

    pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), oldset);

    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);

    if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset))
        return 1;

    while (1)
    {
        if ((sched_cpu = sched_getcpu()) == -1)
            return 1;
        else if (sched_cpu == core)
            break;
    }

    return 0;
}
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

FILE*open_file_log(){
	FILE*file=fopen("file_log.txt","a");//modalità append
	if(file==NULL){
		handle_error_with_exit("error in fopen,impossible open file\n");
	}
	return file;
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
void handle_error_with_exit(char*error_string){//uccide il processo dopo essersi accorto di un errore
    if(error_string==NULL){
        printf("error string is NULL\n");
        perror("");
        exit(EXIT_FAILURE);
    }
    printf("%s",error_string);
    exit(EXIT_FAILURE);
}