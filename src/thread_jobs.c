#include "thread_jobs.h"
extern int cardinality_factor_base;
extern long B;
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
int calculate_start_factor_base(int id_thread){
    long remainder=reduce_int_mod_n_v2(B,NUM_THREAD_FACTOR_BASE+1);
    long length=(B-remainder)/(NUM_THREAD_FACTOR_BASE+1);
    int start=id_thread*length+1;
    return start;
}