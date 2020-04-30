
#include "thread_jobs.h"

__thread unsigned int tid;
extern int cardinality_factor_base;
extern long B;
extern long M;
extern long num_elem_array_number;
extern struct row_factorization r;
extern mpz_t a_old,a_new;//valore del coefficiente a del polinomio
extern mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
extern int num_thread_job;
extern mpz_t *array_bi;
extern mpz_t thresold_large_prime;
extern struct thread_data*thread_polynomial_data;


void find_list_square_relation(struct thread_data *pthread_data, int *num_B_smooth,int*num_semi_B_smooth, int *num_potential_B_smooth, long M,
                               struct node_square_relation **head_square, struct node_square_relation **tail_square,struct node_square_relation **head_residuos, struct node_square_relation **tail_residuos,
                               const mpz_t n,const mpz_t a,const struct a_struct*array_a_struct,int s) {
    mpz_t num;
    struct thread_data thread_data=*pthread_data;
    struct node_factorization*head_factor=NULL;
    char is_B_smooth=0;
    char is_semi_B_smooth=0;
    mpz_t residuos;
    #if DEBUG==1
    if(num_B_smooth==NULL || num_semi_B_smooth==NULL || num_potential_B_smooth==NULL || M<=0 ||
       head_square==NULL || tail_square==NULL || head_residuos==NULL
       || tail_residuos==NULL || (array_a_struct==NULL && s>0)){
        handle_error_with_exit("error in find_list_square_relation");
    }
    #endif
    struct square_relation square_relation;
    mpz_init(num);
    mpz_init(residuos);
    const float log_thresold=thread_data.log_thresold;
    mpz_t b;
    mpz_init(b);
    mpz_set(b,thread_data.b);
    long i_minus_M=-M;
    const long first_index_f_base=1;
    const long last_index_f_base=cardinality_factor_base-1;
    int*numbers_i_dot_sum_log=&(thread_data.numbers[0].sum_log);//===thread_data.numbers[i].sum_log
    int *numbers_i_dot_j=&(thread_data.numbers[0].j);//===thread_data.numbers[i].j
    for(long i=0;i<num_elem_array_number;i++,i_minus_M+=1,numbers_i_dot_sum_log+=NUM_INT_IN_STRUCT_NUMBER,numbers_i_dot_j+=NUM_INT_IN_STRUCT_NUMBER){

        if(*numbers_i_dot_sum_log>=log_thresold){
            //possibile B_smooth trovato
            (*num_potential_B_smooth)++;
            create_num(num,a,b,n,*numbers_i_dot_j);
            if(mpz_cmp_si(num,0)==0){
                head_factor=NULL;
                mpz_set_si(residuos,0);
                is_B_smooth=0;
                is_semi_B_smooth=0;
                (*num_potential_B_smooth)--;
                continue;
            }
            else {
                head_factor = factorize_num(num, *numbers_i_dot_j,first_index_f_base,
                                               last_index_f_base, &is_B_smooth,
                                               &is_semi_B_smooth, residuos, array_a_struct, s, &thread_data);
            }
            #if DEBUG==1
            if(head_factor==NULL && (is_B_smooth==1 || is_semi_B_smooth==1)){
                handle_error_with_exit("error invalid factorize_num\n");
            }
            if(verify_factorization(num,residuos,head_factor,a)==0){
                handle_error_with_exit("error in factorization\n");
            }
            #endif
            if(head_factor==NULL){
                continue;
            }
            if(is_B_smooth==1){
                (*num_B_smooth)++;
                square_relation.head_factorization=head_factor;
                mpz_init(square_relation.square);
                mpz_init(square_relation.num);
                mpz_init(square_relation.residuos);
                mpz_set(square_relation.num,num);
                mpz_set(square_relation.residuos,residuos);
                calculate_square(square_relation.square,a,i_minus_M,b,n);
                insert_at_tail_square_relation(square_relation,head_square,tail_square);
                #if DEBUG==1
                if(verify_square_relation(square_relation,n)==0){
                    handle_error_with_exit("error in create square relation\n");
                }
                #endif
            }
            else if(is_semi_B_smooth==1){
                if(mpz_cmp(residuos,thresold_large_prime)<=0) {
                    (*num_semi_B_smooth)++;
                    square_relation.head_factorization = head_factor;
                    mpz_init(square_relation.square);
                    mpz_init(square_relation.num);
                    mpz_init(square_relation.residuos);
                    mpz_set(square_relation.num, num);
                    mpz_set(square_relation.residuos, residuos);
                    calculate_square(square_relation.square, a, i_minus_M,b, n);
                    insert_at_tail_square_relation(square_relation, head_residuos, tail_residuos);
                }
                else{
                    free_memory_list_factor(head_factor);
                }
            }
            else{
                handle_error_with_exit("error in find list square relation\n");
            }
        }
    }
    mpz_clear(num);
    mpz_clear(residuos);
    mpz_clear(b);
    return;
}
void divide_all_by_p_to_k(int rad,long p,int index_of_prime,long k,long M,struct thread_data *pthread_data,const mpz_t n,const mpz_t a,const mpz_t b,const struct a_struct*array_a_struct,int*index_array_a_struct,int s){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p

    //j1=r-x0;
    //j2=-r-x0;
    struct thread_data thread_data=*pthread_data;
    #if DEBUG==1
    if((p<=1 && p!=-1 ) || n==NULL || a==NULL || b==NULL || k<=0 || M<=0 || (array_a_struct==NULL && s>0) || index_array_a_struct==NULL || s<0){
        handle_error_with_exit("invalid parameter divide all by p to k\n");
    }
    #endif
    struct timespec timer_factor;
    gettime(&timer_factor);
    long indexv,j1,j2,j_temp_long;
    mpz_t p_to_k,r2,a_temp,inverse_a,inverse2b,j1t,j2t,c,b_temp,v,l_temp,j_temp,index,h_temp;
    mpz_init(p_to_k);
    mpz_init(r2);
    mpz_init(a_temp);
    mpz_init(inverse_a);
    mpz_init(inverse2b);
    mpz_init(j1t);
    mpz_init(j2t);
    mpz_init(c);
    mpz_init(b_temp);
    mpz_init(v);
    mpz_init(l_temp);
    mpz_init(j_temp);
    mpz_init(h_temp);
    mpz_init(index);
    mpz_set_si(p_to_k,p);//p^k=p
    mpz_set_si(r2,r.root2_n_mod_p[index_of_prime]);//r2=p^k-r seconda radice quadrata di n modulo p^k
    long j=0;//indici dell'array divisibile per p^k,se j!=0 j2 non esiste
    if(mpz_cmp_si(a,1)==0){//se a=1,infatti a^-1 mod p =1
        mpz_neg(j1t,b);//j1=-b
        mpz_neg(j2t,b);//j2=-b
        mpz_add_ui(j1t,j1t,rad);//j1=-b+r
        mpz_add(j2t,j2t,r2);//j2=-b+r2
    }

    else if((*index_array_a_struct<s && array_a_struct[*index_array_a_struct].index_prime_a!=index_of_prime) || *index_array_a_struct>=s){//p non divide a
        #if DEBUG==1
        if(r.inverse_a_mod_p[index_of_prime]==-1){
            printf("p=%ld,index_of_prime=%d,index_array_a_struct=%d\n",p,index_of_prime,array_a_struct[*index_array_a_struct].index_prime_a);
            handle_error_with_exit("error in factorize_matrix,inverse not found\n");
        }
        #endif
        mpz_set_si(inverse_a,r.inverse_a_mod_p[index_of_prime]);////inverse_a=a^-1 mod p^k
        mpz_neg(j1t,b);//j1=-b
        mpz_neg(j2t,b);//j2=-b
        mpz_add_ui(j1t,j1t,rad);//j1=-b+r
        mpz_add(j2t,j2t,r2);//j2=-b+r2
        mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
        mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
    }
    else if(array_a_struct[*index_array_a_struct].index_prime_a==index_of_prime){//p divide a,esiste solo una soluzione,si usa solo j1,succede solamente poche volte(circa 10)
        //calcolo di c
        #if DEBUG==1
        if(mpz_divisible_ui_p(a,p)==0){
            handle_error_with_exit("errore p non divide a\n");
        }
        #endif
        (*index_array_a_struct)++;
        mpz_mul(v,b,b);//v=b^2
        mpz_sub(v,v,n);//v=b^2-n
        #if DEBUG==1
        if(mpz_divisible_p(v,a)==0){//v non è divisibile per a
            handle_error_with_exit("error in mpz_divisible divide all by p_to_k\n");
        }
        #endif
        mpz_divexact(c,v,a);//c=(b^2-n)/a

        if(k==1){//p_to_k=p,esiste solo una soluzione mod p
            mpz_set(b_temp,b);
            mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
            mpz_invert(inverse2b,b_temp,p_to_k);//(2b)^-1 mod p^k
            mpz_neg(j1t,c);//j1=-c
            mpz_mul(j1t,j1t,inverse2b);//j1=-c*inverse2b;
            j=1;
        }
        else{//k>1,p_to_k>p,gcd(a,p^k)=p,gcd(a/p,p^k)=1
            handle_error_with_exit("error in divide all by p_to_k k>1\n");
            calculate_root_poly_second_degree_mod_p_to_k(j1t,p_to_k,p,k,a,b,n);
            j=1;//solo 1 soluzione
        }
    }
    else{
        handle_error_with_exit("caso non gestito in divide all by p_to_k\n");
    }
    mpz_mod(j1t,j1t,p_to_k);//j1t ridotto modulo p^k
    mpz_mod(j2t,j2t,p_to_k);//j1_t ridotto modulo p^k
    j1=mpz_get_si(j1t);//j1=j1t
    j2=mpz_get_si(j2t);//j2=j2t
    thread_data.j1_mod_p[index_of_prime]=j1;
    if(j!=1) {
        thread_data.j2_mod_p[index_of_prime] = j2;
    }
    else{//j==1
        thread_data.j2_mod_p[index_of_prime]= NO_INDEX;
    }
    //fine calcolo delle radici,ora ciclare sull'array e sommmare i log
    //print_time_elapsed_local("time to calculate root",&timer_factor,tid);
    const int log_p=r.log_prime[index_of_prime];
    j_temp_long=j1;
    indexv=j_temp_long+M;
    const long double_M=2*M;
    //&& 0 perché la somma dei logaritmi viene fatta successivamente
    while(indexv<=double_M && 0){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        thread_data.numbers[indexv].sum_log+=log_p;
        indexv+=p;
    }
    j_temp_long=-p+j1;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp_long+M;
    //&& 0 perché la somma dei logaritmi viene fatta successivamente
    while(indexv>=0 && 0){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        thread_data.numbers[indexv].sum_log+=log_p;
        indexv-=p;
    }
    if(j==1){
        //print_time_elapsed_local("time to sieve array numbers............",&timer_factor,tid);
        mpz_clear(p_to_k);
        mpz_clear(r2);
        mpz_clear(a_temp);
        mpz_clear(inverse_a);
        mpz_clear(inverse2b);
        mpz_clear(j1t);
        mpz_clear(j2t);
        mpz_clear(c);
        mpz_clear(b_temp);
        mpz_clear(v);
        mpz_clear(l_temp);
        mpz_clear(j_temp);
        mpz_clear(h_temp);
        mpz_clear(index);
        return ;//ritorna 0 se non ci sono state divisioni nell'array
    }
    j_temp_long=j2;//si passa alla seconda radice
    indexv=j_temp_long+M;
    //&& 0 perché la somma dei logaritmi viene fatta successivamente
    while(indexv<=double_M && 0){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        thread_data.numbers[indexv].sum_log+=log_p;
        indexv+=p;
    }
    j_temp_long=-p+j2;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp_long+M;
    //&& 0 perché la somma dei logaritmi viene fatta successivamente
    while(indexv>=0 && 0){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        thread_data.numbers[indexv].sum_log+=log_p;
        indexv-=p;
    }
    //print_time_elapsed_local("time to sieve array numbers...............",&timer_factor,tid);
    mpz_clear(p_to_k);
    mpz_clear(r2);
    mpz_clear(a_temp);
    mpz_clear(inverse_a);
    mpz_clear(inverse2b);
    mpz_clear(j1t);
    mpz_clear(j2t);
    mpz_clear(c);
    mpz_clear(b_temp);
    mpz_clear(v);
    mpz_clear(l_temp);
    mpz_clear(j_temp);
    mpz_clear(h_temp);
    mpz_clear(index);
    return ;//ritorna 0 se non ci sono state divisioni nell'array
}


void divide_all_by_2_log(long M,struct thread_data *pthread_data){//divide gli elementi dell'array di numeri per 2
    #if DEBUG==1
    if(M<=0){
        handle_error_with_exit("invalid parameter divide_all_by_2_f\n");
    }
    #endif
    struct thread_data thread_data=*pthread_data;
    long i=-1;
    long count=0;
    count=M & 1;//se M dispari count =1 altrimenti count=0
    count+=1;//a^2*j^2+2bj+c mod 2 =M+1+b,primo elemento ha j=-M=M mod 2
    //trovare il primo numero divisibile per 2 e poi a salti di 2 dividere gli altri
    if(mpz_divisible_2exp_p(thread_data.b,1)!=0){//b è divisibile per 2
    }
    else{
        count+=1;
    }
    //i=1 indice dispari,il primo elemento è dispari,il secondo è pari,i=0 indice pari,il primo elemento è pari
    if(count%2==0){
        i=0;//primo elemento divisibile per 2
    }
    else{
        i=1;//secondo elemento divisibile per 2
    }
    thread_data.j1_mod_p[1]=i;
    thread_data.j2_mod_p[1]=NO_INDEX;
    return;
}
/*void factor_matrix(const mpz_t n,long M,struct thread_data *thread_data,const int cardinality_factor_base,const mpz_t a,
                     const struct a_struct*array_a_struct,int s){
    #if DEBUG==1
    if(n==NULL || a==NULL || mpz_sgn(n)<=0 || M<=0 || cardinality_factor_base<=0 || (array_a_struct==NULL && s>0) || s<0){
        handle_error_with_exit("error in factor matrix_f\n");
    }
    #endif
    long p;
    int index_array_a_struct=0;
    for (int i=0;i<cardinality_factor_base;i++){//per ogni elemento della factor base
        p=r.prime[i];//primo iesimo della factor base
        if(p==-1){
            continue;
        }
        if(p==2){
            divide_all_by_2_log(M,thread_data);
        }
        else{//p>2 e dispari
            divide_all_by_p_to_k_f(r.root_n_mod_p[i],p,i,1,M,thread_data,n,a,thread_data->b,array_a_struct,&index_array_a_struct,s);
        }
    }
    return;
}*/

void factor_matrix(const mpz_t n,long M,struct thread_data *thread_data,const int cardinality_factor_base,const mpz_t a,
                     const struct a_struct*array_a_struct,int s){
    #if DEBUG==1
    if(n==NULL || a==NULL || mpz_sgn(n)<=0 || M<=0 || cardinality_factor_base<=0 || (array_a_struct==NULL && s>0) || s<0){
        handle_error_with_exit("error in factor matrix_f\n");
    }
    #endif
    struct timespec timer_sieving,timer_roots;
    gettime(&timer_roots);
    long p;
    int index_array_a_struct=0;
    long prime;
    int next_index_rad1,next_index_rad2;
    for (int i=0;i<cardinality_factor_base;i++){//per ogni elemento della factor base
        p=r.prime[i];//primo iesimo della factor base
        if(p==-1){
            continue;
        }
        if(p==2){
            divide_all_by_2_log(M,thread_data);
        }
        else{//p>2 e dispari
            divide_all_by_p_to_k(r.root_n_mod_p[i],p,i,1,M,thread_data,n,a,thread_data->b,array_a_struct,&index_array_a_struct,s);
        }
    }
    print_time_elapsed_local("time to calculate roots...............",&timer_roots,tid);
    gettime(&timer_sieving);
    for(int index_prime=1;index_prime<cardinality_factor_base;index_prime++){
        prime=r.prime[index_prime];//primo iesimo della factor base
        next_index_rad1=thread_data->j1_mod_p[index_prime]+M;
        next_index_rad1%=prime;
        if(thread_data->j2_mod_p[index_prime]!=NO_INDEX){
            next_index_rad2=thread_data->j2_mod_p[index_prime]+M;
            next_index_rad2%=prime;
        }
        else{//esiste una sola soluzione
            next_index_rad2=NO_INDEX;
        }
        int log_prime=r.log_prime[index_prime];
        for(;next_index_rad1<2*M+1;next_index_rad1+=prime){
            thread_data->numbers[next_index_rad1].sum_log+=log_prime;
        }
        if(next_index_rad2!=NO_INDEX){
            for(;next_index_rad2< 2*M+1;next_index_rad2+=prime){
                thread_data->numbers[next_index_rad2].sum_log+=log_prime;
            }
        }
    }
    print_time_elapsed_local("time to sieving...............",&timer_sieving,tid);
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
void print_thread_data(struct thread_data *pthread_data,long M,int cardinality_factor_base){
    #if DEBUG==1
    if(M<=0 || cardinality_factor_base<=0){
        handle_error_with_exit("error in print_thread_data\n");
    }
    #endif
    struct thread_data thread_data=*pthread_data;
    if(M>THRESOLD_PRINT_ARRAY/2){
        return;
    }
    for(int i=0;i<2*M+1;i++){
        gmp_printf("b=%Zd,",thread_data.b);
        printf("sum_log=%d,j=%d\n",thread_data.numbers[i].sum_log,thread_data.numbers[i].j);
    }
    print_array_long(thread_data.j1_mod_p,cardinality_factor_base,THRESOLD_PRINT_ARRAY);
    print_array_long(thread_data.j2_mod_p,cardinality_factor_base,THRESOLD_PRINT_ARRAY);
}

void clear_struct_thread_data(struct thread_data *pt_data) {
    struct number*numbers=pt_data->numbers;
    const long end= num_elem_array_number;
    int*sum_log_address=(int*)((unsigned long)(numbers)+OFFSET_FROM_STRUCT_NUMBER_TO_SUM_LOG);
    for (int i = 0; i < end; i++) {
        *sum_log_address=0;
        sum_log_address+=NUM_INT_IN_STRUCT_NUMBER;
    }
    return;
}
struct thread_data* alloc_array_polynomial_thread_data(int length_array_thread_data,long M){
    #if DEBUG==1
    if(length_array_thread_data<=0 || M<=0){
        handle_error_with_exit("error in alloc_array_thread_data\n");
    }
    if((sizeof(struct thread_data)& 63)!=0){
        handle_error_with_exit("struct thread data is not aligned to cache line\n");
    }
    #endif
    struct thread_data*t_data=malloc(sizeof(struct thread_data)*length_array_thread_data);
    if(t_data==NULL){
        handle_error_with_exit("error in malloc alloc array thread_data\n");
    }
    float log_thresold=calculate_log_thresold(n,M);
    for(int i=0;i<length_array_thread_data;i++){//for each thread
        mpz_init(t_data[i].b);
        t_data[i].log_thresold=log_thresold;
        t_data[i].num_potential_B_smooth=0;
        t_data[i].num_B_smooth=0;
        t_data[i].num_semi_B_smooth=0;
        t_data[i].head_square=NULL;
        t_data[i].tail_square=NULL;
        t_data[i].head_residuos=NULL;
        t_data[i].tail_residuos=NULL;
        t_data[i].numbers=malloc(sizeof(struct number)*(2*M+1));//array dei 2M+1 numeri
        if(t_data[i].numbers==NULL){
            handle_error_with_exit("error in malloc alloc numbers");
        }
        t_data[i].j1_mod_p=alloc_array_long(cardinality_factor_base);
        t_data[i].j2_mod_p=alloc_array_long(cardinality_factor_base);
    }
    for(int i=0;i<2*M+1;i++){
        t_data[0].numbers[i].j=i-M;
        t_data[0].numbers[i].sum_log=0;
    }
    for(int i=0;i<length_array_thread_data;i++){//for each thread
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
    long remainder=reduce_int_mod_n_v2(B,NUM_THREAD_FACTOR_BASE+1);//rem=b mod num_thread
    long length=(B-remainder)/(NUM_THREAD_FACTOR_BASE+1);
    int start=id_thread*length+1;//se id=0 start=1
    return start;
}
int thread_job_to_create_factor_base(int id_thread){
    long remainder=reduce_int_mod_n_v2(B,NUM_THREAD_FACTOR_BASE+1);//rem=b mod num_thread
    long length=(B-remainder)/(NUM_THREAD_FACTOR_BASE+1);
    int start=id_thread*length+1;//se id=0 start=1
    int end=start+length-1;
    //es remainder=0 thread=5 B=500.000 -> len=100.000 start=0*100000+1,end=1+100000-1=100000,start2=100001,end2=200000
    thread_factor_base_data[id_thread].last_prime_factor_base=start;
    create_factor_base_f(&(thread_factor_base_data[id_thread].cardinality_factor_base),end,&thread_factor_base_data[id_thread].head,&thread_factor_base_data[id_thread].tail,n,&(thread_factor_base_data[id_thread].last_prime_factor_base));
    return 0;
}
void* thread_job_factor_base(void*arg){
    int*id=arg;
    thread_job_to_create_factor_base(*id);
    pthread_exit(NULL);
    return NULL;
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

int thread_job_criv_quad(int id_thread){//id inizia da 0,il lavoro di un thread rimane uguale anche se non si riesce a fattorizzare n
    tid=id_thread;
    int cpu = tid;
    cpu_set_t oldset;
    if (pin_thread_to_core(cpu,&oldset))
    {
        handle_error_with_exit("impossible pinning thread to core\n");
    }
    int local_num_thread_job=num_thread_job;
    const int end_condition=local_num_thread_job-2;
    if(id_thread>end_condition){//l'indice del thread eccede il numero di job da fare
        return 0;
    }
    struct timespec timer_thread;//istante di tempo
    struct timespec timer_thread_start;
    struct node_square_relation*head_squares=NULL,*tail_squares=NULL;
    struct node_square_relation*head_residuoss=NULL,*tail_residuoss=NULL;

    //gettime
    gettime(&timer_thread);
    gettime(&timer_thread_start);
    struct thread_data *local_thread_data=&(thread_polynomial_data[id_thread]);
    mpz_t *local_b=&(array_bi[id_thread]);
    
    for(int count=id_thread;count<=end_condition;count+=NUM_THREAD_POLYNOMIAL,local_b+=NUM_THREAD_POLYNOMIAL){//ogni thread prende un sottoinsieme di compiti,il thread con id 0 farà i compiti 0,NUM_THREAD,2*NUM_THREAD,il thread 1 farà 1,NUM_THREAD+1,2*NUM_THREAD+1 ecc
        //fattorizzazione,alla fine ogni thread ha una lista di relazioni quadratiche

        //fattorizza array di 2m+1 elementi e memorizza la somma dei logaritmi per ogni posizione e
        // indici last e first che ci dicono il primo elemento divisibile per num e l'ultimo(questo facilita la trial division)
        //mpz_set(local_thread_data->b,array_bi[count]);//imposta ad ogni ciclo il valore di b,array_bi è globale ma non viene più acceduto,thread_polynomial_data[id_thread].b=array_bi[count] è globale ma viene acceduto solo in lettura
        mpz_set(local_thread_data->b,*local_b);
        factor_matrix(n,M,(local_thread_data),cardinality_factor_base,a_old,array_a_struct,s);//fattorizza una nuova matrice,tutte variabili locali o globali (ma accedute in sola lettura)
        //print_time_elapsed_local("time to factor_matrix",&timer_thread,tid);
        //ricerca dei B_smooth potenziali,reali e fattorizzazione dei B_smooth reali
        find_list_square_relation((local_thread_data),&(local_thread_data->num_B_smooth),&(local_thread_data->num_semi_B_smooth),&(local_thread_data->num_potential_B_smooth),M,&head_squares,&tail_squares,&head_residuoss,&tail_residuoss,n,a_old,array_a_struct,s);
        //pulisci struttura dati del thread per ricominciare con un altro polinomio
        clear_struct_thread_data((local_thread_data));
        //print_time_elapsed_local("time to find_list_square_relationnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn",&timer_thread,tid);
        //unisci la lista dei quadrati trovata con il polinomio con la lista dei quadrati del thread,alla fine ogni thread ha un unica lista dei quadrati
        union_list_square(&(local_thread_data->head_square),&(local_thread_data->tail_square),head_squares,tail_squares);
        union_list_residuos(&(local_thread_data->head_residuos),&(local_thread_data->tail_residuos),head_residuoss,tail_residuoss);
        //print_time_elapsed_local("time to union list",&timer_thread,tid);
        head_squares=NULL;//resetta la lista locale delle relazioni quadratiche
        tail_squares=NULL;//resetta la lista locale delle relazioni quadratiche
        head_residuoss=NULL;
        tail_residuoss=NULL;
    }
    print_time_elapsed_local("time to finish thread job",&timer_thread_start,tid);
    return 0;
}

struct node_factorization* factorize_num(const mpz_t num,int j_of_num,int first_index_f_base,int last_index_f_base,
                                           char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,const struct a_struct*array_a_struct,int s,struct thread_data *pthread_data){
    #if DEBUG==1
    if((first_index_f_base<0)|| (last_index_f_base<0 && last_index_f_base!=-1)
       || is_B_smooth==NULL || is_semi_B_smooth==NULL
       || first_index_f_base>last_index_f_base || (array_a_struct==NULL && s>0)){
        handle_error_with_exit("error in factorize_num\n");
    }
    #endif
    //j_of_num va da -M a +M
    struct thread_data thread_data=*pthread_data;
    int prime=0;
    int exp=0;
    int index=0;
    long j1,j2,j_of_num_mod_p;
    mpz_t temp;
    mpz_init(temp);
    mpz_set(temp,num);
    struct node_factorization*head=NULL;
    struct node_factorization*tail=NULL;
    //bisogna aggiungere i fattori di a alla fattorizzazione

    if(mpz_cmp_si(num,0)<0){//valore negativo->divisibile per 1
        mpz_neg(temp,temp);//rendilo positivo
        alloc_and_insert_at_tail(-1,1,0,&head,&tail);//inserisci nodo -1
    }
    if(s!=0) {
        //se i numeri di a sono ad un indice più piccolo del first_index_f_base allora aggiungili prima
        //questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
        while (index!=s && array_a_struct[index].index_prime_a < first_index_f_base) {
            alloc_and_insert_at_tail(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
                                  &tail);//inserisci nodo con esponente 1
            #if DEBUG==1
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
            #endif
            index++;
        }
    }
    exp=0;
    int*prime_i=&(r.prime[first_index_f_base]);
    long*j1_i=&(thread_data.j1_mod_p[first_index_f_base]);
    long*j2_i=&(thread_data.j2_mod_p[first_index_f_base]);
    for(int i=first_index_f_base;i<=last_index_f_base;i++,j1_i++,j2_i++,prime_i++){//scorri i primi da fist index a last index compresi
        #if DEBUG==1
        if(exp!=0){
            handle_error_with_exit("error in factorize num v2,invalid exponent\n");
        }
        #endif
        if(s!=0 && index!=s && array_a_struct[index].index_prime_a==i){//se un fattore di a ha un indice compreso
            // tra first e last index metti ad uno l'esponente anche se non è divisibile per quel numero
            // (sappiamo che il polinomio è divisibile per quel numero)
            #if DEBUG==1
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
            #endif
            exp+=1;//poni esponente a 1
            index++;//aumenta l'indice
            //n.b. i fattori di a possono comparire nella fattorizzazione con un esponente maggiore di 1
        }
        //prime=r.prime[i];
        prime=*prime_i;
        if(i==1 && (mpz_divisible_2exp_p(temp,1)!=0) ){//divisibile per 2
                mpz_divexact_ui(temp, temp,prime);//dividilo la prima volta per prime
                exp+=1;
                while(mpz_divisible_2exp_p(temp,1)!=0){//finquando è divisibie per 2
                    mpz_divexact_ui(temp, temp,prime);//dividilo
                    exp+=1;//aumenta esponente
                }
                if(exp>0) {//è stato diviso almeno 1 volta
                    alloc_and_insert_at_tail(prime, exp, i, &head, &tail);
                }
                exp=0;//resetta esponente
                if(mpz_cmp_si(temp,1)==0){//se il residuo della divisione è 1 allora è B-smooth
                    last_index_f_base=i;//goto add exponent of a and exit
                    break;
                }
                continue;
        }
        j1=*j1_i;
        j2=*j2_i;
        //j1=thread_data.j1_mod_p[i];
        //j2=thread_data.j2_mod_p[i];
        j_of_num_mod_p=reduce_int_mod_n_v2(j_of_num,prime);
        //TODO invece di calcolare il modulo rispetto all'indice attuale è 
        //possibile calcolarlo rispetto al modulo dell'indice precedente con l'operazione mod(mod(indice precedente)+(indice attuale-indice precedente)) )
        //esempio se indice precendente è 9 e modulo 5 veniva 4 e ora indice è 15 modulo 5 viene 0 ==mod5(4+(15-9)) In questo modo i numeri sono tutti più trattabili
        //per implementarlo bisogna avere un array tampone che memorizza il modulo precedente per ogni prime della factor base(o per entrambe le radici di un primo della factor base,vedere come è più comodo)
        //j_of_num ridotto modulo prime se è uguale a j1 o j2 allora è divisibile per p
        if(j_of_num_mod_p==j1 || j_of_num_mod_p==j2){
            mpz_divexact_ui(temp, temp,prime);//dividilo la prima volta per prime
            exp+=1;//aumenta esponente
            //verifica se è ulteriormente divisibile per prime facendo la divisione
            while(mpz_divisible_ui_p(temp,prime)!=0) {//se il numero è divisibile per un primo della fattor base
                mpz_divexact_ui(temp, temp,prime);//dividilo
                exp+=1;//aumenta esponente
            }
            if(exp>0) {//è stato diviso almeno 1 volta
                alloc_and_insert_at_tail(prime, exp, i, &head, &tail);
            }
            exp=0;//resetta esponente
            if(mpz_cmp_si(temp,1)==0){//se il residuo della divisione è 1 allora è B-smooth
                last_index_f_base=i;//goto add exponent of a and exit
                break;
            }
            continue;
        }
        if(exp>0) {//è un fattore di a
            alloc_and_insert_at_tail(prime, exp, i, &head, &tail);
        }
        exp=0;//resetta esponente
        continue;
    }
    if(s!=0) {
        //se i numeri di a sono ad un indice più grande del last_index_f_base allora aggiungili dopo
        //questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
        while (index!=s && array_a_struct[index].index_prime_a > last_index_f_base && index < s) {
            alloc_and_insert_at_tail(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
                                  &tail);//inserisci nodo
            #if DEBUG==1
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
            #endif
            index++;
        }
    }
    #if DEBUG==1
    if(s!=0 && index!=s){
        handle_error_with_exit("error in factorize num invalid s\n");
    }
    #endif
    if(head==NULL){
        *is_B_smooth=0;
        *is_semi_B_smooth=0;
        mpz_set(residuos, temp);//imposta il residuo
        mpz_clear(temp);
        return NULL;
    }
    if(mpz_cmp_si(temp,1)==0){//se il residuo della divisione è 1 allora è B-smooth
        *is_B_smooth=1;
        *is_semi_B_smooth=0;
        mpz_set(residuos,temp);//imposta il residuo
        mpz_clear(temp);
        return head;
    }
    else{//non è B-smooth è semi_B_smooth
        *is_B_smooth=0;
        *is_semi_B_smooth=1;
        mpz_set(residuos,temp);//imposta il residuo
        mpz_clear(temp);
        return head;
    }
    handle_error_with_exit("never reached\n");
    return NULL;
}

void* thread_job_polynomial(void*arg){
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