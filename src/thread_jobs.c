#include "thread_jobs.h"
extern int cardinality_factor_base;
extern long B;
extern long M;
extern struct row_factorization r;
extern mpz_t a_old,a_new;//valore del coefficiente a del polinomio
extern mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
char divide_all_by_p_to_k_with_thread(int rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b,struct a_struct*array_a_struct){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p

    //j1=r-x0;
    //j2=-r-x0;
    if((p<=1 && p!=-1 ) || n==NULL || a==NULL || b==NULL || k<=0 || M<=0 || array_a_struct==NULL){
        handle_error_with_exit("invalid parameter divide all by p to k with thread\n");
    }
    char array_divided=0;//0 se nessuna divisione effettuata,1 se sono state effettutate divisioni per p^k
    long indexv,j1,j2,j_temp2;
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
    else if(!value_is_in_sorted_array(index_of_prime,array_a_struct,s)){//p non divide a
        if(r.inverse_a_mod_p[index_of_prime]==-1){
            print_array_a_struct(array_a_struct,s);
            printf("index_of_prime=%d\n",index_of_prime);
            handle_error_with_exit("error in factorize_matrix,inverse not found\n");
        }
        mpz_set_si(inverse_a,r.inverse_a_mod_p[index_of_prime]);////inverse_a=a^-1 mod p^k
        mpz_neg(j1t,b);//j1=-b
        mpz_neg(j2t,b);//j2=-b
        mpz_add_ui(j1t,j1t,rad);//j1=-b+r
        mpz_add(j2t,j2t,r2);//j2=-b+r2
        mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
        mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
    }
    else{//p divide a,esiste solo una soluzione,si usa solo j1,succede solamente poche volte(circa 10)
        //calcolo di c
        mpz_mul(v,b,b);//v=b^2
        mpz_sub(v,v,n);//v=b^2-n
        if(mpz_divisible_p(v,a)==0){//v non divide a
            handle_error_with_exit("error in mpz_divisible divide all by p_to_k\n");
        }
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
    mpz_mod(j1t,j1t,p_to_k);//j1t ridotto modulo p^k
    mpz_mod(j2t,j2t,p_to_k);//j1_t ridotto modulo p^k
    j1=mpz_get_si(j1t);//j1=j1t
    j2=mpz_get_si(j2t);//j2=j2t
    j_temp2=j1;
    indexv=j_temp2+M;//l'indice deve essere positivo
    while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        //indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
        thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
        //se il first_index_f_base non è stato modificato setta first e last
        if(thread_data.numbers[indexv].first_index_f_base==-1){
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
            //index minore di first-> diventa il nuovo first
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
        }
        else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
            //index maggiore di last-> diventa il nuovo last
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
            handle_error_with_exit("error in index\n");
        }
        array_divided=1;//una divisione è stata effettuata
        j_temp2+=p;//ad ogni ciclo aggiungo p
        indexv+=p;
    }
    j_temp2=-p+j1;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp2+M;
    while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        //indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
        thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
        //se il first_index_f_base non è stato modificato setta first e last
        if(thread_data.numbers[indexv].first_index_f_base==-1){
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
            //index minore di first-> diventa il nuovo first
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
        }
        else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
            //index maggiore di last-> diventa il nuovo last
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
            handle_error_with_exit("error in index\n");
        }
        array_divided=1;//una divisione è stata effettuata
        j_temp2-=p;//ad ogni ciclo tolgo p
        indexv-=p;
    }
    if(j==1){
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
        if(array_divided==1){
            return 1;//ritorna 1 se c'è stata almeno 1 divisione
        }
        return 0;//ritorna 0 se non ci sono state divisioni nell'array
    }
    j_temp2=j2;//si passa alla secodna radice
    indexv=j_temp2+M;
    while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        //indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
        thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
        //se il first_index_f_base non è stato modificato setta first e last
        if(thread_data.numbers[indexv].first_index_f_base==-1){
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
            //index minore di first-> diventa il nuovo first
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
        }
        else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
            //index maggiore di last-> diventa il nuovo last
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
            handle_error_with_exit("error in index\n");
        }
        array_divided=1;//una divisione è stata effettuata
        j_temp2+=p;//ad ogni ciclo aggiungo p
        indexv+=p;
    }
    j_temp2=-p+j2;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp2+M;
    while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
        //indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
        thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
        //se il first_index_f_base non è stato modificato setta first e last
        if(thread_data.numbers[indexv].first_index_f_base==-1){
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
            //index minore di first-> diventa il nuovo first
            thread_data.numbers[indexv].first_index_f_base=index_of_prime;
        }
        else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
            //index maggiore di last-> diventa il nuovo last
            thread_data.numbers[indexv].last_index_f_base=index_of_prime;
        }
        if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
            handle_error_with_exit("error in index\n");
        }
        array_divided=1;//una divisione è stata effettuata
        j_temp2-=p;//ad ogni ciclo tolgo p
        indexv-=p;
    }
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
    if(array_divided==1){
        return 1;//ritorna 1 se c'è stata almeno 1 divisione
    }
    return 0;//ritorna 0 se non ci sono state divisioni nell'array
}
char divide_all_by_2_log(long M,struct thread_data thread_data){//divide gli elementi dell'array di numeri per 2
    if(M<=0){
        handle_error_with_exit("invalid parameter divide_all_by_2_f\n");
    }
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
    for(;i<2*M+1;i=i+2){
        thread_data.numbers[i].last_index_f_base=1;
        thread_data.numbers[i].first_index_f_base=1;
        thread_data.numbers[i].sum_log=r.log_prime[1];
    }
    thread_data.j1_mod_p[1]=-1;
    thread_data.j2_mod_p[1]=-1;
    return 1;
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
void print_thread_data(struct thread_data thread_data,long M,int cardinality_factor_base){
    if(M<=0 || cardinality_factor_base<=0){
        handle_error_with_exit("error in print_thread_data\n");
    }
    if(M>THRESOLD_PRINT_ARRAY/2){
        return;
    }
    for(int i=0;i<2*M+1;i++){
        gmp_printf("b=%Zd,",thread_data.b);
        printf("first_index=%d,last_index=%d,sum_log=%d,j=%d\n",thread_data.numbers[i].first_index_f_base,thread_data.numbers[i].last_index_f_base,thread_data.numbers[i].sum_log,thread_data.numbers[i].j);
    }
    print_array_long(thread_data.j1_mod_p,cardinality_factor_base);
    print_array_long(thread_data.j2_mod_p,cardinality_factor_base);
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
int thread_job_to_create_factor_base(int id_thread){
    long remainder=reduce_int_mod_n_v2(B,NUM_THREAD_FACTOR_BASE+1);//rem=b mod num_thread
    long length=(B-remainder)/(NUM_THREAD_FACTOR_BASE+1);
    int start=id_thread*length+1;//se id=o start=1
    int end=start+length-1;
    //es remainder=0 thread=5 B=500.000 -> len=100.000 start=0*100000+1,end=1+100000-1=100000,start2=100001,end2=200000
    thread_factor_base_data[id_thread].last_prime_factor_base=start;
    create_factor_base_f(&(thread_factor_base_data[id_thread].cardinality_factor_base),end,&thread_factor_base_data[id_thread].head,&thread_factor_base_data[id_thread].tail,n,&(thread_factor_base_data[id_thread].last_prime_factor_base));
    return 0;
}
void*thread_job_factor_base(void*arg){
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
void*thread_factorization_job(void*arg){
    long p;
    struct factorization_thread_data*factorization_thread_data=arg;
    int start=(*factorization_thread_data).start;
    int end=(*factorization_thread_data).end;
    mpz_t a1;
    mpz_init(a1);
    if((*factorization_thread_data).id_thread==1){
    }
    if((*factorization_thread_data).is_a_default==1){
        mpz_set_si(a1,1);
    }
    else{
        mpz_set(a1,a_old);
    }
    for(int i=start;i<=end;i++){
        if(i==0){//skippa indice 0 e indice 1
            continue;
        }
        if(i==1){
            divide_all_by_2_log(M, (*factorization_thread_data).thread_data);
            continue;
        }
        p=r.prime[i];//primo iesimo della factor base
        divide_all_by_p_to_k_with_thread(r.root_n_mod_p[i],p,i,1,M,(*factorization_thread_data).thread_data,n,a1,(*factorization_thread_data).thread_data.b,array_a_struct);
    }
    mpz_clear(a1);
    return NULL;
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
