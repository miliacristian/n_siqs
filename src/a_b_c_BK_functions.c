
#include "a_b_c_BK_functions.h"

void print_array_Bk(mpz_t*array_Bk,long s){
    if(s<0){
        handle_error_with_exit("error in print_array_Bk\n");
    }
    if(array_Bk==NULL){
        printf("array_Bk is empty\n");
        return;
    }
    if(not_print_array(s)==1){
        return;
    }
    printf("array Bk:\n");
    print_array_mpz(array_Bk,s);
    return;
}
void print_array_bi(mpz_t* array_bi,long s){
    if(s<0){
        handle_error_with_exit("error in print_array_bi\n");
    }
    if(array_bi==NULL){
        printf("array_bi is empty\n");
        return;
    }
    long length=(long)pow(2,s-1);
    if(not_print_array(length)==1){
        return;
    }
    printf("array_bi:\n");
    print_array_mpz(array_bi,length);
    return;

}
void print_array_chosen_for_a(long*array_of_prime_chosen_for_a,int card_factor_base,int s){
    if(card_factor_base<=0){
        handle_error_with_exit("error in card factor base print_array_chosen_for_a\n");
    }
    if(array_of_prime_chosen_for_a==NULL){
        if(s==0){
            return;
        }
        else{
            handle_error_with_exit("error in print_array_chosen_for_a\n");
        }
    }
    printf("array of prime chosen for a:\n");
    print_array_long(array_of_prime_chosen_for_a,card_factor_base*2);
    return;
}
void create_num(mpz_t num,const mpz_t a,const mpz_t b,const mpz_t n,long j){
    //num=a*j^2+2*bj+c mod n
    mpz_t c,double_b_mul_j,a_mul_square_j;

    mpz_init(c);
    mpz_init(double_b_mul_j);
    mpz_init(a_mul_square_j);

    //a_mul_square_j=a
    mpz_set(a_mul_square_j,a);//a_mul_suqare_j=a
    mpz_mul_si(a_mul_square_j,a_mul_square_j,j);//a_mul_suqare_j=a*j
    mpz_mul_si(a_mul_square_j,a_mul_square_j,j);//a_mul_suqare_j=a*j^2

    //double_b_mul_j
    mpz_set(double_b_mul_j,b);//double_b_mul_j=b
    mpz_mul_2exp(double_b_mul_j,double_b_mul_j,1);//double_b_mul_j=2*b
    mpz_mul_si(double_b_mul_j,double_b_mul_j,j);//double_b_mul_j=2*b*j

    //c
    mpz_mul(c,b,b);//c=b^2
    mpz_sub(c,c,n);//c=b^2-n
    if(mpz_divisible_p(c,a)==0){
        handle_error_with_exit("error in create_num\n");
    }
    mpz_divexact(c,c,a);//c=b^2-n/a

    mpz_add(num,a_mul_square_j,double_b_mul_j);
    mpz_add(num,num,c);//num=a*j^2+2*bj+c

    mpz_clear(c);
    mpz_clear(double_b_mul_j);
    mpz_clear(a_mul_square_j);
}
void square_root_mod_p_to_k(mpz_t rootpk,const mpz_t x1,long p1,const mpz_t n,int k){//root r1,root r2=p^k-r1,calcola radici quadrate di n modulo p alla k
    //parametri:p1 è il primo,k è la potenza a cui elevare p, x1 è la radice quadrata di n mod p
    // calcola rootpk^2===n mod p^k,sapendo x1: x1^2==n mod p
    if(mpz_sgn(x1)<=0 || mpz_sgn(n)<=0){
        handle_error_with_exit("root zero or negative\n");
    }
    if(rootpk==NULL || n==NULL || k<=0 || p1<=1){
        handle_error_with_exit("error in parameter\n");
    }
    if(p1==2){
        handle_error_with_exit("invalid prime\n");
    }
    mpz_t q,r,r2,p,v,w,exp,x;

    mpz_init(q);
    mpz_init(r);
    mpz_init(r2);
    mpz_init(p);
    mpz_init(v);
    mpz_init(w);
    mpz_init(exp);
    mpz_init(x);

    mpz_set_si(p,p1);//p=p
    mpz_set(x,x1);//x=x
    mpz_pow_ui(q,p,k);//q=p^k
    if(mpz_divisible_p(q,p)==0){
        handle_error_with_exit("error in mpz_divisible square_root_mod_p_to_k\n");
    }
    mpz_divexact(r,q,p);//r=q/p

    mpz_mul_ui(r2,r,2);//r2=2r
    mpz_neg(r2,r2);//r2=-2r
    mpz_set(exp,q);//exp=q=p^k
    mpz_add(exp,exp,r2);//exp=q-2r
    mpz_add_ui(exp,exp,1);//exp=q-2r+1
    if(mpz_divisible_ui_p(exp,2)==0){
        handle_error_with_exit("error in mpz_divisible square_root_mod_p_to_k\n");
    }
    mpz_divexact_ui(exp,exp,2);//exp=(q-2r+1)/2
    mpz_powm(v,x,r,q);//v=x^r mod p^k=x^p^(k-1) mod p^k
    mpz_powm(w,n,exp,q);//w=n^exp mod p^k
    mpz_mul(rootpk,v,w);//rootpk=w*v
    mpz_mod(rootpk,rootpk,q);//rootpk=r1 mod p^k

    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(r2);
    mpz_clear(p);
    mpz_clear(v);
    mpz_clear(w);
    mpz_clear(exp);
    mpz_clear(x);

    return;
}

void calculate_root_poly_second_degree_mod_p_to_k(mpz_t j1t,const mpz_t p_to_k,long p,long k,const mpz_t a,const mpz_t b,const mpz_t n){
//calcola la radice del polinomio ax^2+2bx+c mod p^k con gcd(a,p^k)=p,gcd(a/p,p^k)=1
    //x=(-ko)*(a/p)^-1 mod p^k ko=(b-r)/p r^2==n mod p^(k+1),x1^2==n mod p
    if(j1t==NULL || p_to_k==NULL || p<=0 || k<=0 || a==NULL || b==NULL || n==NULL || p%2==0 || (mpz_divisible_ui_p(a,2)!=0) ){
        handle_error_with_exit("calculate_root_poly_second_degree_mod_p_to_k\n");
    }
    int t;
    mpz_t r,x1,p_temp,k0,ptok1,a_temp,inv;
    mpz_init(r);
    mpz_init(x1);
    mpz_init(p_temp);
    mpz_init(k0);
    mpz_init(ptok1);
    mpz_init(a_temp);
    mpz_init(inv);

    mpz_set_si(p_temp,p);//p_temp==p
    t=quadratic_residue(x1,n,p_temp);//x1^2==n mod p
    if(t==-1 || t==0){
        handle_error_with_exit("error in calculate_root_poly_second_degrre_mod_p_to_k\n");
    }
    square_root_mod_p_to_k(r,x1,p,n,k+1);//r^2==n mod p^k+1
    mpz_sub(k0,r,b);//-ko=r-b;
    if(mpz_divisible_ui_p(k0,p)==0){//r-b non è divisibile per p
        mpz_mul_ui(ptok1,p_to_k,p);//ptok=p^(k+1)
        mpz_sub(r,ptok1,r);//r=p^(k+1)-r
        mpz_sub(k0,r,b);//-ko=r-b;
        mpz_divexact_ui(k0,k0,p);//-k0=r-b/p
    }
    else{//r-b divisibile per p
        mpz_divexact_ui(k0,k0,p);//-k0=r-b/p
    }
    mpz_set(a_temp,a);//a_temp=a
    if(mpz_divisible_ui_p(a_temp,p)==0){
        handle_error_with_exit("error in mpz_divisible calculate_root_poly_second_degree_mod_p_to_k\n");
    }
    mpz_divexact_ui(a_temp,a_temp,p);//a_temp=a/p
    mpz_invert(inv,a_temp,p_to_k);//inv=(a/p)^-1 mod p^k
    mpz_mul(j1t,inv,k0);//j1t=inv*ko=(a/p)^-1 mod p^k *(r-b/p)

    mpz_clear(p_temp);
    mpz_clear(r);
    mpz_clear(x1);
    mpz_clear(k0);
    mpz_clear(ptok1);
    mpz_clear(a_temp);
    mpz_clear(inv);
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

void print_array_a_struct(struct a_struct*array_a_struct,int length){
    if(array_a_struct==NULL || length<=0){
        handle_error_with_exit("error in print_aray_a_struct\n");
    }
    for(int i=0;i<length;i++){
        printf("number=%d,index=%d\n",array_a_struct[i].number_prime_a,array_a_struct[i].index_prime_a);
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
void calculate_square(mpz_t square,const mpz_t a,int index,const mpz_t b,const mpz_t n){
    // -M<index<M
    mpz_mul_si(square,a,index);//a*j
    mpz_add(square,square,b);//a*j+b
    mpz_mod(square,square,n);// a*j+b mod n
    return;
}
void adjust_array_bi(mpz_t *array_bi,int s,const mpz_t a){//if 2*b>a ->b=a-b
    if(array_bi==NULL || s<0 || a==NULL){
        handle_error_with_exit("error in adjust array_bi\n");
    }
    mpz_t temp;
    mpz_init(temp);
    long length=pow(2,s-1);
    for(long i=0;i<length;i++){
        mpz_mul_ui(temp,array_bi[i],2);//temp=2*b[i]
        if(mpz_cmp(temp,a)>0){//if 2*b[1]>a
            mpz_sub(array_bi[i],a,array_bi[i]);//b[1]=a-b[1]
        }
    }
    mpz_clear(temp);
    return;
}

void calculate_thresold_a(mpfr_t thresold_a,const mpz_t n,long M){
    if(M<=0 || thresold_a==NULL || n==NULL){
        handle_error_with_exit("error in calculate thresold_a\n");
    }
    mpfr_t rad2,radnn,m;
    mpfr_t t;

    mpfr_init(t);
    mpfr_init(radnn);
    mpfr_init(rad2);
    mpfr_init(m);

    mpfr_set_si(m,M,MPFR_RNDN);//m=M
    mpfr_set_d(rad2,RAD2,MPFR_RNDN);//rad2=1.4142
    mpfr_set_z(radnn,n,MPFR_RNDN);//radn=n
    mpfr_sqrt(radnn,radnn,MPFR_RNDN);//radn=rad(n)
    mpfr_mul(t,radnn,rad2,MPFR_RNDN);//t=rad(n)*rad2
    mpfr_div(t,t,m,MPFR_RNDN);//t2=rad(n)*rad2/M

    mpfr_set(thresold_a,t,MPFR_RNDN);
    mpfr_clear(t);
    mpfr_clear(rad2);
    mpfr_clear(radnn);
    mpfr_clear(m);
    return;
}

void calculate_x0(mpz_t x0,const mpz_t n,int k,char *factorized){
    if(x0 ==NULL || n==NULL || (k!=1 && k!=3 && k!=5 && k!=7) || factorized==NULL){
        handle_error_with_exit("error in calculate x0\n");
    }
    mpz_t s,n_temp;
    mpz_init(s);
    mpz_init(n_temp);

    if(mpz_divisible_ui_p(n,k)==0){
        handle_error_with_exit("error in mpz_divisible calculate x0,n not divisible by k\n");
    }
    mpz_divexact_ui(n_temp,n,k);//n_temp=n/k
    mpz_sqrt(x0,n_temp);//x0=rad(n)
    mpz_mul(s,x0,x0);//s=x0^2
    if(mpz_cmp(s,n_temp)==0){//se s=x0^2=n_temp
        gmp_printf("factorization of n=%Zd*%Zd\n",x0,x0);
        /*fprintf(file_log,"p=");
        mpz_out_str(file_log,10,x0);
        fprintf(file_log," ");
        fprintf(file_log,"q=");
        mpz_out_str(file_log,10,x0);
        fprintf(file_log," ");*/
        *factorized=1;
        mpz_clear(s);
        mpz_clear(n_temp);
        return;
    }

    mpz_sqrt(x0,n);//x0=rad(n)
    mpz_mul(s,x0,x0);//s=x0^2
    if(mpz_cmp(s,n)==0){//se s=x0^2=n
        if(mpz_divisible_ui_p(n,k*k)==0){
            handle_error_with_exit("error in calculate x0,n not divisible by k^2\n");
        }
        mpz_divexact_ui(s,n,k*k);//s=/k
        mpz_sqrt(s,s);
        gmp_printf("factorization of n=%Zd*%Zd*%d\n",s,s,k);
        /*fprintf(file_log,"p=");
        mpz_out_str(file_log,10,s);
        fprintf(file_log," ");
        fprintf(file_log,"q=");
        mpz_out_str(file_log,10,s);
        fprintf(file_log," ");
        fprintf(file_log,"z=%d ",k);*/
        *factorized=1;
        mpz_clear(s);
        mpz_clear(n_temp);
        return;
    }

    mpz_clear(s);
    mpz_clear(n_temp);
    return;
}
void calculate_p_min_p_max_i_f(long*p_min_i,long*p_max_i,struct node_factor_base*head_f_base_f,long cardinality_factor_base){
    if(p_min_i==NULL || p_max_i==NULL || head_f_base_f==NULL || cardinality_factor_base<=0){
        handle_error_with_exit("error in calculate _p_min_p_max_i\n");
    }
    struct node_factor_base*list=head_f_base_f;
    mpz_t temp;
    mpz_init(temp);
    long i=0;
    while(list!=NULL){
        mpz_set_si(temp,list->prime);//valore del primo della factor base
        if(*p_min_i==0 && mpz_cmp_si(temp,SIQS_MIN_PRIME_POLYNOMIAL)>=0){
            *p_min_i=i;
        }
        if(*p_max_i==0 && mpz_cmp_si(temp,SIQS_MAX_PRIME_POLYNOMIAL)>0){
            *p_max_i=i-1;
            break;
        }
        list=list->next;//vai al prossimo elemento
        i++;//incrementa indice
    }
    if(*p_max_i==0){//primi della factor base<SIQS_MAX_PRIME_POLYNOMIAL
        *p_max_i=cardinality_factor_base-1;
    }
    if(*p_min_i<=1){//primi della factor base<SIQS_MIN_PRIME_POLYNOMIAL
        *p_min_i=2;//salta -1 e 2
    }
    if((*p_max_i-*p_min_i)<20){//se la distanza tra p_max e p_min è minore di 20(s_max=20) allora imposta p_min al minimo tra 2 e 5
        *p_min_i=min(*p_min_i,5);
    }
    mpz_clear(temp);
    return;
}

void calculate_target_a1_f(mpfr_t target_a1,const mpfr_t target_a,struct node_factor_base*head_f_base_f,
                           long p_min_i,long p_max_i,int 	cardinality_factor_base){
    if(target_a==NULL || head_f_base_f==NULL || p_min_i<0 || p_max_i<=0 || p_min_i>p_max_i || target_a1==NULL ||
       cardinality_factor_base<=p_max_i){
        handle_error_with_exit("error in calculate_target_a1\n");
    }
    mpz_t prime_min;
    mpz_t prime_max;
    mpz_t prime_sum;
    mpfr_t prime_avg;
    mpfr_t temp;
    mpfr_t sqrt_prime_avg;

    mpz_init(prime_sum);
    mpz_init(prime_min);
    mpz_init(prime_max);
    mpfr_init(prime_avg);
    mpfr_init(sqrt_prime_avg);
    mpfr_init(temp);
    int count=0;
    struct node_factor_base*p=head_f_base_f;
    while(count<p_min_i){
        p=p->next;
        count++;
    }
    if(count!=p_min_i){
        handle_error_with_exit("error in calculate target_a1\n");
    }

    mpz_set_si(prime_min,p->prime);//ottieni prime_min
    while(count<p_max_i){
        p=p->next;
        count++;
    }
    if(count!=p_max_i){
        handle_error_with_exit("error in calculate target_a1\n");
    }
    mpz_set_si(prime_max,p->prime);//ottieni prime_max

    mpz_add(prime_sum,prime_max,prime_min);//prime_sum=prime_min+prime_max
    mpfr_set_z(prime_avg,prime_sum,MPFR_RNDN);//prime_avg=prime_sum
    mpfr_div_2ui(prime_avg,prime_avg,1,MPFR_RNDN);//prime_avg=(pmax+pmin)/2
    mpfr_sqrt(sqrt_prime_avg,prime_avg,MPFR_RNDN);//sqrt_prime=rad((pmax+pmin)/2)
    mpfr_set(temp,target_a,MPFR_RNDN);//temp=target=rad(2*n)/M
    mpfr_div(temp,temp,sqrt_prime_avg,MPFR_RNDN);//temp=target/rad((pmax+pmin)/2)
    mpfr_set(target_a1,temp,MPFR_RNDN);//target_a1=temp=target/rad((pmax+pmin)/2)
    mpz_clear(prime_min);
    mpz_clear(prime_max);
    mpz_clear(prime_sum);
    mpfr_clear(prime_avg);
    mpfr_clear(sqrt_prime_avg);
    mpfr_clear(temp);
    return;
}

void calculate_a_f2(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,int**best_q,int**best_q_number){
    if(s==NULL || target_a==NULL || mpfr_sgn(target_a)<0 || head_f_base_f==NULL || cardinality_factor_base<=0 || a==NULL || best_q==NULL || best_q_number==NULL){
        handle_error_with_exit("error in calculate_a\n");
    }
    long p_min_i=0;//indice del minimo primo da scegliere rispetto alla factor base
    long p_max_i=0;//indice del massimo primo da scegliere rispetto alla factor base
    //p_max_i>p_min_i
    char initialized=0;//se è zero la prima volta che si trova un a,anche se il rapporto a/target è pessimo si
    //prende quell'a
    int value;
    int p_temp2;
    long p_i;//indici dei primi scelti nella factor base per rappresentare a
    long s_max;//valore massimo di s=pmax-pmin+1,numero massimo di primi da scegliere
    mpz_t v;
    int iter=0,iter2=0;
    mpz_init(v);
    if(cardinality_factor_base<3){
        mpz_set_si(a,0);//poni a=0
        *s=0;
        mpz_clear(v);
        *best_q=NULL;
        *best_q_number=NULL;
        return;//ritorna array_of_prime_chosen_for_a==NULL e a=0
    }
    if(cardinality_factor_base==3){
        get_element_linked_list_f(&value,head_f_base_f,2);//a=valore del 3 primo della factor base
        mpz_set_si(a,value);
        *s=1;
        //fprintf(file_log,"p_min=p_max=2 ");
        *best_q_number=alloc_array_int(1);//array che contiene i dei primi scelti
        *best_q=alloc_array_int(1);
        (*best_q_number)[0]=value;
        (*best_q)[0]=2;
        mpz_clear(v);
        return;
    }
    mpz_t p_temp;
    mpfr_t target_a1;
    mpz_t best_a;
    mpfr_t best_ratio;
    double ratio_double=0,best_ratio_double=0;
    mpfr_t ratio;
    mpz_t a2_int;
    mpfr_t a2_double;
    mpfr_t p_rational;
    double dist_best,dist_ratio;
    mpfr_init(p_rational);
    mpz_init(a2_int);
    mpfr_init(a2_double);
    mpz_init(p_temp);
    mpfr_init(target_a1);

    mpz_init(best_a);
    mpfr_init(best_ratio);//double
    mpfr_init(ratio);//double


    mpz_set_si(best_a,0);//best_a=0
    mpfr_set_si(best_ratio,0,MPFR_RNDN);//best_target=0

    calculate_p_min_p_max_i_f(&p_min_i,&p_max_i,head_f_base_f,cardinality_factor_base);
    s_max=p_max_i-p_min_i+1;//nel caso peggiore s=p_max-p_min+1
    if(s_max>S_MAX){
        s_max=S_MAX;
    }
    calculate_target_a1_f(target_a1,target_a,head_f_base_f,p_min_i,p_max_i,cardinality_factor_base);
    if(mpfr_cmp_si(target_a1,1)<=0){//se target_a1<0 poni s=0 e ritorna
        mpz_set_si(a,0);//poni a=0
        *s=0;
        mpz_clear(v);
        mpfr_clear(p_rational);
        mpz_clear(a2_int);
        mpfr_clear(a2_double);
        mpz_clear(p_temp);
        mpfr_clear(target_a1);
        mpz_clear(best_a);
        mpfr_clear(best_ratio);//double
        mpfr_clear(ratio);//double
        return;
    }
    //printf("target_a1=");
    //mpfr_out_str(stdout,10,0,target_a1,MPFR_RNDN);
    //printf("\n");
    //print_time_elapsed("time to calculate target_a1");

    int count=0;//conta quante volte while(mpz_cmp(a,target_a1)<0) è verificata
    int*q=alloc_array_int(s_max);//array che contiene gli indici dei primi scelti
    int*q_number=alloc_array_int(s_max);//array che contiene i dei primi scelti
    *best_q_number=alloc_array_int(s_max);//array che contiene i numeri dei primi scelti
    *best_q=alloc_array_int(s_max);//array che contiene gli indici dei primi scelti
    int length_best_q=0;//inizialmente la lunghezza della lista migliore è zero
    for(int i=0;i<NUM_ITER_FOR_CALCULATE_A;i++){//iterazioni per cercare di migliorare a
        mpz_set_si(a2_int,1);//a2=1 azzera a2
        memset(q,0,sizeof(int)*s_max);//azzera q
        memset(q_number,0,sizeof(int)*s_max);//azzera q
        count=0;//azzera le append a q
        iter=0;
        mpfr_set_z(a2_double,a2_int,MPFR_RNDN);
        while(mpfr_cmp(a2_double,target_a1)<0 && iter<=MAX_ITER){//finquando a2<target_a1 oppure sono state raggiunte tot iterazioni
            p_i=2;
            iter2=0;
            while((p_i==2 || is_in_array_int(q,s_max,p_i)) && iter2<=MAX_ITER2){//scegli un p_i che non è stato ancora scelto
                p_i=rand_long(p_min_i,p_max_i);
                if(p_i>p_max_i || p_i<p_min_i){
                    handle_error_with_exit("error random long generator\n");
                }
                iter2++;
            }
            //numero diverso da 2 e che non è nell'array di q
            if(iter2>=MAX_ITER2){
                break;//ricomincia da capo e setta a2=1
            }
            get_element_linked_list_f(&p_temp2,head_f_base_f,p_i);//prendi l'elemento p_i-esimo dalla factor base
            mpz_set_si(p_temp,p_temp2);
            mpfr_set_z(p_rational,p_temp,MPFR_RNDN);//p_rational=p_temp
            mpfr_mul(a2_double,a2_double,p_rational,MPFR_RNDN);//a2=a2*p
            mpz_mul(a2_int,a2_int,p_temp);//a2=a2*p

            if(count>=s_max){
                handle_error_with_exit("error in calculate a,index out of bounds\n");
            }
            q_number[count]=mpz_get_si(p_temp);
            q[count]=p_i;
            count++;//aumenta il conto delle append fatte
            iter++;
            if(count>=s_max){
                break;
            }
        }
        //una sequenza di numeri è stata scelta vediamo quanto è buona(vedendo i rapporti della migliore attualmente),se è buona la 			aggiorniamo
        mpfr_div(ratio,a2_double,target_a,MPFR_RNDN);//ratio=a2/target_a

        ratio_double=mpfr_get_d(ratio,MPFR_RNDN);
        if(initialized==0){
            initialized=1;
            mpz_set(best_a,a2_int);//best_a=a2
            mpfr_set(best_ratio,ratio,MPFR_RNDN);//best_ratio=ratio
            best_ratio_double=ratio_double;
            memcpy(*best_q,q,sizeof(int)*s_max);//best_q=q
            memcpy(*best_q_number,q_number,sizeof(int)*s_max);//best_q=q
            length_best_q=count;//la lunghezza di best_q coincide con il numero di append fatte
        }
        dist_best=fabs(RATIO_A-best_ratio_double);//distanza di best_ratio da 1
        dist_ratio=fabs(RATIO_A-ratio_double);//distanza di ratio_da 1
        if(dist_ratio<dist_best){
            mpz_set(best_a,a2_int);//best_a=a2
            mpfr_set(best_ratio,ratio,MPFR_RNDN);//best_ratio=ratio
            best_ratio_double=ratio_double;
            memcpy(*best_q,q,sizeof(int)*s_max);//best_q=q
            memcpy(*best_q_number,q_number,sizeof(int)*s_max);//best_q=q
            length_best_q=count;//la lunghezza di best_q coincide con il numero di append fatte
        }
    }

    mpz_set_si(a,1);//a=1
    for(int i=0;i<length_best_q;i++){//moltiplica tutti i fattori di a
        mpz_mul_si(a,a,(*best_q_number)[i]);
    }
    if(mpz_cmp(best_a,a)!=0){
        gmp_printf("a=%Zd\n",a);
        gmp_printf("best_a=%Zd\n",best_a);
        handle_error_with_exit("error in calculate a\n");
    }
    *s=length_best_q;//imposta il valore di s
    if(*s>s_max){
        handle_error_with_exit("error in calculate a,s must be minor or equal to s_max\n");
    }
    if(*s==0 && mpz_cmp(a,0)!=0){
        handle_error_with_exit("error in function calculate_a s or a\n");
    }
    free(q);
    free(q_number);
    mpfr_clear(p_rational);
    mpz_clear(p_temp);
    mpfr_clear(target_a1);
    mpfr_clear(a2_double);
    mpz_clear(a2_int);
    mpz_clear(best_a);
    mpfr_clear(best_ratio);
    mpfr_clear(ratio);
    mpz_clear(v);
    return;
}

mpz_t*calculate_array_Bk_f(int*number_prime_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1){
    mpz_t*array_Bk=NULL;
    long vpk;
    char t;

    mpz_t ak,root,ak_inverse;
    mpz_t pk,n_temp;

    mpz_init(ak);
    mpz_init(root);
    mpz_init(ak_inverse);
    mpz_init(n_temp);
    mpz_init(pk);
    if(n==NULL || a==NULL || b1==NULL || card_factor_base<=0){
        handle_error_with_exit("error in calculate_array_Bk\n");
    }
    if(mpz_get_si(a)==0){//se a=0
        mpz_clear(ak);
        mpz_clear(root);
        mpz_clear(ak_inverse);
        mpz_clear(n_temp);
        mpz_clear(pk);
        return NULL;//ritorna NULL
    }
    if(mpz_sgn(n)<=0 || s<0 || mpz_sgn(a)<=0){
        handle_error_with_exit("error in calculate_array_Bk\n");
    }
    long index=0;
    mpz_set_si(b1,0);//b1==0
    array_Bk=alloc_array_mpz(s);
    for(int k=0;k<s;k++){
        vpk=number_prime_a[index];
        index++;
        mpz_set_si(pk,vpk);//pk=vpk,pk=primo iesimo scelto per a
        if(mpz_divisible_p(a,pk)==0){
            handle_error_with_exit("error in mpz_divisible_calculate array_BK\n");
        }
        mpz_divexact(ak,a,pk);//ak=a/pk
        if(mpz_invert(ak_inverse,ak,pk)==0){
            handle_error_with_exit("error in mpz_invert calculate array_Bk\n");
        }//ak_inverse=(ak)^-1 mod pk
        mpz_set(n_temp,n);//n_temp=n
        mpz_mod(n_temp,n,pk);//n_temp=n mod pk
        if(vpk==2){//la radice di n mod 2 è la riduzione di n mod 2,in realtà i numeri primi devono essere dispari
            mpz_set(root,n_temp);//root =n_temp
        }
        else{
            t=quadratic_residue(root,n_temp,pk);//root=radice quadrata di n modulo pk
            if(t==-1 || t==0){
                handle_error_with_exit("error in calculate quadratic_residue\n");
            }
        }
        mpz_mul(root,root,ak);//root=root*ak
        mpz_mul(root,root,ak_inverse);//tk*ak*(ak)^-1=root*ak*(ak)^-1
        mpz_mod(root,root,a);//Bk può essere ridotto modulo a
        mpz_set(array_Bk[k],root);//array_Bk[k]=root*ak*ak_inverse;
        mpz_add(b1,b1,array_Bk[k]);//b1=(b1)+array_Bk[k],alla fine del for b1=somma di tutti i Bk
    }

    mpz_clear(pk);
    mpz_clear(n_temp);
    mpz_clear(ak);
    mpz_clear(root);
    mpz_clear(ak_inverse);
    return array_Bk;
}
void multiply_n_for_k(mpz_t n,int *k,char*factorized){//n=n*k in modo tale che n è un quadrato modulo 8,cioè n==1 mod 8
    if(n==NULL || mpz_sgn(n)<=0 || n==NULL || k==NULL || factorized==NULL){
        handle_error_with_exit("error in multiply_n_for_k\n");
    }
    if(mpz_divisible_2exp_p(n,1)){
        mpz_divexact_ui(n,n,2);
        gmp_printf("factorization of n=2*%Zd\n",n);
        /*fprintf(file_log,"p=2 ");
        fprintf(file_log,"q=");
        mpz_out_str(file_log,10,n);
        fprintf(file_log," ");*/
        *factorized=1;
        return;
    }
    return;
    /*mpz_t v;//valore temporaneo
    mpz_init_set_si(v,8);
    mpz_mod(v,n,v);
    *k=mpz_get_si(v);//k congruo a n mod 8,se n congruo a 1 mod 8 allora k=1
    mpz_mul(n,n,v);//n=n*k
    mpz_clear(v);
    return;*/
}

int calculate_v(int i){
    int v=1;
    while(i%2==0){//finquando è divisibile per 2
        i=i/2;//dividilo per 2
        v++;
    }
    return v;
}

mpz_t *calculate_bi(mpz_t *array_Bk,const mpz_t b1,int s){
    if(array_Bk==NULL){
        return NULL;
    }
    if(mpz_sgn(b1)<=0 || s<=0){//se array_Bk diverso da NULL s deve essere almeno 1
        handle_error_with_exit("error in calculate_bi\n");
    }
    long length=pow(2,s-1);//length=2^(s-1)
    mpz_t*array_bi=alloc_array_mpz(length);
    int v;
    int exp;
    mpz_t value;
    mpz_init(value);
    mpz_set(array_bi[0],b1);//prima posizione dell'array è b1
    for(int i=1;i<length;i++){
        v=calculate_v(i);
        exp=ceil((double)i/pow(2,v));// parte intera superiore i/2^v
        mpz_set_si(value,2);//value=2
        if((exp & 1)==0){//exp è pari
            //2*(-1)^[i/2^v]==2
        }
        else{
            mpz_set_si(value,-2);//2*(-1)^[i/2^v]==-2
        }
        mpz_mul(value,value,array_Bk[v-1]);//2*(-1)^[i/2^v]*Bv
        //value=2*(pow(-1,exp))*array_Bk[v-1];//2*(-1)^[i/2^v]*Bv
        mpz_add(array_bi[i],array_bi[i-1],value);//bi=bi-1+2*(-1)^[i/2^v]*Bv
    }
    mpz_clear(value);
    return array_bi;
}

void calculate_a_and_b(int*solution,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base,mpz_t a,mpz_t b,const mpz_t n){
    if(solution==NULL || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0 || card_f_base<=0 || num_B_smooth<card_f_base || a==NULL || b==NULL || n==NULL){
        handle_error_with_exit("error in parameter calculate a and b\n");
    }
    //il vettore solution ci dice quali relazioni vanno moltiplicate,ogni elemento di solution che è 1 è una relazione da moltiplicare con le altre
    mpz_t temp,v_temp;
    mpz_t num_square;

    mpz_init(num_square);
    mpz_init(temp);
    mpz_init(v_temp);

    mpz_set_si(a,1);//a=1 temporaneo
    mpz_set_si(b,1);//b=1 temporaneo
    mpz_t*v=create_array_temp_factorization(card_f_base,matrix_B_smooth[0]);//-1 0 2 0 3 0 ecc primi della factor base con esponente zero
    for(int i=0;i<num_B_smooth;i++){
        if(solution[i]==1){//se solution==1 allora bisogna moltiplicare le relazioni
            mpz_set(num_square,matrix_B_smooth[i][card_f_base*2]);//moltiplica le radici quadrate dei numeri B-smooth mod n
            mpz_mul(a,a,num_square);//a=a*square
            mpz_mod(a,a,n);//a=a*square mod n
            sum_elem_multiple_of_2_mpz(v,matrix_B_smooth[i],card_f_base*2,card_f_base*2);//moltiplica le relazioni sommando opportunamente gli esponenti tra i 2 numeri
        }
    }
    divide_elem_multiple_of_2_by_x(v,card_f_base*2,2);//divide gli esponenti per 2 per calcolare successivamente b
    for(int j=2;j<card_f_base*2;j=j+2){//per ogni primo della factor base...
        mpz_set(v_temp,v[j]);//v_temp=v[j]
        mpz_powm(v_temp,v_temp,v[j+1],n);//v_temp=pow(v[j],v[j+1]) mod n
        mpz_mul(b,b,v_temp);//b=primo factor base*esponente del primo della factor base
        mpz_mod(b,b,n);//b mod n
    }

    mpz_clear(temp);
    mpz_clear(num_square);
    mpz_clear(v_temp);
    free_memory_array_mpz(v,card_f_base*2);
    return;
}
struct a_struct*create_array_a_struct(int*number_prime_a,int*index_number_a,int length){
    if(number_prime_a==NULL || index_number_a==NULL || length<=0){
        handle_error_with_exit("error in create_array_a_struct\n");
    }
    struct a_struct*array_a_struct=malloc(sizeof(struct a_struct)*length);
    if(array_a_struct==NULL){
        handle_error_with_exit("error in malloc create array a struct\n");
    }
    for(int i=0;i<length;i++){
        array_a_struct[i].number_prime_a=number_prime_a[i];
        array_a_struct[i].index_prime_a=index_number_a[i];
    }
    return array_a_struct;
}

void calculate_index_min_max_a(int*number_prime_a,int*index_prime_a,int length,int*index_min_a,int*index_max_a){
    if(number_prime_a==NULL || index_max_a==NULL || length<=0 || index_min_a==NULL || index_max_a==NULL){
        handle_error_with_exit("error in calculate min_max_a\n");
    }
    int min_a=number_prime_a[0];
    *index_min_a=index_prime_a[0];
    int max_a=number_prime_a[0];
    *index_max_a=index_prime_a[0];
    for(int i=1;i<length;i++){
        if(min_a>number_prime_a[i]){
            min_a=number_prime_a[i];
            *index_min_a=index_prime_a[i];
        }
        if(max_a<number_prime_a[i]) {
            max_a=number_prime_a[i];
            *index_max_a=index_prime_a[i];
        }
    }
    return;
}
void calculate_thresold_large_prime(mpz_t thresold_large_prime,int max_prime){
    if(max_prime<=0){
        handle_error_with_exit("error in calculate_thresold_prime\n");
    }
    mpz_t temp;
    mpz_init(temp);
    mpz_set_si(temp,max_prime);//temp=max_prime
    mpz_mul(temp,temp,temp);//temp?max_prime^2
    mpz_set(thresold_large_prime,temp);
    mpz_clear(temp);
    return;
}
void add_exponent_of_a(long**matrix_B_smooth,int num_B_smooth,int s,const mpz_t a,long* array_of_prime_chosen_for_a){
    long index;
    if(a==NULL || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0 || s<0 || mpz_sgn(a)<=0){
        handle_error_with_exit("error in add exponent_of_a\n");
    }
    if(mpz_cmp_si(a,1)==0){//se a=1
        if(array_of_prime_chosen_for_a==NULL){//e se array==NULL
            return;
        }
        else{
            handle_error_with_exit("error in parameter add exponent_of_a\n");
        }
    }
    for(int i=0;i<s;i++){//per tutta la lunghezza dell'array dei primi scelti
        index=array_of_prime_chosen_for_a[2*i];//ritorna l'indice del numero primo nella factor base
        for(int j=0;j<num_B_smooth;j++){//per tutte le righe della matrice dei B_smooth
            matrix_B_smooth[j][2*index+1]+=1;//aumenta l'esponente di 1 sui primi scelti
        }
    }
    return;
}

void dg(mpz_t dig,const mpz_t y){//dg==g' derivata di g calcolata nel punto y
    if(dig==NULL || y==NULL){
        handle_error_with_exit("invalid parameter dg\n");
    }
    mpz_mul_2exp(dig,y,1);//dig=2*y
    mpz_add_ui(dig,dig,1);//dig=2y+1
    return ;
}

void g(mpz_t gy,const mpz_t y1,const mpz_t a){//funzione per calcolare le radici quadrate modulo 2^k
    mpz_t a_temp,y_temp;
    if(gy==NULL || y1==NULL || a==NULL){
        handle_error_with_exit("invalid parameter g\n");
    }
    mpz_init(y_temp);
    mpz_init_set(a_temp,a);//a_temp=a

    mpz_mul_2exp(a_temp,a_temp,1);//a_temp=a*2
    mpz_set(y_temp,y1);//y_temp=y1
    mpz_mul(y_temp,y_temp,y_temp);//y_temp=y1^2
    mpz_add(y_temp,y_temp,y1);//y_temp=y^2+y
    mpz_sub(y_temp,y_temp,a_temp);//y=y*y+y-2*a;
    mpz_set(gy,y_temp);

    mpz_clear(y_temp);
    mpz_clear(a_temp);
    return;
}
void calculate_new_y(mpz_t y1,const mpz_t n,long k){//calcola yt dove t=k-2 sapendo yt si può ricavare yk e calcolare la radice quadrata modulo 2^k
    //y y corrisponde a yl dove l=k-3 da yl si puo risalire a yt e quindi a yk
    if(y1==NULL || n==NULL || k<=0 || y1==NULL || mpz_sgn(n)<=0){
        handle_error_with_exit("error in calculate new y\n");
    }
    long t=k-2;//t=k-2

    mpz_t coff_g,two_to_t,gy1,dig;

    mpz_init(coff_g);//a
    mpz_init(two_to_t);
    mpz_init(gy1);//g(y1)
    mpz_init(dig);//g'(y1)

    mpz_set(coff_g,n);//coff_g=n
    mpz_sub_ui(coff_g,coff_g,1);//coff_g=n-1
    if(mpz_divisible_2exp_p(coff_g,3)==0){//verifica divisibiità per 8
        handle_error_with_exit("error in mpz_divisible calculate new y\n");
    }
    mpz_divexact_ui(coff_g,coff_g,8);//coff_g=(n-1)/8=a

    mpz_ui_pow_ui(two_to_t,2,t);//two_to_t=2^t
    dg(dig,y1);//dig=g'(y1)=2y+1
    mpz_invert(dig,dig,two_to_t);//dig=(g'(y1))^-1 mod 2^t
    g(gy1,y1,coff_g);//y=g(y1);
    mpz_mul(dig,gy1,dig);//dig=g(y1)*(g'(y1))^-1
    mpz_sub(y1,y1,dig);//y1-g(y1)*(g'(y1))^-1
    mpz_mod(y1,y1,two_to_t);//y1==y1-g(y1)*(g'(y1))^-1 mod 2^t

    mpz_clear(coff_g);
    mpz_clear(two_to_t);
    mpz_clear(gy1);
    mpz_clear(dig);
    return;
}