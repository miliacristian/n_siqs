//
// Created by cristian on 25/11/19.
//

#include "a_b_c_BK_functions.h"
#include "properties.h"
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