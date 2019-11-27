#include <stdio.h>
#include <stdlib.h>
#include "factorization_functionsv2.h"
#include <math.h>
#include <gmp.h>
#include "a_b_c_BK_functions.h"
#include "factor_base_functions.h"
#include "mpz_functions.h"

//scriverla in modo parallelizzato
struct row_factorization r;

mpz_t*create_array_temp_factorization(int card_f_base,mpz_t*row_matrix_b_smooth){// crea array -1 0 2 0 30...primo factor base 0...di 2*cardinalità 		della factor base elementi
    if(card_f_base<=0 || row_matrix_b_smooth==NULL){
        handle_error_with_exit("error in parameter create_array_temp_factorization\n");
    }
    mpz_t*v;
    v=alloc_array_mpz(card_f_base*2);//alloca memoria
    for(int i=0;i<card_f_base*2;i=i+2){//ciclo a salti di 2
        mpz_set(v[i],row_matrix_b_smooth[i]);//sui posti dispari (indice pari) viene messo il primo della factor base
    }
    return v;
}
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s){
    if(head_f_base_f==NULL || card_f_base<=0 || (array_a_struct==NULL && s>0) || s<0){
        handle_error_with_exit("error in create row factorization\n");
    }
    mpz_t temp;
    int index=0;
    mpz_init(temp);
    struct node_factor_base*p=head_f_base_f;
    for(int i=0;i<card_f_base;i++){
        r.prime[i]=p->prime;//metti il primo della factor base in posizione pari
        r.root_n_mod_p[i]=p->root_n_mod_prime;
        r.root2_n_mod_p[i]= r.prime[i]-r.root_n_mod_p[i];//rad2=p-rad1 mod p
        if((i==0 && s>0) || (i==1 && s>0) ||
           (index<s && i==array_a_struct[index].index_prime_a) ) {//p divide a non esiste inverso modulo p
            r.inverse_a_mod_p[i]=-1;
            if(i!=0 && i!=1) {
                index++;
            }
        }
        else if(s==0){//a=1
            r.inverse_a_mod_p[i]=1;
        }
        else{//p non divide a ->esiste l'inverso
            mpz_set_si(temp,r.prime[i]);//temp=p
            if(mpz_invert(temp,a,temp)==0){//temp=inverse of a mod p
                if(i==array_a_struct[index].index_prime_a){
                    handle_error_with_exit("error in inverse\n");
                }
                printf("i=%d,s=%d\n",i,s);
                handle_error_with_exit("error in mpz_invert create_row_factorization\n");
            }
            r.inverse_a_mod_p[i]=mpz_get_si(temp);
        }
        if(p->prime!=-1){
            r.log_prime[i]=(int)round(log2f((float)p->prime));
        }
        p=p->next;
    }
    if(index!=s){
        handle_error_with_exit("error in initialize inverse mod p create_row_factorization\n");
    }
    mpz_clear(temp);
    return;
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
void print_factorization(const mpz_t num,struct node_factorization*head_factor){
    if(head_factor==NULL){
        printf("no simple factorization found\n");
        return;
    }
    gmp_printf("%Zd=",num);
    struct node_factorization*p=head_factor;
    while(p!=NULL){
        printf("%d^%d,",p->number,p->exp_of_number);
        p=p->next;
    }
    printf("\n");
    return;
}
void print_matrix_factorization(mpz_t**matrix_factorization,int M,int cardinality_factor_base){
    if(matrix_factorization==NULL || *matrix_factorization==NULL){
        printf("matrix factorization is empty\n");
    }
    if(M <=0 || cardinality_factor_base <=0){
        handle_error_with_exit("error in print_matrix_factorization\n");
    }
    printf("matrix_factorization:\n");
    print_matrix_mpz(matrix_factorization,2*M+1,cardinality_factor_base*2+2);
    return;
}