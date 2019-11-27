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
char verify_factorization(const mpz_t num,mpz_t residuos,struct node_factorization*head_factor,const mpz_t a){
    //num*a deve essere uguale a fattorizzazione*residuo
    mpz_t temp,num_temp;
    long exp_of_factor,factor;
    mpz_t factor_raise_to_exp;
    if(TEST==0){
        return 1;
    }
    mpz_init(num_temp);
    mpz_init(temp);
    mpz_init(factor_raise_to_exp);
    mpz_set(num_temp,num);
    if(mpz_cmp_si(num_temp,0)<0){
        mpz_neg(num_temp,num_temp);
    }
    mpz_set_si(temp,1);//temp=1
    if(head_factor==NULL){
        if(mpz_cmp(num,residuos)!=0){
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            mpz_clear(num_temp);
            printf("residuo !=num\n");
            return 0;
        }
        mpz_clear(temp);
        mpz_clear(factor_raise_to_exp);
        mpz_clear(num_temp);
        return 1;
    }
    struct node_factorization*p=head_factor;
    while(p!=NULL){
        mpz_set_si(factor_raise_to_exp,1);
        factor=p->number;
        //factor
        if(factor<0 && factor!=-1){
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            mpz_clear(num_temp);
            printf("fattore minore di 0\n");
            return 0;
        }
        if(factor==0){
            printf("fattore 0\n");
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            mpz_clear(num_temp);
            return 0;
        }
        if(factor==-1){
            factor=1;
        }
        //exp of factor
        exp_of_factor=p->exp_of_number;
        if(exp_of_factor<=0){
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            mpz_clear(num_temp);
            printf("esponente minore o uguale a 0\n");
            return 0;
        }
        mpz_ui_pow_ui (factor_raise_to_exp,factor,exp_of_factor);
        mpz_mul(temp,temp,factor_raise_to_exp);
        p=p->next;
    }
    if(mpz_cmp_si(residuos,0)<=0){
        printf("residuo minore o uguale a 0\n");
        mpz_clear(temp);
        mpz_clear(factor_raise_to_exp);
        mpz_clear(num_temp);
        return 0;
    }
    mpz_mul(temp,temp,residuos);//moltiplica temp per il residuo alla fine
    mpz_mul(num_temp,num_temp,a);//num*a
    if(mpz_cmp_si(temp,0)<=0){
        handle_error_with_exit("error in function verify factorization\n");
    }
    if(mpz_cmp(temp,num_temp)!=0){//se non sono uguali nega temp
        gmp_printf("temp e num diversi,residuos=%Zd,num=%Zd,temp=%Zd,a=%Zd\n",residuos,num,temp,a);
        mpz_clear(temp);
        mpz_clear(factor_raise_to_exp);
        mpz_clear(num_temp);
        print_factorization(num,head_factor);
        return 0;
    }
    mpz_clear(temp);
    mpz_clear(factor_raise_to_exp);
    mpz_clear(num_temp);
    return 1;
}
char calculate_num_from_factorization(mpz_t num_temp,struct node_factorization*head_factor){
    long exp_of_factor,factor;
    mpz_t factor_raise_to_exp;
    mpz_t temp;
    mpz_init(temp);
    mpz_init(factor_raise_to_exp);

    mpz_set_si(temp,1);//temp=1
    if(head_factor==NULL){
        mpz_clear(temp);
        mpz_clear(factor_raise_to_exp);
        printf("fattorizzazione vuota\n");
        return 0;
    }
    struct node_factorization*p=head_factor;
    while(p!=NULL){
        mpz_set_si(factor_raise_to_exp,1);
        factor=p->number;
        //factor
        if(factor<0 && factor!=-1){
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            printf("fattore minore di 0\n");
            return 0;
        }
        if(factor==0){
            printf("fattore 0\n");
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            return 0;
        }
        //exp of factor
        exp_of_factor=p->exp_of_number;
        if(exp_of_factor<=0){
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
            printf("esponente minore o uguale a 0\n");
            return 0;
        }
        if(factor==-1) {
            if(exp_of_factor!=1) {
                handle_error_with_exit("exponent of -1 greater than 1\n");
            }
            else {
                mpz_neg(temp,temp);
                p=p->next;
                continue;
            }
        }
        mpz_ui_pow_ui (factor_raise_to_exp,factor,exp_of_factor);
        mpz_mul(temp,temp,factor_raise_to_exp);
        p=p->next;
    }
    mpz_set(num_temp,temp);
    mpz_clear(temp);
    mpz_clear(factor_raise_to_exp);
    return 1;
}