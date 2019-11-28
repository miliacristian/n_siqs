
#include "n.h"

void adjust_n(mpz_t n,int *k){
    //divide n per k
    if(n==NULL || k==NULL || *k<=0 || (*k!=1 && *k!=3 && *k!=5 && *k!=7)){
        handle_error_with_exit("error in adjust_n\n");
    }
    if(mpz_divisible_ui_p(n,*k)==0){
        handle_error_with_exit("error in mpz_divisible adjust_n\n");
    }
    mpz_divexact_ui(n,n,*k);//n=n/k
    *k=1;
    return;
}