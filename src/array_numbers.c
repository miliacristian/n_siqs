#include "array_numbers.h"

void print_array_number(mpz_t*array,int M){
    if(M<=0){
        handle_error_with_exit("error in print_array_number\n");
    }
    if(array==NULL){
        printf("array_number is empty\n");
        return;
    }
    if(not_print_array(2*M+1)==1){
        return;
    }
    printf("array number:\n");
    print_array_mpz(array,2*M+1);
    return;
}
char divide_all_by_min1(mpz_t*array_of_number,long M,mpz_t**matrix_factorization){//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione
    if(array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || M<=0){
        handle_error_with_exit("invalid parameter divide all by min1\n");
    }
    for(int i=0;i<2*M+1;i++){
        if(mpz_sgn(array_of_number[i])<0){//finquando gli elementi sono negativi
            mpz_neg(array_of_number[i],array_of_number[i]);//inverti il segno array_of_number[i]=-array_of_number[i];
            mpz_add_ui(matrix_factorization[i][1],matrix_factorization[i][1],1);//imposta la molteplicità di meno uno ad 1
        }
    }
    return 0;
}

