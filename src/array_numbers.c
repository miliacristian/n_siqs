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

