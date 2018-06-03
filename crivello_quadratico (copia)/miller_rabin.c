#include "miller_rabin.h"
#include <stdio.h>
#include "basic.h"
#include <gmp.h>

void test_n(const mpz_t n,int num_test){
	if(n==NULL || num_test<=0){
		handle_error_with_exit("error in test_n\n");
	}
	int test=mpz_probab_prime_p(n,num_test);
	if(test==1){
		handle_error_with_exit("n is probable prime\n");
	}
	else if(test==2){
		handle_error_with_exit("n is prime\n");
	}
	else{
		printf("n is not prime\n");
	}
	return;
}
