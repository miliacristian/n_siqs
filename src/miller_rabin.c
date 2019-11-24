#include "miller_rabin.h"
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

char get_and_check_n(mpz_t n,FILE*file_number){
	if(file_number==NULL || n==NULL){
		handle_error_with_exit("error in get and check n,invalid filename\n");
	}
	int digit=mpz_inp_str(n,file_number,10);//leggi dal file n in base 10
	if(mpz_sgn(n)<0){
		mpz_neg(n,n);//n=-n;//inverte il segno di n
	}
	if(mpz_cmp_si(n,0)==0 || mpz_cmp_si(n,1)==0){//se n==0 o n==1 errore
		handle_error_with_exit("n must be different than 0,1 and -1\n");
	}
	test_n(n,NUM_TEST_MILLER_RABIN);//verifica che n non è un numero primo
	return digit;
}
