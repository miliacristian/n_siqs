#include "properties.h"
long k=-1;
double thresold_relation;
int num_increment_M_and_B;
char combined;
mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
struct factor_base_data*thread_factor_base_data=NULL;
mpz_t a_old,a_new;//valore del coefficiente a del polinomio
struct a_struct*array_a_struct=NULL;
int s=-1;//numero di primi della factor base distinti che compongono a

void handle_error_with_exit(char*error_string){//uccide il processo dopo essersi accorto di un errore
    if(error_string==NULL){
	    printf("error string is NULL\n");
	    perror("");
        exit(EXIT_FAILURE);
    }
    printf("%s",error_string);
    exit(EXIT_FAILURE);
}
void check_variable_in_defines(){
    if(NUM_TEST_MILLER_RABIN<=0){
		handle_error_with_exit("error in value num_test_miller_rabin\n");
    }
    if(RAD2<=0){
		handle_error_with_exit("error in value rad2\n");
    }
    if(THRESOLD_PRINT_ARRAY<=0){
		handle_error_with_exit("error in value thresold_print_array");
    }
	if(THRESOLD_PRINT_MATRIX<=0){
		handle_error_with_exit("error in value thresold_print_matrix\n");
	}
	if(THRESOLD_PRINT_LIST<=0){
		handle_error_with_exit("error in value thresold print_list\n");
	}
	if(SIQS_MIN_PRIME_POLYNOMIAL<=0){
		handle_error_with_exit("error in value siqs_min_prime_polynomial\n");
	}
	if(SIQS_MAX_PRIME_POLYNOMIAL<=0){
		handle_error_with_exit("error in value siqs_max_prime_polynomial\n");
	}
	if(NUM_ITER_FOR_CALCULATE_A<=0){
		handle_error_with_exit("error in value num_iter_for_calculate_a\n");
    }
    if(MAX_ITER<=0){
		handle_error_with_exit("error in value max_iter\n");
    }
	if(MAX_ITER2<=0){
		handle_error_with_exit("error in value max_iter2\n");
	}
	if(RATIO_A<=0){
		handle_error_with_exit("error in value ratio_a\n");
	}
	if(NUM_THREAD_FACTOR_BASE<0){
		handle_error_with_exit("error in value num_thread_factor_base\n");
	}
	if(NUM_THREAD_POLYNOMIAL<0){
		handle_error_with_exit("error in value num_thread_polynomial\n");
	}
	if(S_MAX<=0){
		handle_error_with_exit("error in value s_max\n");
	}
	if(MAX_DIM_SOL<=0){
		handle_error_with_exit("error in value max_dim_sol\n");
	}
	if(PERC_INCREMENT_M<=0){
		handle_error_with_exit("error in value perc_increment_m\n");
	}
	if(PERC_INCREMENT_B<=0){
		handle_error_with_exit("error in value perc_increment b\n");
	}
	if(NUM_OF_N_TO_FACTORIZE<1){
		handle_error_with_exit("error in value num_of_n_to_factorize\n");
	}
	if(ENOUGH_RELATION<0){
		handle_error_with_exit("error in value enough relation\n");
	}
	if(THRESOLD_B<=0){
		handle_error_with_exit("error in value thresold_b\n");
	}
	if(BIT_OF_UNSIGNED_LONG<8 || BIT_OF_UNSIGNED_LONG%8!=0){
		handle_error_with_exit("error in value bit_of_unsigned_long\n");
	}
	if(TEST!=0 && TEST!=1){
        handle_error_with_exit("error in value TEST\n");
    }
	return;
}
