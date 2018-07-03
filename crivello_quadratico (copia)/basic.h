#define NUM_TEST_MILLER_RABIN 10 //numero di test di miller_rabin per vedere se il numero è primo
#define RAD2 1.41421356 //radice di 2
#define THRESOLD_PRINT_ARRAY 1000//dopo il valore di soglia l'array non viene stampato
#define THRESOLD_PRINT_MATRIX 1000//dopo il valore di soglia la matrice non viene stampata
#define THRESOLD_PRINT_LIST 1000 //dopo il valore di soglia la lista non viene stampata

#define SIQS_MIN_PRIME_POLYNOMIAL 400 //parametri per il calcolo di a
#define SIQS_MAX_PRIME_POLYNOMIAL 4000 //parametri per il calcolo di a
#define NUM_ITER_FOR_CALCULATE_A 30 //parametri per il calcolo di a,a diversi che vengono generati,si sceglie quello migliore
#define MAX_ITER 3000 //iterazioni massime per calcolare a
#define MAX_ITER2 3000 //iterazioni massime per calcolare a,
#define RATIO_A 1.0 //rapporto tra "a" ideale e "a" trovato,se ratio=1 si vuole a trovato molto molto vicino ad a ideale

#define NUM_THREAD_FACTOR_BASE 8//numero thread supplementari per costruire factor base
#define NUM_THREAD_POLYNOMIAL 10 //numero di thread supplementari per calcolare i polinomi di siqs
#define S_MAX 11//corrisponde a 2^(S_MAX-1) polinomi diversi,limita il numero di polinomi possibili
#define MAX_DIM_SOL 16 //dimensione soluzione massima del sistema lineare,si hanno a disposizione 2^(max_dim_sol)-1 soluzioni diverse
#define PERC_INCREMENT_M 5 //quanto aumenta in percentuale M se non si riesce a fattorizzare n
#define PERC_INCREMENT_B 5 //quanto aumenta in percentuale B se non si riesce a fattorizzare n

#define NUM_OF_N_TO_FACTORIZE 1
#define THRESOLD_RELATION 0.4
#define ENOUGH_RELATION 1.0 //numero maggiore o uguale a 1 indica quante relazioni
//vanno trovate in più rispetto alla cardinalità della factor base num_b_smooth>cardinality*enough_relation
#define ERROR_LOG 25//errore del logaritmo,aumentare per trovare più numeri B_smooth potenziali ma maggior computazione,
//diminuire per trovare meno numeri B_smooth potenziali ma minor computazione,valore default=25
#define THRESOLD_B 20000 //se B è minore di thresold b non dividere il processo di creazione factor base in più thread
#define BIT_OF_UNSIGNED_LONG (8*sizeof(unsigned long))//numero di bit di una variabile unsigned long

#include <gmp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef LINE_H
#define LINE_H

//strutture dati
struct row_factorization{
    int length;
	int*prime;//contiene primo factor base e
	int*log_prime;//valore del suo logaritmo
    int*root_n_mod_p;//prima radice di n mod p
    int*root2_n_mod_p;//seconda radice di n mod p
    int*inverse_a_mod_p;//a^-1 mod p
};

struct thread_data {
	float log_thresold;//valore della soglia
	mpz_t b;//coefficiente b
    int num_potential_B_smooth;
    int num_B_smooth;
    int num_semi_B_smooth;
	struct number*numbers;//array di struct di numberi
	struct node_square_relation*head_square;//ogni thread dopo il suo compito ha una lista di relazioni quadratiche
	struct node_square_relation*tail_square;//ogni thread dopo il suo compito ha una lista di relazioni quadratiche
	struct node_square_relation*head_residuos;
    struct node_square_relation*tail_residuos;
    long*j1_mod_p;//j1=-b+r)*a^-1 mod p
    long*j2_mod_p;//j2=(-b+r2)*a^-1 mod p

	//ogni volta che un thread analizza un dato appende la lista delle relazioni quadratiche a quella
	// precedentemente calcolata
};
struct a_struct {
    int number_prime_a;
    int index_prime_a;
};

struct number {
	int j;//indice j va da -M a M
	int sum_log;//somma del logaritmo
	int first_index_f_base;//primo indice del primo della factor base che divide number,-1 non è considerato
	int last_index_f_base;//ultimo indice del primo della factor base che divide number,-1 non è considerato
	//number=a^2j^2+2*a*b*j+b^2-n,b^2-n=a*c
};
struct square_relation {
    mpz_t square;
    struct node_factorization*head_factorization;
    mpz_t num;
    mpz_t residuos;
};
struct factor_base_data {
	struct node_factor_base*head;
	struct node_factor_base*tail;
	int cardinality_factor_base;
	int last_prime_factor_base;
};
struct factorization_thread_data{
    int id_thread;
    int start;
    int end;
    char is_a_default;
    pthread_mutex_t *mtx;
    void*pointer;
    struct thread_data thread_data;
};

struct timespec print_time_elapsed(char*string);
void handle_error_with_exit(char*error_string);
char array_is_fill_of_value(char*combination,int length,char value);
void free_memory_matrix_int(int **matrix,int num_row,int num_col);
void free_memory_matrix_long(long **matrix,int num_row,int num_col);
void free_memory_array_mpz(mpz_t*array_number,long length);
void gettime(struct timespec*timer);
char get_and_check_n(mpz_t n,FILE*file_number);
void free_memory_matrix_mpz(mpz_t**matrix,int num_row,int num_col);
int* create_threads(pthread_t*array_tid,int num_thread);
pthread_t *alloc_array_tid(int num_thread);
void join_all_threads(pthread_t*array_tid,int length_array);
void test();
FILE*open_file_log();
void print_time_elapsed_on_file_log(char*string);
FILE*open_file(char*path);
struct thread_data*alloc_array_polynomial_thread_data(int length_array_thread_data,long M);
void print_time_elapsed_local(char*string,struct timespec*timer_thread);
void clear_struct_thread_data(struct thread_data t_data,int M);
void free_array_thread_data(struct thread_data*thread_data,int length_array_thread_data);
void free_memory_list_square_relation(struct node_square_relation*head);
struct factor_base_data*alloc_array_factor_base_data(int length);
int* create_factor_base_threads(pthread_t*array_tid,int num_thread);
struct factorization_thread_data* create_factorization_threads(pthread_t*array_tid,struct thread_data thread_data,const mpz_t a,int num_thread);
void free_memory_matrix_char(char **matrix,int num_row,int num_col);
void initialize_mtx(pthread_mutex_t *mtx);
void lock_mtx(pthread_mutex_t *mtx);
void unlock_mtx(pthread_mutex_t *mtx);
void destroy_mtx(pthread_mutex_t *mtx);
void free_memory_matrix_unsigned_long(unsigned long **matrix,int num_row,int num_col);
void free_list_factorization(struct node_factorization*head_factorization);
int calculate_start_factor_base(int id_thread);
void check_variable_in_defines();
int compare_a_struct( const void* a, const void* b);
void print_estimated_time(int cardinality_factor_base,int num_B_smooth);
#endif
