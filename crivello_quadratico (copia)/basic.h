#define NUM_TEST_MILLER_RABIN 10
#define RAD2 1.41421356
//#define THRESOLD_PROBABLE_PRIME 1000000 //per i numeri più piccoli della soglia si vede se sono primi dividendoli per ogni numero minore di essi che risiede nella factor base
#define THRESOLD_PROBABLE_PRIME 100
#define THRESOLD_PRINT_ARRAY 1000000000
#define THRESOLD_PRINT_MATRIX 100000000
#define SIQS_MIN_PRIME_POLYNOMIAL 400 //parametri per il calcolo di a
#define SIQS_MAX_PRIME_POLYNOMIAL 8000 //parametri per il calcolo di a
#define NUM_ITER_FOR_CALCULATE_A 100 //parametri per il calcolo di a
#define THRESOLD_RELATION 0 
#define NUM_THREAD 4 //numero di thread
#define S_MAX 10//corrisponde a 2^(S_MAX-1) polinomi diversi
#define MAX_DIM_SOL 15 //dimensione soluzione massima
#define PERC_INCREMENT_M 50 
#define PERC_INCREMENT_B 50
#define MAX_ITER 100 //iterazioni massime per calcolare a
#define MAX_ITER2 100 //iterazioni massime per calcolare a
#define MAX_NUM_FOR_DIGIT 1
#define PERCENT_B_SMOOTH 0.95
#include <gmp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef LINE_H
#define LINE_H
struct timespec;
struct row_factorization{
	int length;//card*2
	int*prime;//contiene primo factor base e valore del suo logaritmo
	float*log_prime;
};
struct row{
	int index_first_prime;
	int index_last_prime;
	mpz_t num;
	mpz_t square;
	float sum_log;//somma dei logaritmi
	float log;//logaritmo di num
};
struct matrix_factorization {//matrice che contiene tutte le fattorizazzioni
	struct row*row;//2*M+1
	int num_row;
};
struct row_relation_B_smooth{//riga che contiene relazioni B_smooth
	int*factorization;
	mpz_t square;
	int length;
};
struct matrix_relation_B_smooth{//matrice che contiene le relazioni B_smooth
	struct row_relation_B_smooth*row;
	int num_row;
};

void print_time_elapsed(char*string);
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
#endif