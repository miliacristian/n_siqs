#ifndef PROP_H
#define PROP_H


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
#define PERC_INCREMENT_M 50 //quanto aumenta in percentuale M se non si riesce a fattorizzare n
#define PERC_INCREMENT_B 20 //quanto aumenta in percentuale B se non si riesce a fattorizzare n
#define TEST 0 //test=0 nessuna verifica test=1 molte verifiche

#define NUM_OF_N_TO_FACTORIZE 1
#define ENOUGH_RELATION 1.0 //numero maggiore o uguale a 1 indica quante relazioni
//vanno trovate in più rispetto alla cardinalità della factor base num_b_smooth>cardinality*enough_relation
#define ERROR_LOG 25//errore del logaritmo,aumentare per trovare più numeri B_smooth potenziali ma maggior computazione,
//diminuire per trovare meno numeri B_smooth potenziali ma minor computazione,valore default=25
#define THRESOLD_B 20000 //se B è minore di thresold b non dividere il processo di creazione factor base in più thread
#define BIT_OF_UNSIGNED_LONG (8*sizeof(unsigned long))//numero di bit di una variabile unsigned long

#include <stdio.h>
#include <stdlib.h>
void handle_error_with_exit(char*error_string);

#endif
