#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "basic.h"
#include "math_function.h"
#include "dynamic_list.h"
#include "dyn_list_base.h"
#include "matrix_function.h"
#include "print.h"
#include "miller_rabin.h"
#include "criv_quad.h"
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <mpfr.h>
#include "main.h"
#include <pthread.h>

	//valori globali(presi da altri file)
	extern FILE*file_log;//file di log
	extern struct timespec timer;//istante di tempo 
	extern struct timespec time_start;//istante di tempo iniziale

	//valori globali del file main
	struct row_factorization r;
	int k=-1;//moltiplicatore di n,se n dispari k=1,3,5 o 7,se n pari so fattorizzarlo
	//argv[1]=path del file da leggere per ottenere n da fattorizzare
	struct node_f*head_f_base_f=NULL;//testa della lista dinamica factor base
	struct node_f*tail_f_base_f=NULL;//coda della lista dinamica factor base
	mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
	mpz_t a,thresold_q,q;//valore del coefficiente a del polinomio,soglia q e q
	mpfr_t thresold_a;//soglia per il calcolo di a
	mpz_t temp;
	int s=-1;//numero di primi della factor base distinti che compongono a
	long B=-1;//smoothness_bound
	long M=-1;//metà dimensione array
	mpz_t b1;//valore del primo b=somma di tutti i Bk presi dall'array Bk
	mpz_t *array_bi=NULL;
	mpz_t *array_Bk=NULL;
	int dim_sol=-1;//numero di vettori linearmente indpendenti ottenuti dalla risoluzione del sistema lineare,combinandoli opportunamente 		calcolano tutte le possibili combinazioni/soluzioni del sistema
	int cardinality_factor_base=-1;//cardinalità factor base
	struct matrix_factorization**array_matrix_B_smooth=NULL;
	int length_array_matrix=-1;
	int num_increment_M_and_B=0;
	int*index_prime_a=NULL;//indice dei primi usati per ottenere a,rispetto alla factor base
	int*number_prime_a=NULL;//numeri primi usati per ottenere a

int main(int argc,char*argv[]){
	srand((unsigned int)time(NULL));//imposta seme casuale
	if(argc!=2){
		handle_error_with_exit("usage<path>\n");
	}
	FILE*file_number=open_file(argv[1]);//apri file in cui risiede il numero n da fattorizzare
	double mean_increment_M_and_B=0;
	for(int i=0;i<MAX_NUM_FOR_DIGIT;i++){
		//apertura file e dichiarazione variabili
		file_log=open_file_log();//apri file di log
		num_increment_M_and_B=0;
		struct matrix_factorization*matrix_B_smooth=NULL;//matrice che mememorizza le fattorizzazioni dei soli numeri b_smooth 			presenti 			nell'array
		struct matrix_factorization*mat=NULL;
		mpz_t*array_number=NULL;
		mpz_t a_default,b_default;//a,b sono i coefficienti del polinomio aj^2+2bj+c,thresold a serve per calcolare il valore di a
		int*array_id=NULL;
		int num_B_smooth=-1;//numero di numeri b-smooth trovati nell'array
		char*linear_system=NULL;//sistema lineare da risolvere per trovare a e b
		int**base_matrix=NULL;//matrice che riporta per colonna i vettori che formano una base del sistema lineare
		pthread_t *array_tid=NULL;
		int row_result=-1;//numero di righe risultanti dalla concatenazione di tutte le matrici
		int factorizations_founded=-1;//numero di fattorizazzioni trovate
		char factorized=0;
		char digit=-1;//numero cifre di n
		
		k=1;
		//mpz_init
		mpz_init(n);
		mpz_init(a);
		mpz_init(temp);
		mpz_init(x0);
		mpfr_init(thresold_a);
		mpz_init(thresold_q);
		mpz_init(b_default);
		mpz_init(b1);
		mpz_init(a_default);
		mpz_init(q);
		
		mpz_set_si(b1,-1);//b1=-1
	
		//gettime
		gettime(&timer);
		gettime(&time_start);

		//n
		digit=get_and_check_n(n,file_number);
		gmp_printf("n=%Zd,digit=%d\n",n,digit);//stampa n e numero cifre
		fprintf(file_log,"n=");
		mpz_out_str(file_log,10,n);
		fprintf(file_log," ");
		print_time_elapsed("time to calculate n");
	
		//M e B
		M=-1;
		B=-1;
		calculate_best_M_and_B(n,digit,&M,&B);
		printf("M=%ld\n",M);
		fprintf(file_log,"M_start=%ld ",M);
		print_time_elapsed("time to calculate M");
		fprintf(file_log,"B_start=%ld ",B);
		if(mpz_cmp_si(n,B)<=0){// se n è minore o uguale a B imposta B=n-3
			B=mpz_get_si(n)-2;//b non deve essere maggiore di n-2
		}
		printf("B=%ld\n",B);
		print_time_elapsed("time to calculate B");

		//k,moltiplicatore di n per renderlo un quadrato modulo 8
		multiply_n_for_k(n,&k,&factorized);
		if(factorized==1){
			goto clean_memory;
		}
		printf("k=%d\n",k);
		fprintf(file_log,"k=%d ",k);
		gmp_printf("n*k=%Zd\n",n);
		print_time_elapsed("time to calculate k");
	
		if(B<=0 || M<=0){
			handle_error_with_exit("error in B or M\n");
		}
	
		//x0
		calculate_x0(x0,n,k,&factorized);//x0=n^(1/2),xo radice quadrata di n
		if(factorized==1){
			goto clean_memory;
		}
		gmp_printf("x0=%Zd\n",x0);//stampa x0
		print_time_elapsed("time to calculate x0");
		
		//a_default
		mpz_set_si(a_default,1);
		gmp_printf("a_default=%Zd\n",a_default);
	
		//b_default
		mpz_set(b_default,x0);
		gmp_printf("b_default=%Zd\n",b_default);

		while(factorizations_founded<=0){//finquando non sono stati trovati fattori
			adjust_n(n,&k);//aggiusta n per calcolare i suoi fattori(divide n per k)
			multiply_n_for_k(n,&k,&factorized);//k,moltiplicatore di n per renderlo un quadrato modulo 8
			print_time_elapsed("time to adjust & multiply n for k");

			//thresold_a per applicare siqs
			calculate_thresold_a(thresold_a,x0,M);//a è circa rad(2*n)/M
			printf("thresold_a=");
			mpfr_out_str(stdout,10,0,thresold_a,MPFR_RNDN);
			printf("\n");
			print_time_elapsed("time to calculate thresold_a");
	
			//factor base
			head_f_base_f=create_factor_base_f(&cardinality_factor_base,B,&tail_f_base_f,n,q,thresold_q);//crea 				lista dinamica con tutti i primi
			printf("factor base=");
			print_list_f(head_f_base_f);
			printf("cardinality factor_base %d\n",cardinality_factor_base);
			fprintf(file_log,"card_f_base=%d ",cardinality_factor_base);
			print_time_elapsed("time to calculate factor base");

			//a,per siqs e generare tutti gli altri b,prodotto di primi dispari distinti
			calculate_a_f2(a,thresold_a,&s,head_f_base_f,cardinality_factor_base,&index_prime_a,&number_prime_a);
			gmp_printf("a=%Zd\n",a);
			fprintf(file_log,"a=");
			mpz_out_str(file_log,10,a);
			fprintf(file_log," ");
			printf("s=%d\n",s);

			//number_prime_a e index_prime_a
			if(s>0){
				printf("number_prime_a=");
				print_array_int(number_prime_a,s);
				printf("index_prime_a=");
				print_array_int(index_prime_a,s);
			}
			print_time_elapsed("time to calculate a");

			//array Bk
			array_Bk=calculate_array_Bk_f(index_prime_a,number_prime_a,cardinality_factor_base,n,s,a,b1);
			if(array_Bk!=NULL){
				gmp_printf("b1=%Zd\n",b1);
				//print_array_Bk(array_Bk,s);
			}
			print_time_elapsed("time to calculate array Bk");

			//array bi
			array_bi=calculate_bi(array_Bk,b1,s);
			if(array_Bk!=NULL){
				free_memory_array_mpz(array_Bk,s);
				array_Bk=NULL;
			}
			if(array_bi!=NULL){
				reduce_array_mpz_mod_n(array_bi,(int)pow(2,s-1),a);
				adjust_array_bi(array_bi,s,a);
				//print_array_bi(array_bi,s);
			}
			print_time_elapsed("time to calculate array bi");

			//creazione array_matrix_factorization
			if(array_bi==NULL){
				length_array_matrix=1;
			}
			else{
				length_array_matrix=(int)pow(2,s-1)+1;
			}
			if(NUM_THREAD==0){
				length_array_matrix=1;
			}
			array_matrix_B_smooth=alloc_array_matrix_factorization(length_array_matrix);
			printf("length_array_matrix_factorization=%d\n",length_array_matrix);

			//fattorizza numeri nell'array
			r.prime=alloc_array_int(cardinality_factor_base);
			r.log_prime=alloc_array_float(cardinality_factor_base);
			create_row_factorization(head_f_base_f,cardinality_factor_base);
			//print_array_float(r.log_prime,cardinality_factor_base);
			//print_array_int(r.prime,cardinality_factor_base);
			free_memory_list_f(head_f_base_f);
			head_f_base_f=NULL;
			tail_f_base_f=NULL;
			print_time_elapsed("time to create row factorization");
			//creazione e avvio thread
			if(length_array_matrix!=1){
				array_tid=alloc_array_tid(NUM_THREAD);//alloca memoria per contenere tutti i tid
				array_id=create_threads(array_tid,NUM_THREAD);//crea tutti i thread
			}
			print_time_elapsed("time to create thread");
	
			//creazione matrice fattorizazione main thread
			mat=create_matrix_factorization_f(M,cardinality_factor_base,a_default,b_default,n);

			print_time_elapsed("time_to_create matrix_factorization main thread");
			factor_matrix_f(n,M,mat,cardinality_factor_base,a_default,b_default);
			print_time_elapsed("time_to_factor matrix_factorization main thread");
			
			num_B_smooth=count_number_B_smooth_matrix_unsorted_f(mat,2*M+1);
			printf("num_B_smooth=%d\n",num_B_smooth);
			if(num_B_smooth>0){
					matrix_B_smooth=create_matrix_B_smooth_f(num_B_smooth,
					*mat,2*M+1,cardinality_factor_base*2);
			}
			else{
				num_B_smooth=0;
				matrix_B_smooth=alloc_matrix_factorization(0);
			}
			free_memory_matrix_factorization(mat);
			array_matrix_B_smooth[length_array_matrix-1]=matrix_B_smooth;//aggiungi la matrice fattorizzazioni di default 				come ultima matrice	
			//aspetta tutti i thread e libera memoria
			if(length_array_matrix!=1 && NUM_THREAD>0){
				join_all_threads(array_tid,NUM_THREAD);//aspetta tutti i thread
				
				if(array_tid!=NULL){//libera memoria allocata
					free(array_tid);
					array_tid=NULL;
				}
				if(array_id!=NULL){
					free(array_id);
					array_id=NULL;
				}
			}
			if(index_prime_a!=NULL){
				free(index_prime_a);
				index_prime_a=NULL;
			}
			if(number_prime_a!=NULL){
				free(number_prime_a);
				number_prime_a=NULL;
			}
			if(array_bi!=NULL){
				free_memory_array_mpz(array_bi,(int)pow(2,s-1));
				array_bi=NULL;
			}
			free(r.log_prime);
			r.log_prime=NULL;
			printf("threads ended the job\n");
			print_time_elapsed("time to wait all threads");
			concatenate_all_matrix_B_smooth(array_matrix_B_smooth,length_array_matrix,&row_result);
			printf("row_result=%d\n",row_result);
			fprintf(file_log,"num potential B_smooth=%d ",row_result);
			print_time_elapsed("time to create single matrix factorization");
			print_matrix_factorization_f(*(array_matrix_B_smooth[0]));
			//crea matrice_fattorizzazioni_B_smooth
			print_time_elapsed("time to create matrix_B_smooth");
			linear_system=create_linear_system_f(array_matrix_B_smooth[0],cardinality_factor_base);
			free(linear_system);
			free(r.prime);
			r.prime=NULL;
			print_time_elapsed("time_to_create_linear_system");
			free_memory_matrix_factorization(array_matrix_B_smooth[0]);//libera la memoria 
			//della matrice concatenazione di matrici
			array_matrix_B_smooth[0]=NULL;
			free(array_matrix_B_smooth);//libera la memoria dell'array_matrix_mpz
			array_matrix_B_smooth=NULL;
			break;
			/*//crea sistema lineare dalla matrix_B_smooth
			linear_system=create_linear_system(cardinality_factor_base,matrix_B_smooth,num_B_smooth);
			//print_linear_system(linear_system,cardinality_factor_base,num_B_smooth);
			if(check_if_matrix_is_reduce_mod_n(linear_system,cardinality_factor_base,num_B_smooth,2)==0){
				handle_error_with_exit("error in main,create linear system\n");
			}
			print_time_elapsed("time to calculate linear system");
			//calcola una base della matrice
			base_matrix=calculate_base_linear_system(linear_system,cardinality_factor_base,num_B_smooth,&dim_sol);
			if(base_matrix==NULL){//non ci sono abbastanza soluzioni,ricomincia il crivello quadratico
				free_memory_matrix_mpz(matrix_B_smooth,num_B_smooth,cardinality_factor_base*2+1);
				matrix_B_smooth=NULL;
				calculate_news_M_and_B(&M,&B);
				continue;
			}
			if(check_if_matrix_is_reduce_mod_n(base_matrix,num_B_smooth,dim_sol,2)==0){
				handle_error_with_exit("error in main,calculate base_linear_system\n");
			}
			//print_base_linear_system(base_matrix,num_B_smooth,dim_sol);
			printf("dim_sol=%d\n",dim_sol);
			print_time_elapsed("time to calculate base linear system");
			if(check_solution_base_matrix(linear_system,cardinality_factor_base,num_B_smooth,base_matrix,num_B_smooth,dim_sol)==0){
				handle_error_with_exit("error in main,invalid solution\n");
			}
			//prova a fattorizzare n
			adjust_n(n,&k);//aggiusta n per calcolare i suoi fattori(divide n per k)
			factorizations_founded=find_factor_of_n_from_base_matrix(base_matrix,num_B_smooth,&dim_sol,linear_system,
			cardinality_factor_base,num_B_smooth,n,matrix_B_smooth,num_B_smooth,cardinality_factor_base);
			free_memory_matrix_int(linear_system,cardinality_factor_base,num_B_smooth);
			linear_system=NULL;
			free_memory_matrix_int(base_matrix,num_B_smooth,dim_sol);
			base_matrix=NULL;
			print_time_elapsed("time to calculate all solution linear system");
			free_memory_matrix_mpz(matrix_B_smooth,num_B_smooth,cardinality_factor_base*2+1);
			matrix_B_smooth=NULL;
			printf("factorizations founded=%d\n",factorizations_founded);
			print_time_elapsed("time to calculate factor of n from solutions");
			if(factorizations_founded==0){//nessuna fattorizzazione trovata,ricomincia il crivello quadratico
				calculate_news_M_and_B(&M,&B);
				continue;
			}
			mean_increment_M_and_B+=num_increment_M_and_B;*/
		}
		clean_memory:
		mpz_clear(n);
		mpz_clear(a);
		mpz_clear(temp);
		mpz_clear(x0);
		mpfr_clear(thresold_a);
		mpz_clear(thresold_q);
		mpz_clear(b_default);
		mpz_clear(b1);
		mpz_clear(a_default);
		mpz_clear(q);
		fprintf(file_log,"M_end=%ld ",M);
		fprintf(file_log,"B_end=%ld ",B);
		fprintf(file_log,"num_increment=%d ",num_increment_M_and_B);
		//tempo totale,imposta il tempo iniziale alla struct,tempo totale=get_time-tempo iniziale
		timer.tv_nsec=time_start.tv_nsec;//timer=time_start
		timer.tv_sec=time_start.tv_sec;//timer=time_start
		print_time_elapsed_on_file_log("time_total");
		print_time_elapsed("time_total");
		fprintf(file_log,"\n");
		if(fclose(file_log)!=0){
			handle_error_with_exit("error in close file_log\n");
		}
	}
	file_log=open_file_log();//apri file di log
	fprintf(file_log,"mean_increment=%lf\n",mean_increment_M_and_B/MAX_NUM_FOR_DIGIT);
	if(fclose(file_log)!=0){
			handle_error_with_exit("error in close file_log\n");
		}
	if(fclose(file_number)!=0){
		handle_error_with_exit("error in close file_number\n");
	}
	return 0;
}

int thread_job_criv_quad(int id_thread){//id inizia da 0,il lavoro di un thread rimane uguale anche se non si riesce a fattorizzare n
	if(id_thread+1>length_array_matrix-1){//l'indice del thread eccede il numero di job da fare
		return 0;
	}
	int count=id_thread;//indica quale polinomio deve usare per fare il crivello quadratico
	int num_B_smooth=-1;
	struct matrix_factorization*matrix_B_smooth=NULL;
	struct matrix_factorization *matrix=NULL;
	while(count<=length_array_matrix-2){//ogni thread prende un sottoinsieme di compiti,il thread con id 0 farà i compiti 				0,NUM_THREAD,2*NUM_THREAD,il thread 1 farà 1,NUM_THREAD+1,2*NUM_THREAD+1 ecc

		//creazione e fattorizzazione matrice dei thread secondari
		printf("thread=%d\n",count);
		matrix=create_matrix_factorization_f(M,cardinality_factor_base,a,array_bi[count],n);
		//print_matrix_factorization_f(*matrix);
		print_time_elapsed("time_to_create matrix_factorization");
		factor_matrix_f(n,M,matrix,cardinality_factor_base,a,array_bi[count]);
		//print_matrix_factorization_f(*matrix);
		print_time_elapsed("time to factor matrix_factorization ");
		num_B_smooth=count_number_B_smooth_matrix_unsorted_f(matrix,2*M+1);
		printf("num_B_smooth thread=%d\n",num_B_smooth);
		if(num_B_smooth>0){
			matrix_B_smooth=create_matrix_B_smooth_f(num_B_smooth,
			*matrix,2*M+1,cardinality_factor_base*2);
		}
		else{
			num_B_smooth=0;
			matrix_B_smooth=alloc_matrix_factorization(0);
		}
		array_matrix_B_smooth[count]=matrix_B_smooth;
		free_memory_matrix_factorization(matrix);
		count+=NUM_THREAD;//modulo numero dei thread
	}
	return 0;
}
	
