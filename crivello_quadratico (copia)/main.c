#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "basic.h"
#include "math_function.h"
#include "dynamic_list.h"
#include "list_factor_base.h"
#include "matrix_function.h"
#include "print.h"
#include "miller_rabin.h"
#include "criv_quad.h"
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <mpfr.h>
#include "main.h"
#include "list_factorization.h"
#include "list_square_relation.h"
#include <pthread.h>

	//valori globali(presi da altri file)
	//extern FILE*file_log;//file di log
	extern struct timespec timer;//istante di tempo 
	extern struct timespec time_start;//istante di tempo iniziale
	extern struct timespec timer_test;

	//valori globali del file main
	struct row_factorization r;//contiene tutti i primi della factor base e i relativi log
	int k=-1;//moltiplicatore di n,se n dispari k=1,3,5 o 7,se n pari so fattorizzarlo
	//argv[1]=path del file da leggere per ottenere n da fattorizzare
	struct node_factor_base*head_f_base_f=NULL;//testa della lista dinamica factor base
	struct node_factor_base*tail_f_base_f=NULL;//coda della lista dinamica factor base
	mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
	mpz_t a_old,a_new;//valore del coefficiente a del polinomio
	mpfr_t thresold_a;//soglia per il calcolo di a
	mpz_t temp;//mpz temporaneo
	int s=-1;//numero di primi della factor base distinti che compongono a
	long B=-1;//smoothness_bound
	long M=-1;//metà dimensione array
	mpz_t b1;//valore del primo b=somma di tutti i Bk presi dall'array Bk
	mpz_t *array_bi=NULL;//array dei coefficienti bi
	mpz_t *array_Bk=NULL;//array che serve per calcolare array_bi
	int dim_sol=-1;//numero di vettori linearmente indpendenti ottenuti dalla risoluzione del sistema lineare,combinandoli opportunamente 		calcolano tutte le possibili combinazioni/soluzioni del sistema
	int cardinality_factor_base=-1;//cardinalità factor base
	struct thread_data*thread_polynomial_data=NULL;//struttura dati globale per far svolgere la computazione ai thread
	struct factor_base_data*thread_factor_base_data=NULL;
	int num_thread_job=-1;//lunghezza dell'array di matrici,ogni thread riceve una matrice per fare la computazione
	int num_increment_M_and_B=0;
	int*index_prime_a=NULL;//indice dei primi usati per ottenere a,rispetto alla factor base
	int*number_prime_a=NULL;//numeri primi usati per ottenere a
    struct a_struct*array_a_struct=NULL;
    int combined_relations=0;
    mpz_t thresold_large_prime;


int main(int argc,char*argv[]){
	srand((unsigned int)time(NULL));//imposta seme casuale
	if(argc!=2){//se non c'è esattamente un parametro,termina,serve il path in cui c'è scritto il numero all'interno
		handle_error_with_exit("usage<path>\n");
	}
	check_variable_in_defines();
	FILE*file_number=open_file(argv[1]);//apri file in cui risiede il numero n da fattorizzare
	for(int i=0;i<NUM_OF_N_TO_FACTORIZE;i++){

		//dichiarazione variabili
		num_increment_M_and_B=0;
		cardinality_factor_base=0;
		mpz_t a_default,b_default;//a,b sono i coefficienti del polinomio aj^2+2bj+c,thresold a serve per calcolare il valore di a
		int*array_id=NULL;//array che contiene gli id dei nuovi thread creati a partire da 0
		int num_B_smooth=0,num_semi_B_smooth=0,num_potential_B_smooth=0;//numero di numeri b-smooth potenziali e reali trovati nell'array
		char*linear_system=NULL;//sistema lineare da risolvere per trovare a e b
		unsigned long**binary_linear_system=NULL;//sistema lineare binario da ridurre a scala
		int num_col_binary_matrix,num_col_linear_system;//numero colonne sistema binario e sistema lineare
		int**base_matrix=NULL;//matrice che riporta per colonna i vettori che formano una base del sistema lineare
		pthread_t *array_tid=NULL;//array di tid per fare le join dal main thread
		int factorizations_founded=-1;//numero di fattorizazzioni trovate
		char factorized=0;//indica se il numero è stato fattorizzato o no
		char digit=-1;//numero cifre di n
		struct node_square_relation*head_residuos=NULL,*tail_residuos=NULL;
        struct node_square_relation*head_square=NULL,*tail_square=NULL;
		struct node_square_relation *head_sort_residuos=NULL,*tail_sort_residuos=NULL;//contengono tutte le relazioni quadratiche
		struct node_square_relation *head_sort_square=NULL,*tail_sort_square=NULL;
        int last_prime_factor_base=1;//indica l'ultimo primo della factor base,e quindi da quale primo si inizia a creare la factor base se una factor base è già esistente
                //se last_prime==1 allora si aggiungono alla factor base anche -1 e 2
        char factor_base_already_exist=0;//indica se la factor base è già esistente oppure no
        int start=0,number_cycle=0;
		k=1;//moltiplicatore di n

		//mpz_init
		mpz_init(n);
		mpz_init(a_old);
		mpz_init(a_new);
		mpz_init(temp);
		mpz_init(x0);
		mpfr_init(thresold_a);
		mpz_init(b_default);
		mpz_init(b1);
		mpz_init(a_default);
		mpz_init(thresold_large_prime);
		
		mpz_set_si(b1,-1);//b1=-1
        mpz_set_si(a_old,0);
        mpz_set_si(a_new,0);
	
		//gettime
		gettime(&timer);
		gettime(&timer_test);
		gettime(&time_start);

		//n
		digit=get_and_check_n(n,file_number);
		if(i>0){
		    digit=digit-1;
		}
		gmp_printf("n=%Zd,digit=%d\n",n,digit);//stampa n e numero cifre
		print_time_elapsed("time to calculate n");
	
		//M e B
		M=-1;
		B=-1;
		calculate_best_M_and_B(n,digit,&M,&B);
		printf("M=%ld\n",M);
		print_time_elapsed("time to calculate M");
		if(mpz_cmp_si(n,B)<=0){// se n è minore o uguale a B imposta B=n-3
			B=mpz_get_si(n)-2;//b non deve essere maggiore di n-2
		}
		printf("B=%ld\n",B);
		print_time_elapsed("time to calculate B");
        if(B<=0 || M<=0){
            handle_error_with_exit("error in B or M\n");
        }

		//k,moltiplicatore di n per renderlo un quadrato modulo 8
		multiply_n_for_k(n,&k,&factorized);
		if(factorized==1){//se n viene fattorizzato pulire la memoria e terminare il programma
			goto clean_memory;
		}
		printf("k=%d\n",k);
		gmp_printf("n*k=%Zd\n",n);
		print_time_elapsed("time to calculate k");

		//x0=rad(n)
		calculate_x0(x0,n,k,&factorized);//x0=n^(1/2),xo radice quadrata di n
		if(factorized==1){//se n viene fattorizzato pulire la memoria e terminare il programma
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
		factor_base_already_exist=0;

		while(factorizations_founded<=0){//finquando non sono stati trovati fattori:
		    // calcola a=p1*p2*...ps,unisci le relazioni quadratiche vedi se puoi calcolare il sistema lineare,
            // trova souzioni sistema lineare e trova tutti gli a,b del crivello quadratico
			printf("inizio ciclo numero %d\n",number_cycle);
            gmp_printf("n=%Zd,digit=%d\n",n,digit);//stampa n e numero cifre
			adjust_n(n,&k);//aggiusta n per calcolare i suoi fattori(divide n per k)
			multiply_n_for_k(n,&k,&factorized);//k,moltiplicatore di n per renderlo un quadrato modulo 8
			print_time_elapsed("time to adjust & multiply n for k");

			//thresold_a per applicare siqs
			calculate_thresold_a(thresold_a,n,M);//a è circa rad(2*n)/M
			printf("thresold_a=");
			mpfr_out_str(stdout,10,0,thresold_a,MPFR_RNDN);
			printf("\n");
			print_time_elapsed("time to calculate thresold_a");

			//factor base
			if(factor_base_already_exist==0 && B>THRESOLD_B && NUM_THREAD_FACTOR_BASE>0) {//se la factor base non è mai stata creata
				// e se B è maggiore del valore soglia e numero thread per creare factor base>0
				thread_factor_base_data = alloc_array_factor_base_data(NUM_THREAD_FACTOR_BASE + 1);//alloca struttura condivisa a tutti i thread,
				        // ogni thread ha il propri indice
				array_tid = alloc_array_tid(NUM_THREAD_FACTOR_BASE);//alloca memoria per contenere tutti i tid
				array_id = create_factor_base_threads(array_tid, NUM_THREAD_FACTOR_BASE);//crea tutti i thread
				join_all_threads(array_tid, NUM_THREAD_FACTOR_BASE);//aspetta tutti i thread
				if (array_tid != NULL) {//libera memoria allocata
					free(array_tid);
					array_tid = NULL;
				}
				if (array_id != NULL) {
					free(array_id);
					array_id = NULL;
				}
				start=calculate_start_factor_base(NUM_THREAD_FACTOR_BASE);//calcola lo start del main thread,l'end è pari a B
				last_prime_factor_base=start;//inizia a creare la factor base da start
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);//crea la factor base da lastprime fino a B
				printf("time to create factor base main thread");
				for(int i=1;i<NUM_THREAD_FACTOR_BASE;i++){//unisci tutte le liste nella prima lista tranne la lista del main thread
					union_list_factor_base(&(thread_factor_base_data[0].head),&(thread_factor_base_data[0].tail),&(thread_factor_base_data[0].cardinality_factor_base),&(thread_factor_base_data[0].last_prime_factor_base),
										   (thread_factor_base_data[i].head),(thread_factor_base_data[i].tail),thread_factor_base_data[i].cardinality_factor_base,thread_factor_base_data[0].last_prime_factor_base);
				}
				//unisci la prima lista con la lista del main thread
				union_list_factor_base(&(thread_factor_base_data[0].head),&(thread_factor_base_data[0].tail),&(thread_factor_base_data[0].cardinality_factor_base),&(thread_factor_base_data[0].last_prime_factor_base),
                      head_f_base_f,tail_f_base_f,cardinality_factor_base,last_prime_factor_base);
				//reimposta valori globali
				head_f_base_f=thread_factor_base_data[0].head;
				tail_f_base_f=thread_factor_base_data[0].tail;
				cardinality_factor_base=thread_factor_base_data[0].cardinality_factor_base;
				last_prime_factor_base=thread_factor_base_data[0].last_prime_factor_base;
				print_time_elapsed("time to union all lists factor base");

				//libera memoria
				if(thread_factor_base_data!=NULL){
					free(thread_factor_base_data);
					thread_factor_base_data=NULL;
				}
                factor_base_already_exist=1;
			}
			else if(B<=THRESOLD_B || NUM_THREAD_FACTOR_BASE==0){//se B minore della soglia o numero thread factor base=0
				//crea factor base da 0 a B,-1 e 2 sono già presenti,alle prossime iterazioni calcola la lista appendendo la sottolista
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);//crea lista dinamica con tutti i primi
                //imposta factor base come creata
				factor_base_already_exist=1;
			}
			else if(factor_base_already_exist==1){//factor base è già stata creata parti da start=last_prime_factor_base fino a B e appendi la sottolista alla lista precedente
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);//crea lista dinamica con tutti i primi
			}
			else{
				handle_error_with_exit("caso factor base non gestito\n");
			}
			printf("factor base=");
			print_list_factor(head_f_base_f,cardinality_factor_base);
			printf("cardinality factor_base %d\n",cardinality_factor_base);
			print_time_elapsed("time to calculate factor base");

			//verifica che la factor base è corretta
            if(verify_factor_base(head_f_base_f,cardinality_factor_base,last_prime_factor_base)==0){
                handle_error_with_exit("error in main verify factor base\n");
            }
            if(verify_cardinality_list_factor_base(head_f_base_f,cardinality_factor_base)==0){
            	handle_error_with_exit("invalid dimension list factor base\n");
            }
			print_time_elapsed("time to verify factor base");
			//a,per siqs e generare tutti gli altri b,prodotto di primi dispari distinti
			calculate_a_f2(a_new,thresold_a,&s,head_f_base_f,cardinality_factor_base,&index_prime_a,&number_prime_a);
            if(s==0 && mpz_cmp_si(a_new,0)!=0){
                handle_error_with_exit("error in main calculate a 1\n");
            }
            if(s>0 && mpz_cmp_si(a_new,0)==0){
                handle_error_with_exit("error in main calculate a 2\n");
            }
			while(s>0 && (mpz_cmp(a_old,a_new)==0 && mpz_cmp_si(a_new,0)!=0)){//continua fino a quando non trovi un a diverso
				increment_M_and_B(&M,&B);//aumenta M e B
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);
                if(index_prime_a!=NULL){
                    free(index_prime_a);
                    index_prime_a=NULL;
                }
                if(number_prime_a!=NULL){
                    free(number_prime_a);
                    number_prime_a=NULL;
                }
				calculate_thresold_a(thresold_a,n,M);
				printf("thresold_a=");
				mpfr_out_str(stdout,10,0,thresold_a,MPFR_RNDN);
				printf("\n");
				calculate_a_f2(a_new,thresold_a,&s,head_f_base_f,cardinality_factor_base,&index_prime_a,&number_prime_a);
                if(s==0 && mpz_cmp_si(a_new,0)!=0){
                    handle_error_with_exit("error in main calculate a 1\n");
                }
                if(s>0 && mpz_cmp_si(a_new,0)==0){
                    handle_error_with_exit("error in main calculate a 2\n");
                }
            }
            mpz_set(a_old,a_new);//imposta a_old=a_new
			if(s==0 && mpz_cmp_si(a_new,0)==0){
				increment_M_and_B(&M,&B);//aumenta M e B
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);
				if(verify_factor_base(head_f_base_f,cardinality_factor_base,last_prime_factor_base)==0){
					handle_error_with_exit("error in main verify factor base\n");
				}
				if(verify_cardinality_list_factor_base(head_f_base_f,cardinality_factor_base)==0){
					handle_error_with_exit("invalid dimension list factor base\n");
				}
				if(index_prime_a!=NULL){
					free(index_prime_a);
					index_prime_a=NULL;
				}
				if(number_prime_a!=NULL){
					free(number_prime_a);
					number_prime_a=NULL;
				}
			}
			if(s>0) {
				//crea array a struct e ordina i fattori di a
                array_a_struct = create_array_a_struct(number_prime_a, index_prime_a, s);
                qsort(array_a_struct, (size_t)s, sizeof(struct a_struct), compare_a_struct);
                print_array_a_struct(array_a_struct, s);
            }
			gmp_printf("a_old=%Zd,s=%d\n",a_old,s);
			print_time_elapsed("time to calculate a and array_a_struct");
			if(s>0) {
                printf("index_min_a=%d,index_max_a=%d\n", array_a_struct[0].index_prime_a,
                       array_a_struct[s - 1].index_prime_a);
            }
			print_time_elapsed("time to calculate index_min_a index_max_a");

			//number_prime_a e index_prime_a
			if(s>0){
				print_array_a_struct(array_a_struct,s);
			}

			//array Bk,può essere ridotto modulo a
			array_Bk=calculate_array_Bk_f(number_prime_a,cardinality_factor_base,n,s,a_old,b1);
			if(array_Bk!=NULL){
				gmp_printf("b1=%Zd\n",b1);
				print_array_Bk(array_Bk,s);
			}
			print_time_elapsed("time to calculate array Bk");

			//array bi
			array_bi=calculate_bi(array_Bk,b1,s);
			if(array_Bk!=NULL){
				free_memory_array_mpz(array_Bk,s);
				array_Bk=NULL;
			}
			if(array_bi!=NULL){//bi può essere ridotto modulo a
				reduce_array_mpz_mod_n(array_bi,(int)pow(2,s-1),a_old);
				adjust_array_bi(array_bi,s,a_old);
				print_array_bi(array_bi,s);
			}
			print_time_elapsed("time to calculate array bi");

			//creazione struttura dati dei thread
			if(array_bi==NULL){//se array_bi==NULL non ci sono polinomi,numero dei thread è 1
				num_thread_job=1;
			}
			else{
                num_thread_job=(int)pow(2,s-1)+1;//abbiamo 2^(s-1) polinomi diversi con un a fissato
			}
			if(NUM_THREAD_POLYNOMIAL==0){//se non ci sono thread la lunghezza della matrice è 1
				num_thread_job=1;
			}
			printf("num_thread_job=%d\n",num_thread_job);

			//alloca struttura dati dei thread per fare il sieving
            thread_polynomial_data=alloc_array_polynomial_thread_data(NUM_THREAD_POLYNOMIAL+1,M);
			print_time_elapsed("time to create thread data");

			//creazione della struttura row_factorization
			// row_factorization contiene primi factor base,log per ogni primo radice 1 e radice 2 di n mod p e a^-1 mod p)
			r.prime=alloc_array_int(cardinality_factor_base);
			r.log_prime=alloc_array_int(cardinality_factor_base);
			r.root_n_mod_p=alloc_array_int(cardinality_factor_base);
			r.root2_n_mod_p=alloc_array_int(cardinality_factor_base);
			r.inverse_a_mod_p=alloc_array_int(cardinality_factor_base);
			create_row_factorization(head_f_base_f,cardinality_factor_base,a_old,array_a_struct,s);
			printf("logaritmi(arrotondati) factor base:");
			print_array_int(r.log_prime,cardinality_factor_base);
            printf("factor base:");
			print_array_int(r.prime,cardinality_factor_base);
			printf("radice di n mod p:");
			print_array_int(r.root_n_mod_p,cardinality_factor_base);
            printf("inverse a mod p:");
			print_array_int(r.inverse_a_mod_p,cardinality_factor_base);
			print_time_elapsed("time to create row factorization");

			//thresold_large_prime
			calculate_thresold_large_prime(thresold_large_prime,r.prime[cardinality_factor_base-1]);
			gmp_printf("thresold_large_prime=%Zd\n",thresold_large_prime);
			print_time_elapsed("time to calculate thresold_large_prime");
			//creazione e avvio thread per fase sieving
			if(num_thread_job!=1){
				array_tid=alloc_array_tid(NUM_THREAD_POLYNOMIAL);//alloca memoria per contenere tutti i tid
				array_id=create_threads(array_tid,NUM_THREAD_POLYNOMIAL);//crea tutti i thread
			}
			print_time_elapsed("time to create thread");
            mpz_set(thread_polynomial_data[NUM_THREAD_POLYNOMIAL].b,b_default);//imposta b
			//fattorizza numeri nell'array lungo 2m+1
            printf("main thread\n");
			factor_matrix_f(n,M,(thread_polynomial_data[NUM_THREAD_POLYNOMIAL]),cardinality_factor_base,a_default,array_a_struct,s);//fattorizza numeri
			print_time_elapsed("time_to_factor matrix_factorization main thread");
			//print_thread_data(thread_polynomial_data[NUM_THREAD_POLYNOMIAL],M,cardinality_factor_base);

			//log_thresold main thread
			thread_polynomial_data[NUM_THREAD_POLYNOMIAL].log_thresold=calculate_log_thresold(n,M);
			printf("log_thresold main thread=%f\n",thread_polynomial_data[NUM_THREAD_POLYNOMIAL].log_thresold);

			//trova relazioni quadratiche o semi_B_smooth e ordinale per numero
			find_list_square_relation(thread_polynomial_data[NUM_THREAD_POLYNOMIAL],&num_B_smooth,&num_semi_B_smooth,&num_potential_B_smooth,M,&head_square,&tail_square,&head_residuos,&tail_residuos,n,a_default,NULL,0);
			//printf("square\n");
			//print_list_square_relation(head_square,num_B_smooth);
			//printf("residuos\n");
			//print_list_square_relation(head_residuos,num_B_smooth);
			print_time_elapsed("time_to find_list_square_relation main thread");

			//aspetta tutti i thread e libera memoria
			if(num_thread_job!=1 && NUM_THREAD_POLYNOMIAL>0){
				join_all_threads(array_tid,NUM_THREAD_POLYNOMIAL);//aspetta tutti i thread
				
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
			free(array_a_struct);
			array_a_struct=NULL;

			printf("threads ended the job\n");
			timer_test=print_time_elapsed("time to wait all threads");

			printf("num_potential_B_smooth_main_thread=%d,num_B_smooth_main_thread=%d,num_semi_B_smooth main thread=%d\n",
                   num_potential_B_smooth,num_B_smooth,num_semi_B_smooth);

			for(int i=0;i<NUM_THREAD_POLYNOMIAL;i++){//metti tutte le relazioni in head e tail,somma tutti i numeri B_smooth e semi_B_smooth
				//le liste vengono unite in modo tale che quella finale è ordinata per numero
			    union_list_square(&head_square,&tail_square,
								  thread_polynomial_data[i].head_square,thread_polynomial_data[i].tail_square);
                union_list_residuos(&head_residuos,&tail_residuos,
                                  thread_polynomial_data[i].head_residuos,thread_polynomial_data[i].tail_residuos);
				num_B_smooth+=thread_polynomial_data[i].num_B_smooth;
				num_potential_B_smooth+=thread_polynomial_data[i].num_potential_B_smooth;
				num_semi_B_smooth+=thread_polynomial_data[i].num_semi_B_smooth;
				printf("num_potential_B_smooth=%d,num_B_smooth=%d,num_semi_B_smooth=%d\n",num_potential_B_smooth,num_B_smooth,num_semi_B_smooth);
			}
            printf("num_potential_B_smooth=%d,num_B_smooth=%d,num_semi_B_smooth=%d\n",num_potential_B_smooth,num_B_smooth,num_semi_B_smooth);

            //unisici tutte le relazioni quadrate(solamente quelle B_smooth)alla lista delle relazioni quadratiche,
            // la lista finale conterrà relazioni quadratiche ordinate per square
            add_square_relation_to_list_sorted(&head_sort_square,&tail_sort_square,head_square);
            print_time_elapsed("time to add square relation to list sorted");
            head_square=NULL;
            tail_square=NULL;
            if(verify_sorted_square_rel_list(head_sort_square)==0){
                handle_error_with_exit("error in sort relation by square\n");
            }
			if(verify_cardinality_list_square_relation(head_sort_square,num_B_smooth)==0){
            	handle_error_with_exit("error in cardinality square relation\n");
            }
            //unisici tutte le relazioni semi_B_smooth alla lista delle relazioni semi_B_smooth,
            // la lista finale conterrà relazioni semi_B_smooth ordinate per residuo

            /*add_relation_semi_B_smooth_to_list(&head_sort_residuos,&tail_sort_residuos,head_residuos);//commentare questa riga
            print_time_elapsed("time to add relation semi_B_smooth");
            if(verify_sorted_residuos_square_rel_list(head_sort_residuos)==0){
                handle_error_with_exit("error in sort relation by square\n");
            }
            head_residuos=NULL;
            tail_residuos=NULL;

			factorizations_founded = combine_relation_B_smooth_and_semi_B_smooth(&head_sort_square,
																				 &tail_sort_square, head_sort_residuos, n, &num_B_smooth, &num_semi_B_smooth,&combined_relations);
			print_time_elapsed("time_to_combine relation_B_smooth");
			//riassegna la lista delle relazioni quadratiche a head e tail
			head_sort_residuos = NULL;
			tail_sort_residuos = NULL;
			if (factorizations_founded == 1) {
				break;
			}*///commentare questa riga

            //trova nuove relazioni quadratiche con un nuovo square,una nuova fattorizazzione e imposta num=0
            if(num_B_smooth>=cardinality_factor_base*THRESOLD_RELATION) {
				add_relation_semi_B_smooth_to_list(&head_sort_residuos,&tail_sort_residuos,head_residuos);
				print_time_elapsed("time to add relation semi_B_smooth");
				head_residuos=NULL;
				tail_residuos=NULL;
				//factorizations_founded = combine_relation_B_smooth_and_semi_B_smooth(&head_sort_square,
				//		&tail_sort_square, head_sort_residuos, n, &num_B_smooth, &num_semi_B_smooth,&combined_relations);
				factorizations_founded = combine_relation_B_smooth_and_semi_B_smooth_v2(&head_sort_square,
																					 &tail_sort_square, &head_sort_residuos,&tail_sort_residuos, n, &num_B_smooth, &num_semi_B_smooth,&combined_relations);
				print_time_elapsed("time_to_combine relation_B_smooth");
				//riassegna la lista delle relazioni quadratiche a head e tail
				head_sort_residuos = NULL;
				tail_sort_residuos = NULL;
                if (factorizations_founded == 1) {
                    break;
                }
				if(verify_cardinality_list_square_relation(head_sort_square,num_B_smooth)==0){
					handle_error_with_exit("error in cardinality square relation\n");
				}
			}
            //print_list_square_relation(head,num_B_smooth);
            if(verify_sorted_square_rel_list(head_sort_square)==0){
                handle_error_with_exit("error in sorted list by square\n");
            }
            printf("card_f_base=%d\n",cardinality_factor_base);
            if(num_B_smooth>=cardinality_factor_base*ENOUGH_RELATION){
                remove_same_square(&head_sort_square,&tail_sort_square,&num_B_smooth,&num_semi_B_smooth);
                print_time_elapsed("time to remove same square");
                if(verify_sorted_square_rel_list(head_sort_square)==0){
                    handle_error_with_exit("error in sort list by square\n");
                }
				if(verify_cardinality_list_square_relation(head_sort_square,num_B_smooth)==0){
					handle_error_with_exit("error in cardinality square relation\n");
				}
            }
            //verifica che il numero di relazioni trovate è ancora sufficiente dopo aver rimosso gli square uguali
            if(num_B_smooth<cardinality_factor_base*ENOUGH_RELATION){
                calculate_news_M_and_B(&M,&B);
				free_array_thread_data(thread_polynomial_data,NUM_THREAD_POLYNOMIAL+1);
				thread_polynomial_data=NULL;
				if(r.log_prime!=NULL && r.prime!=NULL && r.root_n_mod_p!=NULL && r.inverse_a_mod_p!=NULL &&
						r.root2_n_mod_p!=NULL) {
					free(r.log_prime);
					r.log_prime = NULL;
					free(r.prime);
					r.prime=NULL;
					free(r.root_n_mod_p);
					r.root_n_mod_p=NULL;
					free(r.root2_n_mod_p);
					r.root2_n_mod_p=NULL;
					free(r.inverse_a_mod_p);
					r.inverse_a_mod_p=NULL;
				}
				if(thread_polynomial_data!=NULL) {
					free_array_thread_data(thread_polynomial_data, NUM_THREAD_POLYNOMIAL + 1);
					thread_polynomial_data = NULL;
				}
				number_cycle++;
                continue;
            }

			//algebra step:sistema lineare
            binary_linear_system=create_binary_linear_system(head_sort_square,cardinality_factor_base,num_B_smooth,&num_col_binary_matrix);
            //print_binary_matrix(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            print_time_elapsed("time_to_create_binary linear_system");
            reduce_echelon_form_binary_matrix(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            //print_binary_matrix(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            print_time_elapsed("time_to_reduce echelon form linear_system");
            linear_system=from_matrix_binary_to_matrix_char(binary_linear_system,cardinality_factor_base,num_col_binary_matrix,&num_col_linear_system);
            print_time_elapsed("time_to_copy matrix binary into matrix char");
			//print_linear_system(linear_system,cardinality_factor_base,num_col_linear_system);
            free_memory_matrix_unsigned_long(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            binary_linear_system=NULL;

            //algebra step:base sistema lineare
            base_matrix=calculate_base_linear_system_char(linear_system,cardinality_factor_base,num_col_linear_system,&dim_sol);
            if(base_matrix==NULL){//non ci sono abbastanza soluzioni,ricomincia il crivello quadratico
				calculate_news_M_and_B(&M,&B);
                free(linear_system);
                linear_system=NULL;
				if(r.log_prime!=NULL && r.prime!=NULL && r.root_n_mod_p!=NULL && r.inverse_a_mod_p!=NULL &&
														 r.root2_n_mod_p!=NULL) {
					free(r.log_prime);
					r.log_prime = NULL;
					free(r.prime);
					r.prime=NULL;
					free(r.root_n_mod_p);
					r.root_n_mod_p=NULL;
					free(r.root2_n_mod_p);
					r.root2_n_mod_p=NULL;
					free(r.inverse_a_mod_p);
					r.inverse_a_mod_p=NULL;
				}
				if(thread_polynomial_data!=NULL) {
					free_array_thread_data(thread_polynomial_data, NUM_THREAD_POLYNOMIAL + 1);
					thread_polynomial_data = NULL;
				}
				number_cycle++;
				continue;
			}
			/*if(check_if_matrix_is_reduce_mod_n(base_matrix,num_B_smooth,dim_sol,2)==0){
				handle_error_with_exit("error in main,calculate base_linear_system\n");
			}*/
			printf("dim_sol=%d\n",dim_sol);
			print_time_elapsed("time to calculate base linear system");
			printf("combined_relations=%d\n",combined_relations);
			/*if(check_solution_base_matrix_char(linear_system,cardinality_factor_base,num_col_linear_system,
            base_matrix,num_B_smooth,dim_sol)==0){
				handle_error_with_exit("error in main,invalid solution\n");
			}
			print_time_elapsed("time to check solution base linear system");*/

            //algebra step:calcolo di tutti gli a,b del crivello quadratico
			factorizations_founded=find_factor_of_n_from_base_matrix_char(base_matrix,num_col_linear_system,&dim_sol,
			linear_system,cardinality_factor_base,num_col_linear_system,n,head_sort_square,num_B_smooth,cardinality_factor_base);
			free(linear_system);
            linear_system=NULL;
			free_memory_matrix_int(base_matrix,num_col_linear_system,dim_sol);
			base_matrix=NULL;
			print_time_elapsed("time to calculate solution from base linear system");
			if(factorizations_founded==0){
                calculate_news_M_and_B(&M,&B);
				if(r.log_prime!=NULL && r.prime!=NULL && r.root_n_mod_p!=NULL && r.inverse_a_mod_p!=NULL &&
														 r.root2_n_mod_p!=NULL) {
					free(r.log_prime);
					r.log_prime = NULL;
					free(r.prime);
					r.prime=NULL;
					free(r.root_n_mod_p);
					r.root_n_mod_p=NULL;
					free(r.root2_n_mod_p);
					r.root2_n_mod_p=NULL;
					free(r.inverse_a_mod_p);
					r.inverse_a_mod_p=NULL;
				}
				if(thread_polynomial_data!=NULL) {
					free_array_thread_data(thread_polynomial_data, NUM_THREAD_POLYNOMIAL + 1);
					thread_polynomial_data = NULL;
				}
				number_cycle++;
                continue;
			}
		}
		clean_memory://pulire memoria rimanente
		if(head_f_base_f!=NULL) {
			free_memory_list_f(head_f_base_f);
			head_f_base_f = NULL;
			tail_f_base_f = NULL;
		}
		if(r.log_prime!=NULL && r.prime!=NULL && r.root_n_mod_p!=NULL && r.inverse_a_mod_p!=NULL &&
												 r.root2_n_mod_p!=NULL) {
			free(r.log_prime);
			r.log_prime = NULL;
			free(r.prime);
			r.prime=NULL;
			free(r.root_n_mod_p);
			r.root_n_mod_p=NULL;
			free(r.root2_n_mod_p);
			r.root2_n_mod_p=NULL;
			free(r.inverse_a_mod_p);
			r.inverse_a_mod_p=NULL;
		}
		if(thread_polynomial_data!=NULL) {
			free_array_thread_data(thread_polynomial_data, NUM_THREAD_POLYNOMIAL + 1);
			thread_polynomial_data = NULL;
		}
		if(head_sort_square!=NULL) {
			free_memory_list_square_relation(head_sort_square);
			head_sort_square = NULL;
		}

		//mpz_clear
		mpz_clear(n);
		mpz_clear(a_old);
		mpz_clear(a_new);
		mpz_clear(temp);
		mpz_clear(x0);
		mpfr_clear(thresold_a);
		mpz_clear(b_default);
		mpz_clear(b1);
		mpz_clear(a_default);
		mpz_clear(thresold_large_prime);

		//tempo totale,imposta il tempo iniziale alla struct,tempo totale=get_time-tempo iniziale
		timer.tv_nsec=time_start.tv_nsec;//timer=time_start
		timer.tv_sec=time_start.tv_sec;//timer=time_start
		print_time_elapsed("time_total");
	}
	if(fclose(file_number)!=0){
		handle_error_with_exit("error in close file_number\n");
	}
	printf("combined_relations=%d\n",combined_relations);
	return 0;
}

int thread_job_to_create_factor_base(int id_thread){
    long remainder=reduce_int_mod_n_v2(B,NUM_THREAD_FACTOR_BASE+1);//rem=b mod num_thread
	long length=(B-remainder)/(NUM_THREAD_FACTOR_BASE+1);
	int start=id_thread*length+1;//se id=o start=1
	int end=start+length-1;
	//es remainder=0 thread=5 B=500.000 -> len=100.000 start=0*100000+1,end=1+100000-1=100000,start2=100001,end2=200000
	thread_factor_base_data[id_thread].last_prime_factor_base=start;
    create_factor_base_f(&(thread_factor_base_data[id_thread].cardinality_factor_base),end,&thread_factor_base_data[id_thread].head,&thread_factor_base_data[id_thread].tail,n,&(thread_factor_base_data[id_thread].last_prime_factor_base));
	return 0;
}
int thread_job_criv_quad(int id_thread){//id inizia da 0,il lavoro di un thread rimane uguale anche se non si riesce a fattorizzare n
	if(id_thread+1>num_thread_job-1){//l'indice del thread eccede il numero di job da fare
		return 0;
	}
    struct timespec timer_thread;//istante di tempo
	int count=id_thread;//indica quale polinomio deve usare per fare il crivello quadratico
    struct node_square_relation*head_squares=NULL,*tail_squares=NULL;
    struct node_square_relation*head_residuoss=NULL,*tail_residuoss=NULL;

    //gettime
    gettime(&timer_thread);
    thread_polynomial_data[id_thread].log_thresold=calculate_log_thresold(n,M);

    //log_thresold
    printf("log_thresold=%f\n",thread_polynomial_data[id_thread].log_thresold);
    print_time_elapsed_local("time to calculate log thresold",&timer_thread);
	while(count<=num_thread_job-2){//ogni thread prende un sottoinsieme di compiti,il thread con id 0 farà i compiti 0,NUM_THREAD,2*NUM_THREAD,il thread 1 farà 1,NUM_THREAD+1,2*NUM_THREAD+1 ecc
		//fattorizzazione,alla fine ogni thread ha una lista di relazioni quadratiche
		printf("thread=%d\n",count);

		//fattorizza array di 2m+1 elementi e memorizza la somma dei logaritmi per ogni posizione e
        // indici last e first che ci dicono il primo elemento divisibile per num e l'ultimo(questo facilita la trial division)
        mpz_set(thread_polynomial_data[id_thread].b,array_bi[count]);//imposta ad ogni ciclo il valore di b
		factor_matrix_f(n,M,(thread_polynomial_data[id_thread]),cardinality_factor_base,a_old,array_a_struct,s);//fattorizza una nuova matrice
        print_time_elapsed_local("time to factor matrix_factorization",&timer_thread);
        //print_thread_data(thread_polynomial_data[id_thread],M);
		//
		//ricerca dei B_smooth potenziali,reali e fattorizzazione dei B_smooth reali
        find_list_square_relation(thread_polynomial_data[id_thread],&(thread_polynomial_data[id_thread].num_B_smooth),&(thread_polynomial_data[id_thread].num_semi_B_smooth),&(thread_polynomial_data[id_thread].num_potential_B_smooth),M,&head_squares,&tail_squares,&head_residuoss,&tail_residuoss,n,a_old,array_a_struct,s);
		//printf("square\n");
		//print_list_square_relation(head_square,thread_polynomial_data[id_thread].num_B_smooth);
		//printf("residuos\n");
		//print_list_square_relation(head_residuos,thread_polynomial_data[id_thread].num_B_smooth);
        printf("num_potential_B_smooth=%d,num_B_smooth=%d,num_semi_B_smooth=%d\n",thread_polynomial_data[id_thread].num_potential_B_smooth,thread_polynomial_data[id_thread].num_B_smooth,thread_polynomial_data[id_thread].num_semi_B_smooth);
		print_time_elapsed_local("time to find_list_square_relation",&timer_thread);
		//pulisci struttura dati del thread per ricominciare con un altro polinomio
		clear_struct_thread_data(thread_polynomial_data[id_thread],M);
		print_time_elapsed_local("time to clear struct thread_data",&timer_thread);
		//unisci la lista dei quadrati trovata con il polinomio con la lista dei quadrati del thread,alla fine ogni thread ha un unica lista dei quadrati
        union_list_square(&(thread_polynomial_data[id_thread].head_square),&(thread_polynomial_data[id_thread].tail_square),head_squares,tail_squares);
        union_list_residuos(&(thread_polynomial_data[id_thread].head_residuos),&(thread_polynomial_data[id_thread].tail_residuos),head_residuoss,tail_residuoss);
		head_squares=NULL;//resetta la lista locale delle relazioni quadratiche
		tail_squares=NULL;//resetta la lista locale delle relazioni quadratiche
        head_residuoss=NULL;
        tail_residuoss=NULL;
        print_time_elapsed_local("time to union list square and list residuos",&timer_thread);
		count+=NUM_THREAD_POLYNOMIAL;//modulo numero dei thread
	}
	return 0;
}
	
