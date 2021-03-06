
#include "main.h"
	extern unsigned long num_times_malloc_called;
	extern unsigned long num_times_free_called;
	//valori globali(presi da altri file)
	extern struct timespec timer;//istante di tempo 
	extern struct timespec time_start;//istante di tempo iniziale
	extern struct timespec timer_test;

	//valori globali del file main
	extern struct row_factorization r;//contiene tutti i primi della factor base e i relativi log
	extern int k;//moltiplicatore di n,se n dispari k=1,3,5 o 7,se n pari so fattorizzarlo
	//argv[1]=path del file da leggere per ottenere n da fattorizzare
	struct node_factor_base*head_f_base_f=NULL;//testa della lista dinamica factor base
	struct node_factor_base*tail_f_base_f=NULL;//coda della lista dinamica factor base
    extern mpz_t n,x0;//dichiarazione di n,n da fattorizzare,deve essere inizializzato a zero,e deve essere sovrascritto con il numero preso da riga 		di comando o da file
	extern mpz_t a_old,a_new;//valore del coefficiente a del polinomio
	mpfr_t thresold_a;//soglia per il calcolo di a
	mpz_t temp;//mpz temporaneo
	extern int s;//numero di primi della factor base distinti che compongono a
	extern long B;//smoothness_bound
	extern long M;//metà dimensione array
	extern long num_elem_array_number;
	mpz_t b1;//valore del primo b=somma di tutti i Bk presi dall'array Bk
	extern mpz_t *array_bi;//array dei coefficienti bi
	mpz_t *array_Bk=NULL;//array che serve per calcolare array_bi
	int dim_sol=-1;//numero di vettori linearmente indpendenti ottenuti dalla risoluzione del sistema lineare,combinandoli opportunamente 		calcolano tutte le possibili combinazioni/soluzioni del sistema
	int cardinality_factor_base=-1;//cardinalità factor base
	extern struct thread_data*thread_polynomial_data;//struttura dati globale per far svolgere la computazione ai thread
	extern struct factor_base_data*thread_factor_base_data;
	extern int num_thread_job;//lunghezza dell'array di matrici,ogni thread riceve una matrice per fare la computazione
	extern int num_increment_M_and_B;
	int*index_prime_a=NULL;//indice dei primi usati per ottenere a,rispetto alla factor base
	int*number_prime_a=NULL;//numeri primi usati per ottenere a
	extern struct a_struct*array_a_struct;
    int num_combined_relations;
    char combined;//TODO farlo diventare bool
    extern mpz_t thresold_large_prime;
    double thresold_relation;

void print_git_commit_hash(){
    int res=system("echo commit hash: `git log --pretty=format:'%H' -n 1`");
    (void)res;
}

void print_macros_enabled(){
	#if DEBUG == 1
    printf("\t- DEBUG mode enabled.\n");
	#endif
	#if LINEAR_PINNING== 1
	printf("\t- LINEAR_PINNING mode enabled.\n");
	#endif
}

void print_statistics(){
	printf("TODO statistics to print\n");
}

int main(int argc,char*argv[]){
    if(argc!=2 && argc!=3){//se non c'è esattamente un parametro,termina,serve il path in cui c'è scritto il numero all'interno
        handle_error_with_exit("usage<path>\n");
    }
    printf("taken %d parameters\n",argc);
    print_git_commit_hash();
    if(argc == 3){//pass compilation line as argument
        const char* compilation_line="make";
        if(strncmp(argv[2],compilation_line,strlen(compilation_line))==0){
            printf("Compilation line: %s\n",argv[2]);
        }
        else{
            printf("Last parameter must be compilation_line and start with \"%s\"\n",compilation_line);
            exit(EXIT_FAILURE);
        }    
    }
    print_macros_enabled();
    set_max_memory_allocable(MAX_MEMORY_ALLOCABLE);
    #if DEBUG==1
    test_memory_limit();
    #endif

    #if LINEAR_PINNING!=1
    set_numa_topology();
    #endif

    cpu_set_t oldset;
    if (pin_thread_to_core(0,&oldset))
    {
        handle_error_with_exit("impossible pinning thread to core\n");
    }
	srand((unsigned int)time(NULL));//imposta seme casuale
	check_variable_in_defines();
	FILE*file_number=open_file(argv[1]);//apri file in cui risiede il numero n da fattorizzare
	for(int i=0;i<NUM_OF_N_TO_FACTORIZE;i++){

		//dichiarazione variabili
		num_increment_M_and_B=0;
		cardinality_factor_base=0;//cardinalità della factor base impostata a 0
		char main_thread_work=1;
		combined=0;
		num_combined_relations=0;
		int removed;
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
                //valore speciale ==1,se last_prime==1 allora si aggiungono alla factor base anche -1 e 2
        char factor_base_already_exist=0;//indica se la factor base è già esistente oppure no
        int number_cycle=0;//numero di cicli necessari per trovare la fattorizzazione
		int start=0;//start della factor base
		k=1;//moltiplicatore di n
		unsigned int num_spawned_threads=0,num_times_spawned_threads=0;//indica il numero di volte che i thread vengono creati
		
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
		
		//mpz_set
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
		calculate_best_M_and_B(n,digit,&M,&B,&num_elem_array_number);
		printf("M=%ld,num_elem_array_number=%ld\n",M,num_elem_array_number);
		print_time_elapsed("time to calculate M");
		if(mpz_cmp_si(n,B)<=0){// se n è minore o uguale a B imposta B=n-3
			B=mpz_get_si(n)-2;//b non deve essere maggiore di n-2
		}
		printf("B=%ld\n",B);
		print_time_elapsed("time to calculate B");
        #if DEBUG==1
        if(B<=0 || M<=0){
            handle_error_with_exit("error in B or M\n");
        }
        #endif

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

		print_total_time_elapsed("time to init main variables",time_start);
		while(factorizations_founded<=0){//finquando non sono stati trovati fattori:
		    gettime(&timer);
		    // calcola a=p1*p2*...ps,unisci le relazioni quadratiche vedi se puoi calcolare il sistema lineare,
            // trova souzioni sistema lineare e trova tutti gli a,b del crivello quadratico
			adjust_n(n,&k);//aggiusta n per calcolare i suoi fattori(divide n per k)
			multiply_n_for_k(n,&k,&factorized);//k,moltiplicatore di n per renderlo un quadrato modulo 8

			//thresold_a per applicare siqs
			calculate_thresold_a(thresold_a,n,M);//a è circa rad(2*n)/M
			print_time_elapsed("time to calculate thresold_a");

			//factor base
			if(factor_base_already_exist==0 && B>THRESOLD_B && NUM_THREAD_FACTOR_BASE>0) {//se la factor base non è mai stata creata
				// e se B è maggiore del valore soglia e numero thread per creare factor base>0
				thread_factor_base_data = alloc_array_factor_base_data(NUM_THREAD_FACTOR_BASE);//alloca struttura condivisa a tutti i thread,
				        // ogni thread ha il proprio indice
				array_tid = alloc_array_tid(NUM_THREAD_FACTOR_BASE);//alloca memoria per contenere tutti i tid
				//TODO array_id non è necessario, non va allocato e non va liberato
				array_id = create_factor_base_threads(array_tid, NUM_THREAD_FACTOR_BASE);//crea tutti i thread e ogni thread crea la sua porzione di factor base
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

			//verifica che la factor base è corretta
			#if DEBUG==1
            if(verify_factor_base(head_f_base_f,cardinality_factor_base,last_prime_factor_base)==0){
                handle_error_with_exit("error in main verify factor base\n");
            }
            if(verify_cardinality_list_factor_base(head_f_base_f,cardinality_factor_base)==0){
            	handle_error_with_exit("invalid dimension list factor base\n");
            }
            #endif
            print_time_elapsed("time to calculate factor base");
			//a,per siqs e generare tutti gli altri b,prodotto di primi dispari distinti
			calculate_a(a_new,thresold_a,&s,head_f_base_f,cardinality_factor_base,&index_prime_a,&number_prime_a);
            #if DEBUG==1
            if(s==0 && mpz_cmp_si(a_new,0)!=0){
                handle_error_with_exit("error in main calculate a 1\n");
            }
            if(s>0 && mpz_cmp_si(a_new,0)==0){
                handle_error_with_exit("error in main calculate a 2\n");
            }
            #endif
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
				calculate_a(a_new,thresold_a,&s,head_f_base_f,cardinality_factor_base,&index_prime_a,&number_prime_a);
                #if DEBUG==1
                if(s==0 && mpz_cmp_si(a_new,0)!=0){
                    handle_error_with_exit("error in main calculate a 1\n");
                }
                if(s>0 && mpz_cmp_si(a_new,0)==0){
                    handle_error_with_exit("error in main calculate a 2\n");
                }
                #endif
            }
            mpz_set(a_old,a_new);//imposta a_old=a_new
			if(s==0 && mpz_cmp_si(a_new,0)==0){
				increment_M_and_B(&M,&B);//aumenta M e B
				create_factor_base_f(&cardinality_factor_base,B,&head_f_base_f,&tail_f_base_f,n,&last_prime_factor_base);
				#if DEBUG==1
				if(verify_factor_base(head_f_base_f,cardinality_factor_base,last_prime_factor_base)==0){
					handle_error_with_exit("error in main verify factor base\n");
				}
				if(verify_cardinality_list_factor_base(head_f_base_f,cardinality_factor_base)==0){
					handle_error_with_exit("invalid dimension list factor base\n");
				}
				#endif
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
            }
            print_time_elapsed("time to calculate a");
			//array Bk,può essere ridotto modulo a
			array_Bk=calculate_array_Bk(number_prime_a,cardinality_factor_base,n,s,a_old,b1);

			//array bi
			array_bi=calculate_bi(array_Bk,b1,s);
			if(array_Bk!=NULL){
				free_memory_array_mpz(array_Bk,s);
				array_Bk=NULL;
			}
			if(array_bi!=NULL){//bi può essere ridotto modulo a
				reduce_array_mpz_mod_n(array_bi,(int)pow(2,s-1),a_old);
				adjust_array_bi(array_bi,s,a_old);
			}

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
			print_time_elapsed("time to calculate array_bi and array_Bk");
			//alloca struttura dati dei thread per fare il sieving
            thread_polynomial_data=alloc_array_polynomial_thread_data(NUM_THREAD_POLYNOMIAL+1,M);//TODO la struttura dati in più serve per aggregare?
			print_time_elapsed("time to alloc polynomial data");
			//creazione della struttura row_factorization
			// row_factorization contiene primi factor base,log per ogni primo radice 1 e radice 2 di n mod p e a^-1 mod p)
			//TODO invece di allocare e deallocare e ricalcolare row_factorization è possibile farlo ogni volta anche non in modo parallelo ma in modo append
			r.prime=alloc_array_int(cardinality_factor_base);
			r.log_prime=alloc_array_int(cardinality_factor_base);
			r.root_n_mod_p=alloc_array_int(cardinality_factor_base);
			r.root2_n_mod_p=alloc_array_int(cardinality_factor_base);
			r.inverse_a_mod_p=alloc_array_int(cardinality_factor_base);
			create_row_factorization(head_f_base_f,cardinality_factor_base,a_old,array_a_struct,s);
			print_time_elapsed("time to create_row_factorization");
			//thresold_large_prime,non è costoso e va ricalcolato ogni volta
			calculate_thresold_large_prime(thresold_large_prime,r.prime[cardinality_factor_base-1]);

			//creazione e avvio thread per fase sieving
			if(num_thread_job!=1){
				array_tid=alloc_array_tid(NUM_THREAD_POLYNOMIAL);//alloca memoria per contenere tutti i tid
				array_id=create_threads(array_tid,NUM_THREAD_POLYNOMIAL);//crea tutti i thread
				num_spawned_threads+=NUM_THREAD_POLYNOMIAL;
				num_times_spawned_threads++;
			}
			//fattorizza numeri nell'array lungo 2m+1
            if(s==0){//se non ci sono thread disponibili riattiva il main thread
                main_thread_work=1;
            }
			if(main_thread_work==1) {
				if(num_thread_job!=1) {//se il numero di job da fare è maggiore di 1
					main_thread_work = 0;
				}
				mpz_set(thread_polynomial_data[NUM_THREAD_POLYNOMIAL].b, b_default);//imposta b
				factor_matrix(n, M, (&(thread_polynomial_data[NUM_THREAD_POLYNOMIAL])), cardinality_factor_base,
								a_default, array_a_struct, s);//fattorizza numeri

				//trova relazioni quadratiche o semi_B_smooth e ordinale per numero
				find_list_square_relation(&(thread_polynomial_data[NUM_THREAD_POLYNOMIAL]), &num_B_smooth,
										  &num_semi_B_smooth, &num_potential_B_smooth, M, &head_square, &tail_square,
										  &head_residuos, &tail_residuos, n, a_default, NULL, 0);
			}
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

			//verifica
			#if DEBUG==1
			if(verify_cardinality_list_square_relation(head_square,num_B_smooth)==0){
				handle_error_with_exit("error in cardinality head_square first of union\n");
			}
			#endif
			for(int i=0;i<NUM_THREAD_POLYNOMIAL;i++){//metti tutte le relazioni in head e tail,somma tutti i numeri B_smooth e semi_B_smooth
				//le liste vengono unite in modo tale che quella finale è ordinata per numero
			    union_list_square(&head_square,&tail_square,
								  thread_polynomial_data[i].head_square,thread_polynomial_data[i].tail_square);
				thread_polynomial_data[i].head_square=NULL;
				thread_polynomial_data[i].tail_square=NULL;
                union_list_residuos(&head_residuos,&tail_residuos,
                                  thread_polynomial_data[i].head_residuos,thread_polynomial_data[i].tail_residuos);
				thread_polynomial_data[i].head_residuos=NULL;
				thread_polynomial_data[i].tail_residuos=NULL;
				num_B_smooth+=thread_polynomial_data[i].num_B_smooth;
				thread_polynomial_data[i].num_B_smooth=0;
				num_potential_B_smooth+=thread_polynomial_data[i].num_potential_B_smooth;
				thread_polynomial_data[i].num_potential_B_smooth=0;
				num_semi_B_smooth+=thread_polynomial_data[i].num_semi_B_smooth;
				thread_polynomial_data[i].num_semi_B_smooth=0;
			}
            printf("num_potential_B_smooth=%d,num_B_smooth=%d,num_semi_B_smooth=%d\n",num_potential_B_smooth,num_B_smooth,num_semi_B_smooth);
            printf("card_f_base=%d,B=%ld,M=%ld\n",cardinality_factor_base,B,M);
			#if DEBUG==1
			if(verify_cardinality_list_square_relation(head_square,num_B_smooth)==0){
				handle_error_with_exit("error in cardinality head_square after union\n");
			}
			#endif
            //trova nuove relazioni quadratiche con un nuovo square,una nuova fattorizazzione e imposta num=0
            if((num_B_smooth>=cardinality_factor_base*thresold_relation) && (combined==0)) {
            	combined=1;
                quickSort_residuos(head_residuos);
                head_sort_residuos=head_residuos;
                tail_sort_residuos=lastNode(head_sort_residuos);
                head_residuos=NULL;
                tail_residuos=NULL;
                #if DEBUG==1
                if(verify_sorted_residuos_square_rel_list(head_sort_residuos)==0){
                    handle_error_with_exit("error in sort relation by square\n");
                }
                #endif
				factorizations_founded = combine_relation_B_smooth_and_semi_B_smooth_v3(&head_square,
																						&tail_square, &head_sort_residuos,&tail_sort_residuos, n, &num_B_smooth, &num_semi_B_smooth,&num_combined_relations);
                printf("num_combined_relations=%d\n",num_combined_relations);
				head_sort_residuos = NULL;
				tail_sort_residuos = NULL;
                if (factorizations_founded == 1) {
                    break;
                }
			}
			else if(combined==1 && head_residuos!=NULL){
				free_memory_list_square_relation(head_residuos);
				head_residuos=NULL;
				tail_residuos=NULL;
				head_sort_residuos = NULL;
				tail_sort_residuos = NULL;
            }
            if(num_B_smooth>=cardinality_factor_base*ENOUGH_RELATION){
				#if DEBUG==1
				if(verify_cardinality_list_square_relation(head_square,num_B_smooth)==0){
					handle_error_with_exit("error in cardinality head_square\n");
				}
				#endif
				quickSort_square(head_square);
				head_sort_square=head_square;
				tail_sort_square=lastNode(head_sort_square);

				//head_square e tail_square rimangono !=null fino alla fine
                #if DEBUG==1
                if(verify_sorted_square_rel_list(head_sort_square)==0){
                    handle_error_with_exit("error in sort relation by square\n");
                }
                if(verify_cardinality_list_square_relation(head_sort_square,num_B_smooth)==0){
                    handle_error_with_exit("error in cardinality head_sort_square\n");
                }
                #endif
                removed=remove_same_square(&head_sort_square,&tail_sort_square,&num_B_smooth,&num_semi_B_smooth);
                printf("num removed relations=%d,num_B_smooth=%d\n",removed,num_B_smooth);
                #if DEBUG==1
                if(verify_sorted_square_rel_list(head_sort_square)==0){
                    handle_error_with_exit("error in sort list by square\n");
                }
				if(verify_cardinality_list_square_relation(head_sort_square,num_B_smooth)==0){
					handle_error_with_exit("error in cardinality square relation\n");
				}
				#endif
				head_square=head_sort_square;
                tail_square=tail_sort_square;
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
			print_total_time_elapsed("time total to finish math steps",time_start);
			printf("num spawned threads=%u,num_times_spawned_threads=%u\n",num_spawned_threads,num_times_spawned_threads);
			printf("num malloc called=%ld\n",num_times_malloc_called);
			printf("num free called=%ld\n",num_times_free_called);
			//handle_error_with_exit("finish math steps\n");







			//algebra step:sistema lineare
			gettime(&timer);//init timer
            binary_linear_system=create_binary_linear_system(head_sort_square,cardinality_factor_base,num_B_smooth,&num_col_binary_matrix);
            print_time_elapsed("time_to_create_binary linear_system");
            reduce_echelon_form_binary_matrix(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            print_time_elapsed("time_to_reduce echelon form linear_system");
            linear_system=from_matrix_binary_to_matrix_char(binary_linear_system,cardinality_factor_base,num_col_binary_matrix,&num_col_linear_system);
            print_time_elapsed("time_to_copy matrix binary into matrix char");
            free_memory_matrix_unsigned_long(binary_linear_system,cardinality_factor_base,num_col_binary_matrix);
            binary_linear_system=NULL;

            //algebra step:base sistema lineare
            base_matrix=calculate_base_linear_system_char(linear_system,cardinality_factor_base,num_col_linear_system,&dim_sol,MAX_DIM_SOL);
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
			#if DEBUG==1
			if(check_if_matrix_is_reduce_mod_n(base_matrix,num_B_smooth,dim_sol,2)==0){
				handle_error_with_exit("error in main,calculate base_linear_system\n");
			}
			#endif
			print_time_elapsed("time to calculate base linear system");
			printf("num_combined_relations=%d\n",num_combined_relations);
			#if DEBUG==1
			if(check_solution_base_matrix_char(linear_system,cardinality_factor_base,num_col_linear_system,
            base_matrix,num_B_smooth,dim_sol)==0){
				handle_error_with_exit("error in main,invalid solution\n");
			}
			#endif
			print_time_elapsed("time to check solution base linear system");

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
		}//fine while factorization_founded>=0




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
			head_square=NULL;
		}
		if(head_square!=NULL) {
			free_memory_list_square_relation(head_square);
			head_square = NULL;
			head_sort_square=NULL;
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
		printf("num_combined_relations=%d,combined=%d\n",num_combined_relations,combined);
		timer.tv_nsec=time_start.tv_nsec;//timer=time_start
		timer.tv_sec=time_start.tv_sec;//timer=time_start
		print_time_elapsed("time_total");
        if(fclose(file_number)!=0){
            handle_error_with_exit("error in close file_number\n");
        }
	}
	printf("num malloc called=%ld\n",num_times_malloc_called);
	printf("num free called=%ld\n",num_times_free_called);
	print_statistics();
	printf("Max allocated space...........................: %lu MB\n", get_memory_allocated()/(1024));
	return 0;
}