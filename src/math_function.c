#include "math_function.h"


extern struct row_factorization r;
extern long M;
extern mpz_t n;
extern mpz_t a_old;
extern struct a_struct*array_a_struct;
extern int cardinality_factor_base;
extern int s;
extern mpz_t thresold_large_prime;












































void calculate_thresold_large_prime(mpz_t thresold_large_prime,int max_prime){
	if(max_prime<=0){
		handle_error_with_exit("error in calculate_thresold_prime\n");
	}
	mpz_t temp;
	mpz_init(temp);
	mpz_set_si(temp,max_prime);//temp=max_prime
	mpz_mul(temp,temp,temp);//temp?max_prime^2
	mpz_set(thresold_large_prime,temp);
	mpz_clear(temp);
	return;
}
struct a_struct*create_array_a_struct(int*number_prime_a,int*index_number_a,int length){
    if(number_prime_a==NULL || index_number_a==NULL || length<=0){
        handle_error_with_exit("error in create_array_a_struct\n");
    }
    struct a_struct*array_a_struct=malloc(sizeof(struct a_struct)*length);
    if(array_a_struct==NULL){
        handle_error_with_exit("error in malloc create array a struct\n");
    }
    for(int i=0;i<length;i++){
        array_a_struct[i].number_prime_a=number_prime_a[i];
        array_a_struct[i].index_prime_a=index_number_a[i];
    }
    return array_a_struct;
}

void calculate_index_min_max_a(int*number_prime_a,int*index_prime_a,int length,int*index_min_a,int*index_max_a){
	if(number_prime_a==NULL || index_max_a==NULL || length<=0 || index_min_a==NULL || index_max_a==NULL){
		handle_error_with_exit("error in calculate min_max_a\n");
	}
	int min_a=number_prime_a[0];
	*index_min_a=index_prime_a[0];
	int max_a=number_prime_a[0];
	*index_max_a=index_prime_a[0];
	for(int i=1;i<length;i++){
		if(min_a>number_prime_a[i]){
			min_a=number_prime_a[i];
			*index_min_a=index_prime_a[i];
		}
		if(max_a<number_prime_a[i]) {
			max_a=number_prime_a[i];
			*index_max_a=index_prime_a[i];
		}
	}
	return;
}
float calculate_log_thresold(const mpz_t n,long M){
	//log_thresold=log(M*rad(n))-error=log(M)+log(rad(n)-error=log(M)+1/2*log(n)-error
    float log_thresold=-1.0;
    mpfr_t log_n;
    mpfr_init(log_n);
    mpfr_set_z(log_n,n,MPFR_RNDN);//log_rad_n=n

	//mpfr_sqrt(log_rad_n,log_rad_n,MPFR_RNDN);//log_rad_n=radice di n
	mpfr_log2(log_n,log_n,MPFR_RNDN);//rad_n=log2(n)
	log_thresold=mpfr_get_flt(log_n,MPFR_RNDN);
	log_thresold=log_thresold*0.5;//log_thresold=log2(n)
    log_thresold=log_thresold+log2f((float)M);
    log_thresold=log_thresold-ERROR_LOG;
    mpfr_free_cache();
    mpfr_clear(log_n);
    return log_thresold;
}
void create_num(mpz_t num,const mpz_t a,const mpz_t b,const mpz_t n,long j){
	//num=a*j^2+2*bj+c mod n
	mpz_t c,double_b_mul_j,a_mul_square_j;

	mpz_init(c);
	mpz_init(double_b_mul_j);
	mpz_init(a_mul_square_j);

	//a_mul_square_j=a
	mpz_set(a_mul_square_j,a);//a_mul_suqare_j=a
	mpz_mul_si(a_mul_square_j,a_mul_square_j,j);//a_mul_suqare_j=a*j
	mpz_mul_si(a_mul_square_j,a_mul_square_j,j);//a_mul_suqare_j=a*j^2

	//double_b_mul_j
	mpz_set(double_b_mul_j,b);//double_b_mul_j=b
	mpz_mul_2exp(double_b_mul_j,double_b_mul_j,1);//double_b_mul_j=2*b
	mpz_mul_si(double_b_mul_j,double_b_mul_j,j);//double_b_mul_j=2*b*j

	//c
	mpz_mul(c,b,b);//c=b^2
	mpz_sub(c,c,n);//c=b^2-n
	if(mpz_divisible_p(c,a)==0){
		handle_error_with_exit("error in create_num\n");
	}
	mpz_divexact(c,c,a);//c=b^2-n/a

	mpz_add(num,a_mul_square_j,double_b_mul_j);
	mpz_add(num,num,c);//num=a*j^2+2*bj+c

	mpz_clear(c);
	mpz_clear(double_b_mul_j);
	mpz_clear(a_mul_square_j);
}
struct node_factorization*factorize_num_v2(const mpz_t num,int j_of_num,int first_index_f_base,int last_index_f_base,
        char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,struct a_struct*array_a_struct,int s,struct thread_data thread_data){
    if((first_index_f_base<0 && first_index_f_base!=-1)|| (last_index_f_base<0 && last_index_f_base!=-1)
        || is_B_smooth==NULL || is_semi_B_smooth==NULL
            || first_index_f_base>last_index_f_base || (array_a_struct==NULL && s>0)){
		handle_error_with_exit("error in factorize_num\n");
	}
    //j_of_num va da -M a +M
	int prime=0;
	int exp=0;
	int index=0;
	long j1,j2,j_of_num_mod_p;
	mpz_t temp;
	mpz_init(temp);
	mpz_set(temp,num);
	struct node_factorization*head=NULL;
	struct node_factorization*tail=NULL;
	//bisogna aggiungere i fattori di a alla fattorizzazione

	if(mpz_cmp_si(num,0)<0){//valore negativo->divisibile per 1
		mpz_neg(temp,temp);//rendilo positivo
		insert_ordered_factor(-1,1,0,&head,&tail);//inserisci nodo -1
	}
	if(s!=0) {
        //se i numeri di a sono ad un indice più piccolo del first_index_f_base allora aggiungili prima
        //questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
        while (index!=s && array_a_struct[index].index_prime_a < first_index_f_base) {
            insert_ordered_factor(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
                                  &tail);//inserisci nodo con esponente 1
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
            index++;
        }
    }
	exp=0;
	for(int i=first_index_f_base;i<=last_index_f_base;i++){//scorri i primi da fist index a last index compresi
		if(exp!=0){
			handle_error_with_exit("error in factorize num v2,invalid exponent\n");
		}
		if(first_index_f_base==-1 || last_index_f_base==-1){//nessun fattore trovato,
			break;//aggiungi solamente i fattori di a
		}
		if(r.prime[i]==-1){//salta il -1,l'abbiamo già considerato
	        continue;
	    }
	    if(s!=0 && index!=s && array_a_struct[index].index_prime_a==i){//se un fattore di a ha un indice compreso
	    	// tra first e last index metti ad uno l'esponente anche se non è divisibile per quel numero
			// (sappiamo che il polinomio è divisibile per quel numero)
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
	    	exp+=1;//poni esponente a 1
	    	index++;//aumenta l'indice
            //n.b. i fattori di a possono comparire nella fattorizzazione con un esponente maggiore di 1
	    }
		prime=r.prime[i];
	    if(i==first_index_f_base || i==last_index_f_base){//se è uguale al primo o all'ultimo indice è sicuramente divisibile per p
			if(i==1){//divisibile per 2
                mpz_divexact_ui(temp, temp,prime);//dividilo la prima volta per prime
                exp+=1;
                while(mpz_divisible_2exp_p(temp,1)!=0){//finquando è divisibie per 2
                    mpz_divexact_ui(temp, temp,prime);//dividilo
                    exp+=1;//aumenta esponente
                }
                if(exp>0) {//è stato diviso almeno 1 volta
                    insert_ordered_factor(prime, exp, i, &head, &tail);
                }
                exp=0;//resetta esponente
                continue;
			}
			else {
                mpz_divexact_ui(temp, temp, prime);//dividilo la prima volta per prime
                exp += 1;
                while (mpz_divisible_ui_p(temp, prime) != 0) {//se il numero è divisibile per un primo della fattor base
                    mpz_divexact_ui(temp, temp, prime);//dividilo
                    exp += 1;//aumenta esponente
                }
                if (exp > 0) {//è stato diviso almeno 1 volta
                    insert_ordered_factor(prime, exp, i, &head, &tail);
                }
                exp = 0;//resetta esponente
                continue;
            }
	    }
	    j1=thread_data.j1_mod_p[i];
	    j2=thread_data.j2_mod_p[i];
        j_of_num_mod_p=reduce_int_mod_n_v2(j_of_num,prime);
        //j_of_num ridotto modulo prime se è uguale a j1 o j2 allora è divisibile per p
	    if(j_of_num_mod_p==j1 || j_of_num_mod_p==j2){
            mpz_divexact_ui(temp, temp,prime);//dividilo la prima volta per prime
            exp+=1;//aumenta esponente
            //verifica se è ulteriormente divisibile per prime facendo la divisione
            while(mpz_divisible_ui_p(temp,prime)!=0) {//se il numero è divisibile per un primo della fattor base
                mpz_divexact_ui(temp, temp,prime);//dividilo
                exp+=1;//aumenta esponente
            }
			if(exp>0) {//è stato diviso almeno 1 volta
				insert_ordered_factor(prime, exp, i, &head, &tail);
			}
			exp=0;//resetta esponente
			continue;
	    }
		if(exp>0) {//è stato diviso almeno 1 volta
			insert_ordered_factor(prime, exp, i, &head, &tail);
		}
		exp=0;//resetta esponente
		continue;
    }
    if(s!=0) {
        //se i numeri di a sono ad un indice più grande del last_index_f_base allora aggiungili dopo
		//questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
        while (index!=s && array_a_struct[index].index_prime_a > last_index_f_base && index < s) {
            insert_ordered_factor(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
                                  &tail);//inserisci nodo
            if(thread_data.j2_mod_p[array_a_struct[index].index_prime_a]!=-1){
                handle_error_with_exit("error in factorize num\n");
            }
            index++;
        }
    }
	if(s!=0 && index!=s){
		handle_error_with_exit("error in factorize num invalid s\n");
	}
	if(head==NULL){
		mpz_set(residuos, temp);//imposta il residuo
		mpz_clear(temp);
		*is_B_smooth=0;
		*is_semi_B_smooth=0;
		return NULL;
	}
    if(mpz_cmp_si(temp,1)==0){//se il residuo della divisione è 1 allora è B-smooth
	    *is_B_smooth=1;
	    *is_semi_B_smooth=0;
	    mpz_set(residuos,temp);//imposta il residuo
	}
	else{//non è B-smooth è semi_B_smooth
    	*is_B_smooth=0;
    	*is_semi_B_smooth=1;
		mpz_set(residuos,temp);//imposta il residuo
    }
    mpz_clear(temp);
	return head;
}
struct node_factorization*factorize_num_v1(const mpz_t num,int first_index_f_base,int last_index_f_base,
										char*is_B_smooth,char*is_semi_B_smooth,mpz_t residuos,struct a_struct*array_a_struct,int s) {
	if ((first_index_f_base < 0 && first_index_f_base != -1) || (last_index_f_base < 0 && last_index_f_base != -1)
		|| is_B_smooth == NULL || is_semi_B_smooth == NULL
		|| first_index_f_base > last_index_f_base || (array_a_struct == NULL && s > 0)) {
		handle_error_with_exit("error in factorize_num\n");
	}
	//j_of_num va da -M a +M
	int prime = 0;
	int exp = 0;
	int index = 0;
	mpz_t temp;
	mpz_init(temp);
	mpz_set(temp, num);
	struct node_factorization *head = NULL;
	struct node_factorization *tail = NULL;
	//bisogna aggiungere i fattori di a alla fattorizzazione

	if (mpz_cmp_si(num, 0) < 0) {//valore negativo->divisibile per 1
		mpz_neg(temp, temp);//rendilo positivo
		insert_ordered_factor(-1, 1, 0, &head, &tail);//inserisci nodo -1
	}
	if (s != 0) {
		//se i numeri di a sono ad un indice più piccolo del first_index_f_base allora aggiungili prima
		//questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
		while (index != s && array_a_struct[index].index_prime_a < first_index_f_base) {
			insert_ordered_factor(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
								  &tail);//inserisci nodo con esponente 1
			index++;
		}
	}
	for (int i = first_index_f_base; i <= last_index_f_base; i++) {//scorri i primi da fist index a last index compresi
		if (first_index_f_base == -1 || last_index_f_base == -1) {//nessun fattore trovato,aggiungi solo i fattori di a
			break;
		}
		if (r.prime[i] == -1) {//salta il -1,l'abbiamo già considerato
			continue;
		}
		if (s != 0 && index != s &&
			array_a_struct[index].index_prime_a == i) {//se un fattore di a ha un indice compreso
			// tra first e last index metti ad uno l'esponente anche se non è divisibile per quel numero
			// (sappiamo che il polinomio è divisibile per quel numero)
			exp = 1;//poni esponente a 1
			index++;//aumneta l'indice
			//n.b. i fattori di a possono comparire nella fattorizzazione con un esponente maggiore di 1
		}
		prime = r.prime[i];
		while (mpz_divisible_ui_p(temp, prime) != 0) {//se il numero è divisibile per un primo della fattor base
			mpz_divexact_ui(temp, temp, prime);//dividilo
			exp += 1;//aumenta esponente
		}
		if (exp > 0) {//è stato diviso almeno 1 volta
			insert_ordered_factor(prime, exp, i, &head, &tail);
		}
		exp = 0;//resetta esponente
		continue;
	}
	if (s != 0) {
		//se i numeri di a sono ad un indice più grande del last_index_f_base allora aggiungili dopo
		//questi fattori di a non dividono il polinomio a*j^2+2bj+c ->avranno sempre esponente 1
		while (index != s && array_a_struct[index].index_prime_a > last_index_f_base && index < s) {
			insert_ordered_factor(array_a_struct[index].number_prime_a, 1, array_a_struct[index].index_prime_a, &head,
								  &tail);//inserisci nodo
			index++;
		}
	}
	if (s != 0 && index != s) {
		handle_error_with_exit("error in factorize num invalid s\n");
	}
	if(head==NULL){
		mpz_set(residuos, temp);//imposta il residuo
		mpz_clear(temp);
		return NULL;
	}
	if (mpz_cmp_si(temp, 1) == 0) {//se il residuo della divisione è 1 allora è B-smooth
		*is_B_smooth = 1;
		*is_semi_B_smooth = 0;
		mpz_set(residuos, temp);//imposta il residuo
	} else {//non è B-smooth è semi_B_smooth
		*is_B_smooth = 0;
		*is_semi_B_smooth = 1;
		mpz_set(residuos, temp);//imposta il residuo
	}
	mpz_clear(temp);
	return head;
}
void calculate_square(mpz_t square,const mpz_t a,int index,const mpz_t b,const mpz_t n){
    // -M<index<M
    mpz_mul_si(square,a,index);//a*j
    mpz_add(square,square,b);//a*j+b
	mpz_mod(square,square,n);// a*j+b mod n
    return;
}
char calculate_num_from_factorization(mpz_t num_temp,struct node_factorization*head_factor){
	long exp_of_factor,factor;
	mpz_t factor_raise_to_exp;
	mpz_t temp;
	mpz_init(temp);
	mpz_init(factor_raise_to_exp);

	mpz_set_si(temp,1);//temp=1
	if(head_factor==NULL){
		mpz_clear(temp);
		mpz_clear(factor_raise_to_exp);
		printf("fattorizzazione vuota\n");
		return 0;
	}
	struct node_factorization*p=head_factor;
	while(p!=NULL){
		mpz_set_si(factor_raise_to_exp,1);
		factor=p->number;
		//factor
		if(factor<0 && factor!=-1){
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			printf("fattore minore di 0\n");
			return 0;
		}
		if(factor==0){
			printf("fattore 0\n");
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			return 0;
		}
		//exp of factor
		exp_of_factor=p->exp_of_number;
		if(exp_of_factor<=0){
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			printf("esponente minore o uguale a 0\n");
			return 0;
		}
		if(factor==-1) {
			if(exp_of_factor!=1) {
				handle_error_with_exit("exponent of -1 greater than 1\n");
			}
			else {
				mpz_neg(temp,temp);
				p=p->next;
				continue;
			}
		}
		mpz_ui_pow_ui (factor_raise_to_exp,factor,exp_of_factor);
		mpz_mul(temp,temp,factor_raise_to_exp);
		p=p->next;
	}
	mpz_set(num_temp,temp);
	mpz_clear(temp);
	mpz_clear(factor_raise_to_exp);
	return 1;
}
char verify_square_relation(struct square_relation square_relation,const mpz_t n){
    if(TEST==0){
        return 1;
    }
	struct node_factorization*head_factor=square_relation.head_factorization;
	mpz_t temp,num_temp,square;
	mpz_init(num_temp);
	mpz_init(temp);
	mpz_init(square);

	if(calculate_num_from_factorization(num_temp,head_factor)==0){
		handle_error_with_exit("error in calculate_num_from_factorization\n");
	}
	mpz_mod(num_temp,num_temp,n);//num modulo n
	mpz_mul(square,square_relation.square,square_relation.square);//square=x^2
	mpz_mod(square,square,n);
	if(mpz_cmp(num_temp,square)!=0){//non sono congrui modulo n
		handle_error_with_exit("error in verify square relation\n");
	}
	mpz_clear(temp);
	mpz_clear(num_temp);
	mpz_clear(square);
	return 1;
}
char verify_factorization(const mpz_t num,mpz_t residuos,struct node_factorization*head_factor,const mpz_t a){
	//num*a deve essere uguale a fattorizzazione*residuo
	mpz_t temp,num_temp;
	long exp_of_factor,factor;
	mpz_t factor_raise_to_exp;
    if(TEST==0){
        return 1;
    }
	mpz_init(num_temp);
	mpz_init(temp);
	mpz_init(factor_raise_to_exp);
	mpz_set(num_temp,num);
	if(mpz_cmp_si(num_temp,0)<0){
		mpz_neg(num_temp,num_temp);
	}
	mpz_set_si(temp,1);//temp=1
	if(head_factor==NULL){
		if(mpz_cmp(num,residuos)!=0){
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			mpz_clear(num_temp);
			printf("residuo !=num\n");
			return 0;
		}
		mpz_clear(temp);
		mpz_clear(factor_raise_to_exp);
		mpz_clear(num_temp);
		return 1;
	}
	struct node_factorization*p=head_factor;
	while(p!=NULL){
		mpz_set_si(factor_raise_to_exp,1);
		factor=p->number;
		//factor
		if(factor<0 && factor!=-1){
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			mpz_clear(num_temp);
			printf("fattore minore di 0\n");
			return 0;
		}
		if(factor==0){
			printf("fattore 0\n");
            mpz_clear(temp);
            mpz_clear(factor_raise_to_exp);
			mpz_clear(num_temp);
			return 0;
		}
		if(factor==-1){
			factor=1;
		}
		//exp of factor
		exp_of_factor=p->exp_of_number;
		if(exp_of_factor<=0){
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			mpz_clear(num_temp);
			printf("esponente minore o uguale a 0\n");
			return 0;
		}
		mpz_ui_pow_ui (factor_raise_to_exp,factor,exp_of_factor);
		mpz_mul(temp,temp,factor_raise_to_exp);
		p=p->next;
	}
	if(mpz_cmp_si(residuos,0)<=0){
		printf("residuo minore o uguale a 0\n");
        mpz_clear(temp);
        mpz_clear(factor_raise_to_exp);
		mpz_clear(num_temp);
		return 0;
	}
	mpz_mul(temp,temp,residuos);//moltiplica temp per il residuo alla fine
	mpz_mul(num_temp,num_temp,a);//num*a
	if(mpz_cmp_si(temp,0)<=0){
		handle_error_with_exit("error in function verify factorization\n");
	}
	if(mpz_cmp(temp,num_temp)!=0){//se non sono uguali nega temp
			gmp_printf("temp e num diversi,residuos=%Zd,num=%Zd,temp=%Zd,a=%Zd\n",residuos,num,temp,a);
			mpz_clear(temp);
			mpz_clear(factor_raise_to_exp);
			mpz_clear(num_temp);
			print_factorization(num,head_factor);
			return 0;
	}
	mpz_clear(temp);
	mpz_clear(factor_raise_to_exp);
	mpz_clear(num_temp);
	return 1;
}
void find_list_square_relation(struct thread_data thread_data, int *num_B_smooth,int*num_semi_B_smooth, int *num_potential_B_smooth, long M,
							   struct node_square_relation **head_square, struct node_square_relation **tail_square,struct node_square_relation **head_residuos, struct node_square_relation **tail_residuos,
							   const mpz_t n,const mpz_t a,struct a_struct*array_a_struct,int s) {
	mpz_t num;
	struct node_factorization*head_factor=NULL;
	char is_B_smooth=0;
	char is_semi_B_smooth=0;
	mpz_t residuos;
	if(num_B_smooth==NULL || num_semi_B_smooth==NULL || num_potential_B_smooth==NULL || M<=0 ||
	   head_square==NULL || tail_square==NULL || head_residuos==NULL
	   || tail_residuos==NULL || (array_a_struct==NULL && s>0)){
		handle_error_with_exit("error in find_list_square_relation");
	}
	struct square_relation square_relation;
	mpz_init(num);
	mpz_init(residuos);
	for(long i=0;i<2*M+1;i++){
        if(thread_data.numbers[i].sum_log>=thread_data.log_thresold){
            //possibile B_smooth trovato
            (*num_potential_B_smooth)++;
            create_num(num,a,thread_data.b,n,thread_data.numbers[i].j);
            if(mpz_cmp_si(num,0)==0){
                head_factor=NULL;
                mpz_set_si(residuos,0);
				is_B_smooth=0;
				is_semi_B_smooth=0;
                continue;
            }
            else {
                head_factor = factorize_num_v2(num, thread_data.numbers[i].j, thread_data.numbers[i].first_index_f_base,
                                               thread_data.numbers[i].last_index_f_base, &is_B_smooth,
                                               &is_semi_B_smooth, residuos, array_a_struct, s, thread_data);
            }
            if(head_factor==NULL && (is_B_smooth==1 || is_semi_B_smooth==1)){
            	handle_error_with_exit("error invalid factorize_num\n");
            }
			if(verify_factorization(num,residuos,head_factor,a)==0){
            	handle_error_with_exit("error in factorization\n");
            }
            if(head_factor==NULL){
            	continue;
            }
            if(is_B_smooth==1){
                (*num_B_smooth)++;
                is_B_smooth=0;
                is_semi_B_smooth=0;
                square_relation.head_factorization=head_factor;
                mpz_init(square_relation.square);
                mpz_init(square_relation.num);
                mpz_init(square_relation.residuos);
                mpz_set(square_relation.num,num);
                mpz_set(square_relation.residuos,residuos);
                calculate_square(square_relation.square,a,i-M,thread_data.b,n);
                insert_at_tail_square_relation(square_relation,head_square,tail_square);
                if(verify_square_relation(square_relation,n)==0){
                	handle_error_with_exit("error in create square relation\n");
                }
            }
            else if(is_semi_B_smooth==1){
            	if(mpz_cmp(residuos,thresold_large_prime)<=0) {
					(*num_semi_B_smooth)++;
					is_B_smooth = 0;
					is_semi_B_smooth = 0;
					square_relation.head_factorization = head_factor;
					mpz_init(square_relation.square);
					mpz_init(square_relation.num);
					mpz_init(square_relation.residuos);
					mpz_set(square_relation.num, num);
					mpz_set(square_relation.residuos, residuos);
					calculate_square(square_relation.square, a, i - M, thread_data.b, n);
					insert_at_tail_square_relation(square_relation, head_residuos, tail_residuos);
				}
				else{
            		free_memory_list_factor(head_factor);
            	}
            }
            else{
            	handle_error_with_exit("error in find list square relation\n");
            }
        }
    }
	mpz_clear(num);
	mpz_clear(residuos);
	return;
}


void adjust_n(mpz_t n,int *k){
	//divide n per k
	if(n==NULL || k==NULL || *k<=0 || (*k!=1 && *k!=3 && *k!=5 && *k!=7)){
		handle_error_with_exit("error in adjust_n\n");
	}
	if(mpz_divisible_ui_p(n,*k)==0){
		handle_error_with_exit("error in mpz_divisible adjust_n\n");
	}
	mpz_divexact_ui(n,n,*k);//n=n/k
	*k=1;
	return;
}

void reduce_array_mpz_mod_n(mpz_t*array,int length,const mpz_t a){
	if(array==NULL){	
		return;
	}
	if(length<=0 || a==NULL){
		handle_error_with_exit("error in reduce array mpz mod n\n");
	}
	
	for(int i=0;i<length;i++){
		mpz_mod(array[i],array[i],a);
	}
	return;
}
void set_to_odd_array_mpz_mod_n(mpz_t *array,int length,const mpz_t a){
	if(array==NULL){	
		return;
	}
	if(length<=0 || a==NULL){
		handle_error_with_exit("error in reduce array mpz mod n\n");
	}
	
	for(int i=0;i<length;i++){
		if(mpz_divisible_ui_p(array[i],2)!=0){//se bi è pari
			mpz_sub(array[i],a,array[i]);//bi=a-bi
		}
	}
	return;
}

 




char divide_all_by_min1(mpz_t*array_of_number,long M,mpz_t**matrix_factorization){//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione
	if(array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || M<=0){
		handle_error_with_exit("invalid parameter divide all by min1\n");
	}
	for(int i=0;i<2*M+1;i++){
		if(mpz_sgn(array_of_number[i])<0){//finquando gli elementi sono negativi
			mpz_neg(array_of_number[i],array_of_number[i]);//inverti il segno array_of_number[i]=-array_of_number[i];
			mpz_add_ui(matrix_factorization[i][1],matrix_factorization[i][1],1);//imposta la molteplicità di meno uno ad 1
		}
	}
	return 0;
}
void dg(mpz_t dig,const mpz_t y){//dg==g' derivata di g calcolata nel punto y
	if(dig==NULL || y==NULL){
		handle_error_with_exit("invalid parameter dg\n");
	}
	mpz_mul_2exp(dig,y,1);//dig=2*y
	mpz_add_ui(dig,dig,1);//dig=2y+1
	return ;
}

void g(mpz_t gy,const mpz_t y1,const mpz_t a){//funzione per calcolare le radici quadrate modulo 2^k
	mpz_t a_temp,y_temp;
	if(gy==NULL || y1==NULL || a==NULL){
		handle_error_with_exit("invalid parameter g\n");
	}
	mpz_init(y_temp);
	mpz_init_set(a_temp,a);//a_temp=a

	mpz_mul_2exp(a_temp,a_temp,1);//a_temp=a*2
	mpz_set(y_temp,y1);//y_temp=y1
	mpz_mul(y_temp,y_temp,y_temp);//y_temp=y1^2
	mpz_add(y_temp,y_temp,y1);//y_temp=y^2+y
	mpz_sub(y_temp,y_temp,a_temp);//y=y*y+y-2*a;
	mpz_set(gy,y_temp);
	
	mpz_clear(y_temp);
	mpz_clear(a_temp);
	return;
}
void calculate_new_y(mpz_t y1,const mpz_t n,long k){//calcola yt dove t=k-2 sapendo yt si può ricavare yk e calcolare la radice quadrata modulo 2^k 
	//y y corrisponde a yl dove l=k-3 da yl si puo risalire a yt e quindi a yk
	if(y1==NULL || n==NULL || k<=0 || y1==NULL || mpz_sgn(n)<=0){
		handle_error_with_exit("error in calculate new y\n");
	}
	long t=k-2;//t=k-2

	mpz_t coff_g,two_to_t,gy1,dig;

	mpz_init(coff_g);//a
	mpz_init(two_to_t);
	mpz_init(gy1);//g(y1)
	mpz_init(dig);//g'(y1)

	mpz_set(coff_g,n);//coff_g=n
	mpz_sub_ui(coff_g,coff_g,1);//coff_g=n-1
	if(mpz_divisible_2exp_p(coff_g,3)==0){//verifica divisibiità per 8
		handle_error_with_exit("error in mpz_divisible calculate new y\n");
	}
	mpz_divexact_ui(coff_g,coff_g,8);//coff_g=(n-1)/8=a

	mpz_ui_pow_ui(two_to_t,2,t);//two_to_t=2^t
	dg(dig,y1);//dig=g'(y1)=2y+1
	mpz_invert(dig,dig,two_to_t);//dig=(g'(y1))^-1 mod 2^t
	g(gy1,y1,coff_g);//y=g(y1);
	mpz_mul(dig,gy1,dig);//dig=g(y1)*(g'(y1))^-1
	mpz_sub(y1,y1,dig);//y1-g(y1)*(g'(y1))^-1
	mpz_mod(y1,y1,two_to_t);//y1==y1-g(y1)*(g'(y1))^-1 mod 2^t

	mpz_clear(coff_g);
	mpz_clear(two_to_t);
	mpz_clear(gy1);
	mpz_clear(dig);
	return;
}

char divide_all_by_2_log(long M,struct thread_data thread_data){//divide gli elementi dell'array di numeri per 2
    if(M<=0){
        handle_error_with_exit("invalid parameter divide_all_by_2_f\n");
    }
    long i=-1;
    long count=0;
    count=M & 1;//se M dispari count =1 altrimenti count=0
    count+=1;//a^2*j^2+2bj+c mod 2 =M+1+b,primo elemento ha j=-M=M mod 2
    //trovare il primo numero divisibile per 2 e poi a salti di 2 dividere gli altri
    if(mpz_divisible_2exp_p(thread_data.b,1)!=0){//b è divisibile per 2
    }
    else{
        count+=1;
    }
    //i=1 indice dispari,il primo elemento è dispari,il secondo è pari,i=0 indice pari,il primo elemento è pari
    if(count%2==0){
        i=0;//primo elemento divisibile per 2
    }
    else{
        i=1;//secondo elemento divisibile per 2
    }
    for(;i<2*M+1;i=i+2){
        thread_data.numbers[i].last_index_f_base=1;
        thread_data.numbers[i].first_index_f_base=1;
        thread_data.numbers[i].sum_log=r.log_prime[1];
    }
    thread_data.j1_mod_p[1]=-1;
    thread_data.j2_mod_p[1]=-1;
    return 1;
}

char divide_all_by_2_to_k_f(long M,struct thread_data thread_data,const mpz_t n,long k,const mpz_t a,const mpz_t b,const mpz_t y){
	//k esponente di 2 a,b coefficienti di a(j)=aj^2+2bj+c y coincide con y(k-2),se si sa y(k-2) allora a=2*y(k-2)+1

	//se n è congruo a 1 mod 8 allora n è congruo a 1 mod 4 e congruo a 1 mod 2
	char count=0;//0 se non ci sono state divisioni 1 se ci sono state divisioni
	
	if(n==NULL || a==NULL || b==NULL || y==NULL || M<=0 || mpz_sgn(n)<=0 || k<=0){
		handle_error_with_exit("invalid parameter_divide_all_by_2_to_k\n");
	}
	if(k==1){//radici quadrate di n modulo 2
		count=divide_all_by_2_log(M,thread_data);//k=1 -> dividi per 2
		return count;
	}
	return 0;//ritorna 0 se non ci sono state divisioni nell'array
}
void square_root_mod_p_to_k(mpz_t rootpk,const mpz_t x1,long p1,const mpz_t n,int k){//root r1,root r2=p^k-r1,calcola radici quadrate di n modulo p alla k
	//parametri:p1 è il primo,k è la potenza a cui elevare p, x1 è la radice quadrata di n mod p
	// calcola rootpk^2===n mod p^k,sapendo x1: x1^2==n mod p
	if(mpz_sgn(x1)<=0 || mpz_sgn(n)<=0){
		handle_error_with_exit("root zero or negative\n");
	}
	if(rootpk==NULL || n==NULL || k<=0 || p1<=1){
		handle_error_with_exit("error in parameter\n");
	}
	if(p1==2){
		handle_error_with_exit("invalid prime\n");
	}
	mpz_t q,r,r2,p,v,w,exp,x;

	mpz_init(q);
	mpz_init(r);
	mpz_init(r2);
	mpz_init(p);
	mpz_init(v);
	mpz_init(w);
	mpz_init(exp);
	mpz_init(x);
	
	mpz_set_si(p,p1);//p=p
	mpz_set(x,x1);//x=x
	mpz_pow_ui(q,p,k);//q=p^k
	if(mpz_divisible_p(q,p)==0){
				handle_error_with_exit("error in mpz_divisible square_root_mod_p_to_k\n");
	}
	mpz_divexact(r,q,p);//r=q/p
	
	mpz_mul_ui(r2,r,2);//r2=2r
	mpz_neg(r2,r2);//r2=-2r
	mpz_set(exp,q);//exp=q=p^k
	mpz_add(exp,exp,r2);//exp=q-2r
	mpz_add_ui(exp,exp,1);//exp=q-2r+1
	if(mpz_divisible_ui_p(exp,2)==0){
				handle_error_with_exit("error in mpz_divisible square_root_mod_p_to_k\n");
	}
	mpz_divexact_ui(exp,exp,2);//exp=(q-2r+1)/2
	mpz_powm(v,x,r,q);//v=x^r mod p^k=x^p^(k-1) mod p^k
	mpz_powm(w,n,exp,q);//w=n^exp mod p^k
	mpz_mul(rootpk,v,w);//rootpk=w*v
	mpz_mod(rootpk,rootpk,q);//rootpk=r1 mod p^k

	mpz_clear(q);
	mpz_clear(r);
	mpz_clear(r2);
	mpz_clear(p);
	mpz_clear(v);
	mpz_clear(w);
	mpz_clear(exp);
	mpz_clear(x);
	
	return; 
}

void calculate_root_poly_second_degree_mod_p_to_k(mpz_t j1t,const mpz_t p_to_k,long p,long k,const mpz_t a,const mpz_t b,const mpz_t n){
//calcola la radice del polinomio ax^2+2bx+c mod p^k con gcd(a,p^k)=p,gcd(a/p,p^k)=1 
	//x=(-ko)*(a/p)^-1 mod p^k ko=(b-r)/p r^2==n mod p^(k+1),x1^2==n mod p
	if(j1t==NULL || p_to_k==NULL || p<=0 || k<=0 || a==NULL || b==NULL || n==NULL || p%2==0 || (mpz_divisible_ui_p(a,2)!=0) ){
		handle_error_with_exit("calculate_root_poly_second_degree_mod_p_to_k\n");
	}
	int t;
	mpz_t r,x1,p_temp,k0,ptok1,a_temp,inv;
	mpz_init(r);
	mpz_init(x1);
	mpz_init(p_temp);
	mpz_init(k0);
	mpz_init(ptok1);
	mpz_init(a_temp);
	mpz_init(inv);

	mpz_set_si(p_temp,p);//p_temp==p
	t=quadratic_residue(x1,n,p_temp);//x1^2==n mod p
	if(t==-1 || t==0){
		handle_error_with_exit("error in calculate_root_poly_second_degrre_mod_p_to_k\n");
	}
	square_root_mod_p_to_k(r,x1,p,n,k+1);//r^2==n mod p^k+1
	mpz_sub(k0,r,b);//-ko=r-b;
	if(mpz_divisible_ui_p(k0,p)==0){//r-b non è divisibile per p
		mpz_mul_ui(ptok1,p_to_k,p);//ptok=p^(k+1)
		mpz_sub(r,ptok1,r);//r=p^(k+1)-r
		mpz_sub(k0,r,b);//-ko=r-b;
		mpz_divexact_ui(k0,k0,p);//-k0=r-b/p
	}
	else{//r-b divisibile per p
		mpz_divexact_ui(k0,k0,p);//-k0=r-b/p
	}
	mpz_set(a_temp,a);//a_temp=a
	if(mpz_divisible_ui_p(a_temp,p)==0){
				handle_error_with_exit("error in mpz_divisible calculate_root_poly_second_degree_mod_p_to_k\n");
	}
	mpz_divexact_ui(a_temp,a_temp,p);//a_temp=a/p
	mpz_invert(inv,a_temp,p_to_k);//inv=(a/p)^-1 mod p^k
	mpz_mul(j1t,inv,k0);//j1t=inv*ko=(a/p)^-1 mod p^k *(r-b/p)

	mpz_clear(p_temp);
	mpz_clear(r);
	mpz_clear(x1);
	mpz_clear(k0);
	mpz_clear(ptok1);
	mpz_clear(a_temp);
	mpz_clear(inv);
	return;
}

char divide_all_by_p_to_k_with_thread(int rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b,struct a_struct*array_a_struct){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p

	//j1=r-x0;
	//j2=-r-x0;
	if((p<=1 && p!=-1 ) || n==NULL || a==NULL || b==NULL || k<=0 || M<=0 || array_a_struct==NULL){
		handle_error_with_exit("invalid parameter divide all by p to k with thread\n");
	}
	char array_divided=0;//0 se nessuna divisione effettuata,1 se sono state effettutate divisioni per p^k
	long indexv,j1,j2,j_temp2;
	mpz_t p_to_k,r2,a_temp,inverse_a,inverse2b,j1t,j2t,c,b_temp,v,l_temp,j_temp,index,h_temp;
	mpz_init(p_to_k);
	mpz_init(r2);
	mpz_init(a_temp);
	mpz_init(inverse_a);
	mpz_init(inverse2b);
	mpz_init(j1t);
	mpz_init(j2t);
	mpz_init(c);
	mpz_init(b_temp);
	mpz_init(v);
	mpz_init(l_temp);
	mpz_init(j_temp);
	mpz_init(h_temp);
	mpz_init(index);
	mpz_set_si(p_to_k,p);//p^k=p
	mpz_set_si(r2,r.root2_n_mod_p[index_of_prime]);//r2=p^k-r seconda radice quadrata di n modulo p^k
	long j=0;//indici dell'array divisibile per p^k,se j!=0 j2 non esiste
	if(mpz_cmp_si(a,1)==0){//se a=1,infatti a^-1 mod p =1
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add_ui(j1t,j1t,rad);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
	}
	else if(!value_is_in_sorted_array(index_of_prime,array_a_struct,s)){//p non divide a
		if(r.inverse_a_mod_p[index_of_prime]==-1){
			print_array_a_struct(array_a_struct,s);
			printf("index_of_prime=%d\n",index_of_prime);
			handle_error_with_exit("error in factorize_matrix,inverse not found\n");
		}
		mpz_set_si(inverse_a,r.inverse_a_mod_p[index_of_prime]);////inverse_a=a^-1 mod p^k
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add_ui(j1t,j1t,rad);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
		mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
	}
	else{//p divide a,esiste solo una soluzione,si usa solo j1,succede solamente poche volte(circa 10)
		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){//v non divide a
			handle_error_with_exit("error in mpz_divisible divide all by p_to_k\n");
		}
		mpz_divexact(c,v,a);//c=(b^2-n)/a

		if(k==1){//p_to_k=p,esiste solo una soluzione mod p
			mpz_set(b_temp,b);
			mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
			mpz_invert(inverse2b,b_temp,p_to_k);//(2b)^-1 mod p^k
			mpz_neg(j1t,c);//j1=-c
			mpz_mul(j1t,j1t,inverse2b);//j1=-c*inverse2b;
			j=1;
		}
		else{//k>1,p_to_k>p,gcd(a,p^k)=p,gcd(a/p,p^k)=1
			handle_error_with_exit("error in divide all by p_to_k k>1\n");
			calculate_root_poly_second_degree_mod_p_to_k(j1t,p_to_k,p,k,a,b,n);
			j=1;//solo 1 soluzione
		}
	}
	mpz_mod(j1t,j1t,p_to_k);//j1t ridotto modulo p^k
	mpz_mod(j2t,j2t,p_to_k);//j1_t ridotto modulo p^k
	j1=mpz_get_si(j1t);//j1=j1t
	j2=mpz_get_si(j2t);//j2=j2t
	j_temp2=j1;
	indexv=j_temp2+M;//l'indice deve essere positivo
	while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
			//index minore di first-> diventa il nuovo first
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
		}
		else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
			//index maggiore di last-> diventa il nuovo last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2+=p;//ad ogni ciclo aggiungo p
		indexv+=p;
	}
	j_temp2=-p+j1;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
	indexv=j_temp2+M;
	while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
			//index minore di first-> diventa il nuovo first
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
		}
		else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
			//index maggiore di last-> diventa il nuovo last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2-=p;//ad ogni ciclo tolgo p
		indexv-=p;
	}
	if(j==1){
		mpz_clear(p_to_k);
		mpz_clear(r2);
		mpz_clear(a_temp);
		mpz_clear(inverse_a);
		mpz_clear(inverse2b);
		mpz_clear(j1t);
		mpz_clear(j2t);
		mpz_clear(c);
		mpz_clear(b_temp);
		mpz_clear(v);
		mpz_clear(l_temp);
		mpz_clear(j_temp);
		mpz_clear(h_temp);
		mpz_clear(index);
		if(array_divided==1){
			return 1;//ritorna 1 se c'è stata almeno 1 divisione
		}
		return 0;//ritorna 0 se non ci sono state divisioni nell'array
	}
	j_temp2=j2;//si passa alla secodna radice
	indexv=j_temp2+M;
	while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
			//index minore di first-> diventa il nuovo first
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
		}
		else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
			//index maggiore di last-> diventa il nuovo last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2+=p;//ad ogni ciclo aggiungo p
		indexv+=p;
	}
	j_temp2=-p+j2;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
	indexv=j_temp2+M;
	while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else if(index_of_prime<thread_data.numbers[indexv].first_index_f_base){
			//index minore di first-> diventa il nuovo first
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
		}
		else if(index_of_prime>thread_data.numbers[indexv].last_index_f_base){
			//index maggiore di last-> diventa il nuovo last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2-=p;//ad ogni ciclo tolgo p
		indexv-=p;
	}
	mpz_clear(p_to_k);
	mpz_clear(r2);
	mpz_clear(a_temp);
	mpz_clear(inverse_a);
	mpz_clear(inverse2b);
	mpz_clear(j1t);
	mpz_clear(j2t);
	mpz_clear(c);
	mpz_clear(b_temp);
	mpz_clear(v);
	mpz_clear(l_temp);
	mpz_clear(j_temp);
	mpz_clear(h_temp);
	mpz_clear(index);
	if(array_divided==1){
		return 1;//ritorna 1 se c'è stata almeno 1 divisione
	}
	return 0;//ritorna 0 se non ci sono state divisioni nell'array
}

char divide_all_by_p_to_k_f(int rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b,struct a_struct*array_a_struct,int*index_array_a_struct,int s){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p

    //j1=r-x0;
    //j2=-r-x0;
    if((p<=1 && p!=-1 ) || n==NULL || a==NULL || b==NULL || k<=0 || M<=0 || (array_a_struct==NULL && s>0) || index_array_a_struct==NULL || s<0){
        handle_error_with_exit("invalid parameter divide all by p to k\n");
    }
    char array_divided=0;//0 se nessuna divisione effettuata,1 se sono state effettutate divisioni per p^k
    long indexv,j1,j2,j_temp2;
    mpz_t p_to_k,r2,a_temp,inverse_a,inverse2b,j1t,j2t,c,b_temp,v,l_temp,j_temp,index,h_temp;
    mpz_init(p_to_k);
    mpz_init(r2);
    mpz_init(a_temp);
    mpz_init(inverse_a);
    mpz_init(inverse2b);
    mpz_init(j1t);
    mpz_init(j2t);
    mpz_init(c);
    mpz_init(b_temp);
    mpz_init(v);
    mpz_init(l_temp);
    mpz_init(j_temp);
    mpz_init(h_temp);
    mpz_init(index);
    mpz_set_si(p_to_k,p);//p^k=p
    mpz_set_si(r2,r.root2_n_mod_p[index_of_prime]);//r2=p^k-r seconda radice quadrata di n modulo p^k
    long j=0;//indici dell'array divisibile per p^k,se j!=0 j2 non esiste
    if(mpz_cmp_si(a,1)==0){//se a=1,infatti a^-1 mod p =1
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add_ui(j1t,j1t,rad);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
    }

    else if((*index_array_a_struct<s && array_a_struct[*index_array_a_struct].index_prime_a!=index_of_prime) || *index_array_a_struct>=s){//p non divide a
        if(r.inverse_a_mod_p[index_of_prime]==-1){
            printf("p=%ld,index_of_prime=%d,index_array_a_struct=%d\n",p,index_of_prime,array_a_struct[*index_array_a_struct].index_prime_a);
            handle_error_with_exit("error in factorize_matrix,inverse not found\n");
        }
        mpz_set_si(inverse_a,r.inverse_a_mod_p[index_of_prime]);////inverse_a=a^-1 mod p^k
        mpz_neg(j1t,b);//j1=-b
        mpz_neg(j2t,b);//j2=-b
        mpz_add_ui(j1t,j1t,rad);//j1=-b+r
        mpz_add(j2t,j2t,r2);//j2=-b+r2
        mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
        mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
    }
    else if(array_a_struct[*index_array_a_struct].index_prime_a==index_of_prime){//p divide a,esiste solo una soluzione,si usa solo j1,succede solamente poche volte(circa 10)
        //calcolo di c
        if(mpz_divisible_ui_p(a,p)==0){
            handle_error_with_exit("errore p non divide a\n");
        }
		(*index_array_a_struct)++;
        mpz_mul(v,b,b);//v=b^2
        mpz_sub(v,v,n);//v=b^2-n
        if(mpz_divisible_p(v,a)==0){//v non è divisibile per a
            handle_error_with_exit("error in mpz_divisible divide all by p_to_k\n");
        }
        mpz_divexact(c,v,a);//c=(b^2-n)/a

        if(k==1){//p_to_k=p,esiste solo una soluzione mod p
            mpz_set(b_temp,b);
            mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
            mpz_invert(inverse2b,b_temp,p_to_k);//(2b)^-1 mod p^k
            mpz_neg(j1t,c);//j1=-c
            mpz_mul(j1t,j1t,inverse2b);//j1=-c*inverse2b;
            j=1;
        }
        else{//k>1,p_to_k>p,gcd(a,p^k)=p,gcd(a/p,p^k)=1
        	handle_error_with_exit("error in divide all by p_to_k k>1\n");
            calculate_root_poly_second_degree_mod_p_to_k(j1t,p_to_k,p,k,a,b,n);
            j=1;//solo 1 soluzione
        }
    }
    else{
        handle_error_with_exit("caso non gestito in divide all by p_to_k\n");
    }
    mpz_mod(j1t,j1t,p_to_k);//j1t ridotto modulo p^k
    mpz_mod(j2t,j2t,p_to_k);//j1_t ridotto modulo p^k
	j1=mpz_get_si(j1t);//j1=j1t
	j2=mpz_get_si(j2t);//j2=j2t
    thread_data.j1_mod_p[index_of_prime]=j1;
    if(j!=1) {
        thread_data.j2_mod_p[index_of_prime] = j2;
    }
    else{
        thread_data.j2_mod_p[index_of_prime]= -1;
    }
    j_temp2=j1;
    indexv=j_temp2+M;//l'indice deve essere positivo
	while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else{//altrimenti modifica solamente last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2+=p;//ad ogni ciclo aggiungo p
		indexv+=p;
	}
	j_temp2=-p+j1;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp2+M;
	while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else{//altrimenti modifica solamente last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2-=p;//ad ogni ciclo tolgo p
		indexv-=p;
	}
    if(j==1){
        mpz_clear(p_to_k);
        mpz_clear(r2);
        mpz_clear(a_temp);
        mpz_clear(inverse_a);
        mpz_clear(inverse2b);
        mpz_clear(j1t);
        mpz_clear(j2t);
        mpz_clear(c);
        mpz_clear(b_temp);
        mpz_clear(v);
        mpz_clear(l_temp);
        mpz_clear(j_temp);
        mpz_clear(h_temp);
        mpz_clear(index);
        if(array_divided==1){
            return 1;//ritorna 1 se c'è stata almeno 1 divisione
        }
        return 0;//ritorna 0 se non ci sono state divisioni nell'array
    }
	j_temp2=j2;//si passa alla secodna radice
    indexv=j_temp2+M;
	while(j_temp2<=M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else{//altrimenti modifica solamente last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2+=p;//ad ogni ciclo aggiungo p
		indexv+=p;
	}
	j_temp2=-p+j2;//all'inizio j_temp=-p+j1t,poi diventa -k*p+j1t(a salti di p negativi)
    indexv=j_temp2+M;
	while(j_temp2>=-M){//all'inizio j_temp=0*p+j1t,poi diventa k*p+j1t(a salti di p)
		//indexv=j_temp2+M;//siccome j_temp è sfalsato di -M si riaggiunge M
		thread_data.numbers[indexv].sum_log+=r.log_prime[index_of_prime];
		//se il first_index_f_base non è stato modificato setta first e last
		if(thread_data.numbers[indexv].first_index_f_base==-1){
			thread_data.numbers[indexv].first_index_f_base=index_of_prime;
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		else{//altrimenti modifica solamente last
			thread_data.numbers[indexv].last_index_f_base=index_of_prime;
		}
		if(thread_data.numbers[indexv].first_index_f_base>thread_data.numbers[indexv].last_index_f_base){
			handle_error_with_exit("error in index\n");
		}
		array_divided=1;//una divisione è stata effettuata
		j_temp2-=p;//ad ogni ciclo tolgo p
		indexv-=p;
	}
    mpz_clear(p_to_k);
    mpz_clear(r2);
    mpz_clear(a_temp);
    mpz_clear(inverse_a);
    mpz_clear(inverse2b);
    mpz_clear(j1t);
    mpz_clear(j2t);
    mpz_clear(c);
    mpz_clear(b_temp);
    mpz_clear(v);
    mpz_clear(l_temp);
    mpz_clear(j_temp);
    mpz_clear(h_temp);
    mpz_clear(index);
    if(array_divided==1){
        return 1;//ritorna 1 se c'è stata almeno 1 divisione
    }
    return 0;//ritorna 0 se non ci sono state divisioni nell'array
}

void*thread_factorization_job(void*arg){
    long p;
    struct factorization_thread_data*factorization_thread_data=arg;
    int start=(*factorization_thread_data).start;
    int end=(*factorization_thread_data).end;
    mpz_t a1;
    mpz_init(a1);
    if((*factorization_thread_data).id_thread==1){
    }
    if((*factorization_thread_data).is_a_default==1){
    	mpz_set_si(a1,1);
    }
    else{
    	mpz_set(a1,a_old);
    }
    for(int i=start;i<=end;i++){
    	if(i==0){//skippa indice 0 e indice 1
    		continue;
    	}
    	if(i==1){
			divide_all_by_2_log(M, (*factorization_thread_data).thread_data);
			continue;
		}
        p=r.prime[i];//primo iesimo della factor base
        divide_all_by_p_to_k_with_thread(r.root_n_mod_p[i],p,i,1,M,(*factorization_thread_data).thread_data,n,a1,(*factorization_thread_data).thread_data.b,array_a_struct);
    }
    mpz_clear(a1);
    return NULL;
}

void factor_matrix_f(const mpz_t n,long M,struct thread_data thread_data,int cardinality_factor_base,const mpz_t a,
                     struct a_struct*array_a_struct,int s){
    if(n==NULL || a==NULL || mpz_sgn(n)<=0 || M<=0 || cardinality_factor_base<=0 || (array_a_struct==NULL && s>0) || s<0){
        handle_error_with_exit("error in factor matrix_f\n");
    }
    long p;
    int index_array_a_struct=0;
    for (int i=0;i<cardinality_factor_base;i++){//per ogni elemento della factor base
        p=r.prime[i];//primo iesimo della factor base
        if(p==-1){
            thread_data.j1_mod_p[i]=-1;
            thread_data.j2_mod_p[i]=-1;
            continue;
        }
        if(p==2){
            divide_all_by_2_log(M,thread_data);
        }
        else{//p>2 e dispari
            divide_all_by_p_to_k_f(r.root_n_mod_p[i],p,i,1,M,thread_data,n,a,thread_data.b,array_a_struct,&index_array_a_struct,s);
        }
    }
    return;
}

void add_exponent_of_a(long**matrix_B_smooth,int num_B_smooth,int s,const mpz_t a,long* array_of_prime_chosen_for_a){
	long index;
	if(a==NULL || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0 || s<0 || mpz_sgn(a)<=0){
		handle_error_with_exit("error in add exponent_of_a\n");
	}
	if(mpz_cmp_si(a,1)==0){//se a=1
		if(array_of_prime_chosen_for_a==NULL){//e se array==NULL
			return;
		}
		else{
			handle_error_with_exit("error in parameter add exponent_of_a\n");
		}
	}
	for(int i=0;i<s;i++){//per tutta la lunghezza dell'array dei primi scelti
		index=array_of_prime_chosen_for_a[2*i];//ritorna l'indice del numero primo nella factor base
		for(int j=0;j<num_B_smooth;j++){//per tutte le righe della matrice dei B_smooth
			matrix_B_smooth[j][2*index+1]+=1;//aumenta l'esponente di 1 sui primi scelti
		}
	}
	return;
}




