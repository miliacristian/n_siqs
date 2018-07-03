#include "math_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "gmp.h"
#include "matrix_function.h"
#include "print.h"
#include "list_factorization.h"
#include "list_square_relation.h"
#include <unistd.h>
#include <mpfr.h>

extern struct row_factorization r;
extern long M;
extern mpz_t n;
extern mpz_t a_old;
extern struct a_struct*array_a_struct;
extern int cardinality_factor_base;
extern int s;
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
/*char factorize_num_B_smooth(int*array_factorization,int len_array,struct row *row){
	mpz_t temp;
	mpz_init(temp);
	if(array_factorization==NULL || len_array<=0 || row==NULL){
		handle_error_with_exit("error in factorize num_B_smooth\n");
	}
	for(int i=0;i<len_array/2;i++){
		mpz_set(temp,row->num);
		array_factorization[2*i]=r.prime[i+row->index_first_prime];
		while(mpz_divisible_ui_p(temp,r.prime[i+row->index_first_prime])!=0){
			mpz_divexact_ui(temp,temp,r.prime[i+row->index_first_prime]);
			array_factorization[2*i+1]+=1;
		}
	}
	if(mpz_cmp_si(temp,1)==0 ||
		return 1;
	}
	if(mpz_cmp_si(temp,-1)==0){
		return 2;
	}
	mpz_clear(temp);
	return 0;
}*/
void calculate_square(mpz_t square,const mpz_t a,int index,const mpz_t b,const mpz_t n){
    // -M<index<M
    mpz_mul_si(square,a,index);//a*j
    mpz_add(square,square,b);//a*j+b
	mpz_mod(square,square,n);// a*j+b mod n
    return;
}
char verify_factorization(const mpz_t num,mpz_t residuos,struct node_factorization*head_factor,const mpz_t a){
	//num*a deve essere uguale a fattorizzazione*residuo
	gmp_printf("inizio verifica num=%Zd\n",num);
	mpz_t temp,num_temp;
	long exp_of_factor,factor;
	mpz_t factor_raise_to_exp;
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
	mpz_init(residuos);
	if(num_B_smooth==NULL || num_semi_B_smooth==NULL || num_potential_B_smooth==NULL || M<=0 ||
	   head_square==NULL || tail_square==NULL || head_residuos==NULL
	   || tail_residuos==NULL || (array_a_struct==NULL && s>0)){
		handle_error_with_exit("error in find_list_square_relation");
	}
	struct square_relation square_relation;
	mpz_init(num);
	for(long i=0;i<2*M+1;i++){
        if(thread_data.numbers[i].sum_log>=thread_data.log_thresold){
            //possibile B_smooth trovato
            (*num_potential_B_smooth)++;
            create_num(num,a,thread_data.b,n,thread_data.numbers[i].j);
            //head_factor=factorize_num_v1(num,thread_data.numbers[i].first_index_f_base,thread_data.numbers[i].last_index_f_base,&is_B_smooth,&is_semi_B_smooth,residuos,array_a_struct,s);
			head_factor=factorize_num_v2(num,thread_data.numbers[i].j,thread_data.numbers[i].first_index_f_base,thread_data.numbers[i].last_index_f_base,&is_B_smooth,&is_semi_B_smooth,residuos,array_a_struct,s,thread_data);
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
            }
            else if(is_semi_B_smooth==1){
				(*num_semi_B_smooth)++;
				is_B_smooth=0;
				is_semi_B_smooth=0;
				square_relation.head_factorization=head_factor;
				mpz_init(square_relation.square);
				mpz_init(square_relation.num);
				mpz_init(square_relation.residuos);
				mpz_set(square_relation.num,num);
				mpz_set(square_relation.residuos,residuos);
				calculate_square(square_relation.square,a,i-M,thread_data.b,n);
                insert_at_tail_square_relation(square_relation,head_residuos,tail_residuos);
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

int max(int i,int j){//calcola il massimo tra i e j
	if(i<=j){
		return j;
	}
	return i;
}

int min(int i,int j){//calcola il minimo tra i e j
	if(i<=j){
		return i;
	}
	return j;
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

 long inverse_mod_n(long x,long n){
	if(x<0){
		handle_error_with_exit("negative x\n");
	}
	if(x==0){
		handle_error_with_exit("impossible inverse zero\n");
	}
	if(n<=0){
		handle_error_with_exit("error in parameter inverse_mod_n\n");
	}
	reduce_mod_n(&x,n);
	if(x==1){
		return 1;
	}
	long inverse;
	long result[3];
	xgcd(result,x,n);
	if(result[0]!=1){
		handle_error_with_exit("impossible inverse\n");
	}
	inverse=result[1];
	reduce_mod_n(&inverse,n);
	return inverse;
}

int quadratic_residue(mpz_t x,const mpz_t q,const mpz_t n)
{
// find x^2 = q mod n
// return
// -1 q is quadratic non-residue mod n
//  1 q is quadratic residue mod n
//  0 q is congruent to 0 mod n
//
	if(x==NULL || q==NULL || n==NULL){
		handle_error_with_exit("error in quadratic_residue\n");
	}
    int                        leg;
    mpz_t                        tmp,ofac,nr,t,r,c,b;
    unsigned int            mod4;
    mp_bitcnt_t                twofac=0,m,i,ix;

    mod4=mpz_tstbit(n,0);
    if(!mod4) // must be odd
        return 0;

    mod4+=2*mpz_tstbit(n,1);

    leg=mpz_legendre(q,n);
    if(leg!=1)
        return leg;

    mpz_init_set(tmp,n);

    if(mod4==3) // directly, x = q^(n+1)/4 mod n
        {
        mpz_add_ui(tmp,tmp,1UL);
        mpz_tdiv_q_2exp(tmp,tmp,2);
        mpz_powm(x,q,tmp,n);
        mpz_clear(tmp);
        }
    else // Tonelli-Shanks
        {
        mpz_inits(ofac,t,r,c,b,NULL);

        // split n - 1 into odd number times power of 2 ofac*2^twofac
        mpz_sub_ui(tmp,tmp,1UL);
        twofac=mpz_scan1(tmp,twofac); // largest power of 2 divisor
        if(twofac)
            mpz_tdiv_q_2exp(ofac,tmp,twofac); // shift right

        // look for non-residue
        mpz_init_set_ui(nr,2UL);
        while(mpz_legendre(nr,n)!=-1)
            mpz_add_ui(nr,nr,1UL);

        mpz_powm(c,nr,ofac,n); // c = nr^ofac mod n

        mpz_add_ui(tmp,ofac,1UL);
        mpz_tdiv_q_2exp(tmp,tmp,1);
        mpz_powm(r,q,tmp,n); // r = q^(ofac+1)/2 mod n

        mpz_powm(t,q,ofac,n);
        mpz_mod(t,t,n); // t = q^ofac mod n

        if(mpz_cmp_ui(t,1UL)!=0) // if t = 1 mod n we're done
            {
            m=twofac;
            do
                {
                i=2;
                ix=1;
                while(ix<m)
                    {
                    // find lowest 0 < ix < m | t^2^ix = 1 mod n
                    mpz_powm_ui(tmp,t,i,n); // repeatedly square t
                    if(mpz_cmp_ui(tmp,1UL)==0)
                        break;
                    i<<=1; // i = 2, 4, 8, ...
                    ix++; // ix is log2 i
                    }
                mpz_powm_ui(b,c,1<<(m-ix-1),n); // b = c^2^(m-ix-1) mod n
                mpz_mul(r,r,b);
                mpz_mod(r,r,n); // r = r*b mod n
                mpz_mul(c,b,b);
                mpz_mod(c,c,n); // c = b^2 mod n
                mpz_mul(t,t,c);
                mpz_mod(t,t,n); // t = t b^2 mod n
                m=ix;
                }while(mpz_cmp_ui(t,1UL)!=0); // while t mod n != 1
            }
        mpz_set(x,r);
        mpz_clears(tmp,ofac,nr,t,r,c,b,NULL);
        }

    return 1;
}


void reduce_mod_n(long *a,long n){//reduce a mod n,rendendo a positivo maggiore o uguale a zero e minore di n,n è positivo
	if(a==NULL || n<=0){
		handle_error_with_exit("error in parameter reduce mod n\n");
	}
	if(*a<0){
		*a=(*a)%n;
		if(*a!=0){
			*a=*a+n;
		}
		return;
	}
	if(*a<n && *a>=0){
		return;
	}
	*a=*a%n;
	return;
}
void reduce_int_mod_n(int *a,int n){//reduce a mod n,rendendo a positivo maggiore o uguale a zero e minore di n,n è positivo
	if(a==NULL || n<=0){
		handle_error_with_exit("error in parameter reduce mod n\n");
	}
	if(*a<0){
		*a=(*a)%n;
		if(*a!=0){
			*a=*a+n;
		}
		return;
	}
	if(*a<n && *a>=0){
		return;
	}
	*a=*a%n;
	return;
}
int reduce_mod_2(int a){
	if(a!=1 && a!=0 && a!=-1){
		handle_error_with_exit("error in reduce_mod_2\n");
	}
	if(a==0){
		return 0;
	}
	//a=-1 o a=1 ->ritorna 1
	return 1;
}
int reduce_int_mod_n_v2(int a,int n){//reduce a mod n,rendendo a positivo maggiore o uguale a zero e minore di n,n è positivo
	if(n<=0){
		handle_error_with_exit("error in parameter reduce mod n\n");
	}
	if(a==0){
		return a;
	}
	if(a<0){
		a=(a)%n;
		if(a<0){
			a=a+n;
		}
		return a;
	}
	if(a<n && a>=0){
		return a;
	}
	a=a%n;
	if(a>=n || a<0){
		handle_error_with_exit("error in reduce int mod n\n");
		return a;
	}
	return a;
}
long calculate_mod_n(long a,long n){//reduce a mod n,rendendo a positivo maggiore o uguale a zero e minore di n,n è positivo
	if(n<=0){
		handle_error_with_exit("error in parameter calculate mod n\n");
	}
	if(a<0){
		a=a%n;
		if(a!=0){
			a=a+n;
		}
		return a;
	}
	if(a<n && a>=0){
		return a;
	}
	a=a%n;
	return a;
}

long power_mod_n(long base,long exponent,long n){// base^exponent mod n con base,esponente ed n interi positivi
	if(base<=0 || exponent<0 || n<=0){
		handle_error_with_exit("error in parameter power_mod_n\n");
	}
	if (n==1){//ogni numero modulo 1 è 0
		return 0;
	}
	if(exponent==1){//se esponent è uguale a 1 basta semplicemente ridurre la base modulo n
		reduce_mod_n(&base,n);//riduce la base mod n
		return base;
	}
	long result=1,temp_exponent;
	reduce_mod_n(&base,n);//per calcolare base^exponent si può ridurre la base modulo n
	while(exponent>0){//finquando l'esponente non è diventato 0
		temp_exponent=exponent;
		reduce_mod_n(&temp_exponent,2);//temp exponent è uguale all'ultima cifra binaria di exponent
		if(temp_exponent==1){//se le cifre binarie sono uguali a 1 va fatta la moltiplicazione tra le quadrature trovate e poi vanno ridotte modulo n
			result=result*base;
			reduce_mod_n(&result,n);
		}
		exponent=exponent>>1;//toglie l'ultima cifra binaria,alla fine exponent sarà uguale a 0
		base=base*base;//quadrature successive di base
		reduce_mod_n(&base,n);//riduzione della quadratura della base
	}
	return result;
}

void xgcd(long result [3], long x, long y){//algoritmo di euclide esteso ax+by=gcd(x,y)
	//result[0]=gcd(x,y),result[1]=a,result[2]=b,l'array result può essere passato vuoto
    if(x<0 || y<0 ){
	handle_error_with_exit("error in parameter xgcd\n");
    }
    long aa[2]={1,0}, bb[2]={0,1}, q;
    while(1) {
        q = x / y;
 	x = x % y;
        aa[0] = aa[0] - q*aa[1];
  	bb[0] = bb[0] - q*bb[1];
        if (x == 0) {
        	result[0] = y;
		result[1] = aa[1];
		result[2] = bb[1];
        	return;
        }
        q = y / x;
	y = y % x;
        aa[1] = aa[1] - q*aa[0];
	bb[1] = bb[1] - q*bb[0];
        if (y == 0) {
        	result[0] = x;
 		result[1] = aa[0];
 		result[2] = bb[0];
        	return;
        }
    }
}

long gcd(long a,long b)//a e b interi positivi
{
	printf("a=%ld b=%ld\n",a,b);
     if (a<0 || b<0){
		handle_error_with_exit("error in parameter gcd\n");
     }
    long temp;
    while (b != 0)//b è il resto della divisione tra a e b,quando il resto è 0 il penultimo resto è il gcd tra a e b
    {
        temp = a % b;

        a = b;
        b = temp;
    }
    return a;
}//
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
/*char divide_all_by_2(mpz_t*array_of_number,long M,mpz_t**matrix_factorization){//divide gli elementi dell'array di numeri per 2
	if(array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || M<=0){
		handle_error_with_exit("invalid parameter divide_all_by_2\n");
	}
	long i;
	//trovare il primo numero divisibile per 2 e poi a salti di 2 dividere gli altri
	if(mpz_divisible_2exp_p(array_of_number[0],1)!=0){//è divisibile per 2
		i=0;//indice pari
	}
	else{
		i=1;//indice dispari
	}
	for(;i<2*M+1;i=i+2){
			if(mpz_divisible_2exp_p(array_of_number[i],1)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 2\n");
			}
			mpz_divexact_ui(array_of_number[i],array_of_number[i],2);//dividi il numero nell'array per 2
			mpz_add_ui(matrix_factorization[i][3],matrix_factorization[i][3],1);//aumenta di 1 la potenza del primo 2
	}
	return 1;
}*/
/*char divide_all_by_2_f(long M,struct thread_data thread_data){//divide gli elementi dell'array di numeri per 2
	if(M<=0){
		handle_error_with_exit("invalid parameter divide_all_by_2_f\n");
	}
	long i;
	//trovare il primo numero divisibile per 2 e poi a salti di 2 dividere gli altri
	if(mpz_divisible_2exp_p(mat->row[0].num,1)!=0){//è divisibile per 2
		i=0;//indice pari
	}
	else{
		i=1;//indice dispari
	}
	for(;i<2*M+1;i=i+2){
			if(mpz_divisible_2exp_p(mat->row[i].num,1)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 2\n");
			}
			mat->row[i].sum_log+=r.log_prime[1];//somma i log
			mat->row[i].index_first_prime=1;
			mat->row[i].index_last_prime=1;
			//mpz_divexact_ui(mat.row[i].num,mat.row[i].num,2);//dividi il numero nell'array per 2
			//mat.row[i].factorization[3]+=1;
	}
	return 1;
}*/
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
/*char divide_all_by_4(mpz_t*array_of_number,long M,mpz_t**matrix_factorization,const mpz_t n,const mpz_t a,const mpz_t b){//ritorna 1 se c'è stata almeno 1 divisione,zero altrimenti
	if(n==NULL || a==NULL || b==NULL || array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || M<=0 || mpz_sgn(n)<=0){
		handle_error_with_exit("invalid parameter divide all_by_4\n");
	}
	char array_divided=0;
	int index_of_prime=2;//2 ha indice 2 
	long r=1;//radice quadrate modulo 4,se n==1 mod 4
	long r2=3;//radice quadrate modulo 4,se n==1 mod 4
	long j1,j2;
	
	long j=0;//elementi dell'array divisibili per 4 + salti di 4 
	long p_to_k=4;//stiamo considerando p=2 k=2
	
	long l,h;//interi per calcolare i multipli di 4
	long index;//indice nell'array degli elementi divisibili per 4
	long inverse_a,inverse2b;

	mpz_t c;//c coefficiente di a(j)
	mpz_t j1t,j2t;
	mpz_t a_temp,b_temp,v,ptok4;

	mpz_init(a_temp);
	mpz_init(b_temp);
	mpz_init(j1t);
	mpz_init(j2t);
	mpz_init(v);
	mpz_init(c);
	mpz_init(ptok4);

	if(mpz_divisible_2exp_p(a,1)==0){//se 2 non divide a vedi note
		mpz_set_si(ptok4,p_to_k);//ptok4=p_to_k
		mpz_invert(a_temp,a,ptok4);//a_temp=a^-1 mod 4
		inverse_a=mpz_get_si(a_temp);
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add_ui(j1t,j1t,r);//j1=-b+r
		mpz_add_ui(j2t,j2t,r2);//j2=-b+r
		mpz_mul_ui(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul_ui(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
	}
	else{//2 divide a,esiste solo una soluzione,si usa solo j1
		handle_error_with_exit("2 is a factor of a!\n");
		mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
		mpz_invert(b_temp,b_temp,ptok4);//(2b)^-1 mod p^k
		inverse2b=mpz_get_si(b_temp);
		
		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 4\n");
		}
		mpz_divexact(c,v,a);//c=(b^2-n)/a
		
		mpz_neg(j1t,c);//j1=-c
		mpz_mul_ui(j1t,j1t,inverse2b);//j1=-c*inverse2b;
		j=1;
	}
	mpz_mod(j1t,j1t,ptok4);
	mpz_mod(j2t,j2t,ptok4);
	j1=mpz_get_si(j1t);
	j2=mpz_get_si(j2t);//x0+j=r mod 2^k -> j=r-xo mod 2^k,restituisce j>=0,però ricordarsi che j va da -M a M
	l=0;
	
	while(j1+l*p_to_k<=M){//soluzioni di r1 con l positivo
		index=j1+l*p_to_k+M;
		if(mpz_divisible_2exp_p(array_of_number[index],1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[index],array_of_number[index],2);//dividi il numero nell'array per 2
		mpz_add_ui(matrix_factorization[index][index_of_prime+1],matrix_factorization[index][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	while(j1-l*p_to_k>=-M){//soluzioni di r1 con l negativo
		index=j1-l*p_to_k+M;
		if(mpz_divisible_2exp_p(array_of_number[index],1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[index],array_of_number[index],2);//dividi il numero nell'array per 2
		mpz_add_ui(matrix_factorization[index][index_of_prime+1],matrix_factorization[index][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
	}
	h=0;
	while(j2+h*p_to_k<=M && j==0){//r2+p^k,soluzioni di r2 con h positivo
		index=j2+h*p_to_k+M;
		if(mpz_divisible_2exp_p(array_of_number[index],1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[index],array_of_number[index],2);//dividi il numero nell'array per 2
		mpz_add_ui(matrix_factorization[index][index_of_prime+1],matrix_factorization[index][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
	}
	h=1;
	while(j2-h*p_to_k>=-M && j==0){//soluzioni di r2 con h negativo
		index=j2-h*p_to_k+M;
		if(mpz_divisible_2exp_p(array_of_number[index],1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[index],array_of_number[index],2);//dividi il numero nell'array per 2
		mpz_add_ui(matrix_factorization[index][index_of_prime+1],matrix_factorization[index][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
	}
	mpz_clear(c);
	mpz_clear(j1t);
	mpz_clear(j2t);
	mpz_clear(a_temp);
	mpz_clear(b_temp);
	mpz_clear(v);
	mpz_clear(ptok4);
	if(array_divided==1){
		return 1;//ritorna 1 se c'è stata almeno 1 divisione
	}
	return 0;//ritorna 0 se non ci sono state divisioni nell'array
}
char divide_all_by_4_f(long M,struct matrix_factorization *mat,const mpz_t n,const mpz_t a,const mpz_t b){//ritorna 1 se c'è stata almeno 1 divisione,zero altrimenti
	if(n==NULL || a==NULL || b==NULL || M<=0 || mpz_sgn(n)<=0 || mat==NULL){
		handle_error_with_exit("invalid parameter divide all_by_4\n");
	}
	char array_divided=0;
	int index_of_prime=2;//2 ha indice 2 
	long rad=1;//radice quadrate modulo 4,se n==1 mod 4
	long r2=3;//radice quadrate modulo 4,se n==1 mod 4
	long j1,j2;
	
	long j=0;//elementi dell'array divisibili per 4 + salti di 4 
	long p_to_k=4;//stiamo considerando p=2 k=2
	
	long l,h;//interi per calcolare i multipli di 4
	long index;//indice nell'array degli elementi divisibili per 4
	long inverse_a,inverse2b;

	mpz_t c;//c coefficiente di a(j)
	mpz_t j1t,j2t;
	mpz_t a_temp,b_temp,v,ptok4;

	mpz_init(a_temp);
	mpz_init(b_temp);
	mpz_init(j1t);
	mpz_init(j2t);
	mpz_init(v);
	mpz_init(c);
	mpz_init(ptok4);

	if(mpz_divisible_2exp_p(a,1)==0){//se 2 non divide a vedi note
		mpz_set_si(ptok4,p_to_k);//ptok4=p_to_k
		mpz_invert(a_temp,a,ptok4);//a_temp=a^-1 mod 4
		inverse_a=mpz_get_si(a_temp);
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add_ui(j1t,j1t,rad);//j1=-b+r
		mpz_add_ui(j2t,j2t,r2);//j2=-b+r
		mpz_mul_ui(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul_ui(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
	}
	else{//2 divide a,esiste solo una soluzione,si usa solo j1
		handle_error_with_exit("2 is a factor of a!\n");
		mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
		mpz_invert(b_temp,b_temp,ptok4);//(2b)^-1 mod p^k
		inverse2b=mpz_get_si(b_temp);
		
		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 4\n");
		}
		mpz_divexact(c,v,a);//c=(b^2-n)/a
		
		mpz_neg(j1t,c);//j1=-c
		mpz_mul_ui(j1t,j1t,inverse2b);//j1=-c*inverse2b;
		j=1;
	}
	mpz_mod(j1t,j1t,ptok4);
	mpz_mod(j2t,j2t,ptok4);
	j1=mpz_get_si(j1t);
	j2=mpz_get_si(j2t);//x0+j=r mod 2^k -> j=r-xo mod 2^k,restituisce j>=0,però ricordarsi che j va da -M a M
	l=0;
	
	while(j1+l*p_to_k<=M){//soluzioni di r1 con l positivo
		index=j1+l*p_to_k+M;
		if(mpz_divisible_2exp_p(mat->row[index].num,1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mat->row[index].sum_log=mat->row[index].sum_log+r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[index].num,mat.row[index].num,2);//dividi il numero nell'array per 2
		//mat.row[index].factorization[index_of_prime+1]+=1; 
		array_divided=1;//una divisione è stata effettuata
		l++;
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	while(j1-l*p_to_k>=-M){//soluzioni di r1 con l negativo
		index=j1-l*p_to_k+M;
		if(mpz_divisible_2exp_p(mat->row[index].num,1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mat->row[index].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[index].num,mat.row[index].num,2);//dividi il numero nell'array per 2
		//mat.row[index].factorization[index_of_prime+1]+=1; 
		array_divided=1;//una divisione è stata effettuata
		l++;
	}
	h=0;
	while(j2+h*p_to_k<=M && j==0){//r2+p^k,soluzioni di r2 con h positivo
		index=j2+h*p_to_k+M;
		if(mpz_divisible_2exp_p(mat->row[index].num,1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mat->row[index].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[index].num,mat.row[index].num,2);//dividi il numero nell'array per 2
		//mat.row[index].factorization[index_of_prime+1]+=1; 
		array_divided=1;//una divisione è stata effettuata
		h++;
	}
	h=1;
	while(j2-h*p_to_k>=-M && j==0){//soluzioni di r2 con h negativo
		index=j2-h*p_to_k+M;
		if(mpz_divisible_2exp_p(mat->row[index].num,1)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mat->row[index].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[index].num,mat.row[index].num,2);//dividi il numero nell'array per 2
		//mat.row[index].factorization[index_of_prime+1]+=1; 
		array_divided=1;//una divisione è stata effettuata
		h++;
	}
	mpz_clear(c);
	mpz_clear(j1t);
	mpz_clear(j2t);
	mpz_clear(a_temp);
	mpz_clear(b_temp);
	mpz_clear(v);
	mpz_clear(ptok4);
	if(array_divided==1){
		return 1;//ritorna 1 se c'è stata almeno 1 divisione
	}
	return 0;//ritorna 0 se non ci sono state divisioni nell'array
}*/
/*char divide_all_by_2_to_k(mpz_t*array_of_number,long M,mpz_t**matrix_factorization,const mpz_t n,long k,const mpz_t a,const mpz_t b,const mpz_t y){
	//k esponente di 2 a,b coefficienti di a(j)=aj^2+2bj+c y coincide con y(k-2),se si sa y(k-2) allora a=2*y(k-2)+1

	//se n è congruo a 1 mod 8 allora n è congruo a 1 mod 4 e congruo a 1 mod 2
	char count=0;//0 se non ci sono state divisioni 1 se ci sono state divisioni
	long j=0;//elementi dell'array divisibili per 2^k + salti di 2^k 

	long l,h,t,z;//interi per generare multipli di 2^k
	long indexv;//indice degli elementi nell'array divisibili per 2^k
	int p=2;
	char array_divided=0;
	int index_of_prime=2;

	mpz_t c,v;//coefficiente di a(j)
	mpz_t inverse_a,inverse2b;
	mpz_t r,r2,r3,r4;//radici quadrate modulo 2^k
	mpz_t ntemp,b_temp,p_to_k;
	mpz_t j1t,j2t,j3t,j4t;
	mpz_t l_temp,j_temp,index,h_temp,t_temp,z_temp;
	
	if(n==NULL || a==NULL || b==NULL || y==NULL || array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || M<=0 || mpz_sgn(n)<=0 || k<=0){
		handle_error_with_exit("invalid parameter_divide_all_by_2_to_k\n");
	}
	mpz_init(c);
	mpz_init(v);
	mpz_init(inverse_a);
	mpz_init(inverse2b);
	mpz_init(r);
	mpz_init(r2);
	mpz_init(r3);
	mpz_init(r4);
	mpz_init_set(ntemp,n);//ntemp=n
	mpz_init(b_temp);
	mpz_init(p_to_k);
	mpz_init(j1t);
	mpz_init(j2t);
	mpz_init(j3t);
	mpz_init(j4t);
	mpz_init(index);
	mpz_init(l_temp);
	mpz_init(j_temp);
	mpz_init(h_temp);
	mpz_init(t_temp);
	mpz_init(z_temp);
	mpz_ui_pow_ui(p_to_k,2,k);//p_to_k=2^k
	if(k==1){//radici quadrate di n modulo 2
		count=divide_all_by_2(array_of_number,M,matrix_factorization);//k=1 -> dividi per 2
			mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(r);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);
		return count;
	}
	else if(k==2){//radici quadrate di n modulo 4
		mpz_mod_ui(ntemp,ntemp,4);//n_temp=n mod 4
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 4
			handle_error_with_exit("n!=1 mod 4\n");//non esistono radici quadrate mod 4
			mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(r);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);
			return 0;
		}
		count=divide_all_by_4(array_of_number,M,matrix_factorization,n,a,b);//k=2 -> dividi per 4
		mpz_clear(c);
		mpz_clear(v);
		mpz_clear(inverse_a);
		mpz_clear(inverse2b);
		mpz_clear(r);
		mpz_clear(r2);
		mpz_clear(r3);
		mpz_clear(r4);
		mpz_clear(ntemp);//ntemp=n
		mpz_clear(b_temp);
		mpz_clear(p_to_k);
		mpz_clear(j1t);
		mpz_clear(j2t);
		mpz_clear(j3t);
		mpz_clear(j4t);
		mpz_clear(index);
		mpz_clear(l_temp);
		mpz_clear(j_temp);
		mpz_clear(h_temp);
		mpz_clear(t_temp);
		mpz_clear(z_temp);
		return count;
	}
	else if(k==3){//radici quadrate di n modulo 8
		mpz_mod_ui(ntemp,ntemp,8);
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 8
			handle_error_with_exit("n!=1 mod 8\n");//non esistono radici quadrate mod 8
			mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(r);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);
			return 0;
		}
		//radici quadrate note modulo 8
		mpz_set_si(r,1);//r=1;
		mpz_set_si(r2,3);//r2=3;
		mpz_set_si(r3,5);//r3=5;
		mpz_set_si(r4,7);//r4=7;
	}
	else{//k>=3
		mpz_mod_ui(ntemp,ntemp,8);
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 8
			handle_error_with_exit("n!=1 mod 8,k>=3\n");//non esistono radici quadrate mod 2^k con k>=3
			return 0;
		}
		//calcolo delle radici quadrate di n,vedi note
		mpz_set(r,y);//r=y
		mpz_mul_ui(r,r,2);//r=2y
		mpz_add_ui(r,r,1);//r=2*y+1==a

		mpz_neg(r2,r);//r2=-r==b
		mpz_set(r3,p_to_k);//r3=2^k
		if(mpz_divisible_ui_p(r3,2)==0){
				handle_error_with_exit("error in mpz_divisible all by 2 to k\n");
		}
		mpz_divexact_ui(r3,r3,2);//r3=2^(k-1)
		mpz_add_ui(r3,r3,1);//r3=(1+pow(2,k-1))
		mpz_mul(r3,r3,r);//r3=(1+pow(2,k-1))*r==c
		mpz_neg(r4,r3);//r4=-r3==d
		//riduzione modulo n delle radici
		mpz_mod(r,r,p_to_k);
		mpz_mod(r2,r2,p_to_k);
		mpz_mod(r3,r3,p_to_k);
		mpz_mod(r4,r4,p_to_k);
	}
	if(mpz_divisible_ui_p(a,2)==0){//se 2 non divide a calcolo i vari j,sono 4
		mpz_invert(inverse_a,a,p_to_k);//a_temp=a^-1 mod p^k
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_neg(j3t,b);//j3=-b
		mpz_neg(j4t,b);//j4=-b
		mpz_add(j1t,j1t,r);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
		mpz_add(j3t,j3t,r3);//j3=-b+r3
		mpz_add(j4t,j4t,r4);//j4=-b+r4
		mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
		mpz_mul(j3t,j3t,inverse_a);//j3=(-b+r)*inverse_a;
		mpz_mul(j4t,j4t,inverse_a);//j1=(-b+r)*inverse_a;
	}
	else{//2 divide a,esiste solo una soluzione,si usa solo j1
		handle_error_with_exit("2 is a factor of a\n");
		mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
		mpz_invert(inverse2b,b_temp,p_to_k);//(2b)^-1 mod p^k
		
		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 2 to k\n");
		}
		mpz_divexact(c,v,a);//c=(b^2-n)/a
		
		mpz_neg(j1t,c);//j1=-c
		mpz_mul(j1t,j1t,inverse2b);//j1=-c*inverse2b;
		j=1;
	}
	//riduzione dei j modulo 2^k
	mpz_mod(j1t,j1t,p_to_k);
	mpz_mod(j2t,j2t,p_to_k);
	mpz_mod(j3t,j3t,p_to_k);
	mpz_mod(j4t,j4t,p_to_k);
	l=0;
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0){//soluzioni di r1 con l positivo
		mpz_set(index,l_temp);//index=l*p_to_k
		mpz_add(index,index,j1t);//index=j1+l*p_to_k
		mpz_add_ui(index,index,M);//index=j1+l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division l positivo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k	
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0){//soluzioni di r1 con l negativo
		mpz_set(index,l_temp);//index=-l*p_to_k
		mpz_add(index,index,j1t);//index=j1-l*p_to_k
		mpz_add_ui(index,index,M);//index=j1-l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division l negativo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	}
	h=0;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,h_temp);//index=h*p_to_k
		mpz_add(index,index,j2t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division h positivo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo htemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	}
	
	h=1;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,h_temp);//index=-h*p_to_k
		mpz_add(index,index,j2t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division h negativo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo jtemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	}
	t=0;
	mpz_set_si(t_temp,t);//h_temp=h
	mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j3t,t_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,t_temp);//index=h*p_to_k
		mpz_add(index,index,j3t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division t positivo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		t++;
		//calcolo del nuovo htemp
		mpz_set_si(t_temp,t);//h_temp=h
		mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j3t,t_temp);//j_temp=j2t+h*p_to_k
	}
	t=1;
	mpz_set_si(t_temp,t);//h_temp=h
	mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(t_temp,t_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j3t,t_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,t_temp);//index=-h*p_to_k
		mpz_add(index,index,j3t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division t negativo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		t++;
		//calcolo del nuovo jtemp
		mpz_set_si(t_temp,t);//h_temp=h
		mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(t_temp,t_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j3t,t_temp);//j_temp=j2t-h*p_to_k
	}
	z=0;
	mpz_set_si(z_temp,z);//h_temp=h
	mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j4t,z_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,z_temp);//index=h*p_to_k
		mpz_add(index,index,j4t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z positivo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		z++;
		//calcolo del nuovo htemp
		mpz_set_si(z_temp,z);//h_temp=h
		mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j4t,z_temp);//j_temp=j2t+h*p_to_k
	}

	z=1;
	mpz_set_si(z_temp,z);//h_temp=h
	mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(z_temp,z_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j4t,z_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,z_temp);//index=-h*p_to_k
		mpz_add(index,index,j4t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		z++;
		//calcolo del nuovo jtemp
		mpz_set_si(z_temp,z);//h_temp=h
		mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(z_temp,z_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j4t,z_temp);//j_temp=j2t-h*p_to_k
	}
	mpz_clear(c);
	mpz_clear(v);
	mpz_clear(inverse_a);
	mpz_clear(inverse2b);
	mpz_clear(r);
	mpz_clear(r2);
	mpz_clear(r3);
	mpz_clear(r4);
	mpz_clear(ntemp);//ntemp=n
	mpz_clear(b_temp);
	mpz_clear(p_to_k);
	mpz_clear(j1t);
	mpz_clear(j2t);
	mpz_clear(j3t);
	mpz_clear(j4t);
	mpz_clear(index);
	mpz_clear(l_temp);
	mpz_clear(j_temp);
	mpz_clear(h_temp);
	mpz_clear(t_temp);
	mpz_clear(z_temp);
	if(array_divided==1){
		return 1;//ritorna 1 se c'è stata almeno 1 divisione
	}
	return 0;//ritorna 0 se non ci sono state divisioni nell'array
}*/
char divide_all_by_2_to_k_f(long M,struct thread_data thread_data,const mpz_t n,long k,const mpz_t a,const mpz_t b,const mpz_t y){
	//k esponente di 2 a,b coefficienti di a(j)=aj^2+2bj+c y coincide con y(k-2),se si sa y(k-2) allora a=2*y(k-2)+1

	//se n è congruo a 1 mod 8 allora n è congruo a 1 mod 4 e congruo a 1 mod 2
	char count=0;//0 se non ci sono state divisioni 1 se ci sono state divisioni
	/*long j=0;//elementi dell'array divisibili per 2^k + salti di 2^k

	long l,h,t,z;//interi per generare multipli di 2^k
	long indexv;//indice degli elementi nell'array divisibili per 2^k
	int p=2;
	char array_divided=0;
	int index_of_prime=2;

	mpz_t c,v;//coefficiente di a(j)
	mpz_t inverse_a,inverse2b;
	mpz_t rad,r2,r3,r4;//radici quadrate modulo 2^k
	mpz_t ntemp,b_temp,p_to_k;
	mpz_t j1t,j2t,j3t,j4t;
	mpz_t l_temp,j_temp,index,h_temp,t_temp,z_temp;*/
	
	if(n==NULL || a==NULL || b==NULL || y==NULL || M<=0 || mpz_sgn(n)<=0 || k<=0){
		handle_error_with_exit("invalid parameter_divide_all_by_2_to_k\n");
	}
	/*mpz_init(c);
	mpz_init(v);
	mpz_init(inverse_a);
	mpz_init(inverse2b);
	mpz_init(rad);
	mpz_init(r2);
	mpz_init(r3);
	mpz_init(r4);
	mpz_init_set(ntemp,n);//ntemp=n
	mpz_init(b_temp);
	mpz_init(p_to_k);
	mpz_init(j1t);
	mpz_init(j2t);
	mpz_init(j3t);
	mpz_init(j4t);
	mpz_init(index);
	mpz_init(l_temp);
	mpz_init(j_temp);
	mpz_init(h_temp);
	mpz_init(t_temp);
	mpz_init(z_temp);
	mpz_ui_pow_ui(p_to_k,2,k);//p_to_k=2^k
	 */
	if(k==1){//radici quadrate di n modulo 2
		count=divide_all_by_2_log(M,thread_data);//k=1 -> dividi per 2
			/*mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(rad);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);*/
		return count;
	}
	/*else if(k==2){//radici quadrate di n modulo 4
		mpz_mod_ui(ntemp,ntemp,4);//n_temp=n mod 4
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 4
			handle_error_with_exit("n!=1 mod 4\n");//non esistono radici quadrate mod 4
			mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(rad);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);
			return 0;
		}
		count=divide_all_by_4_f(M,mat,n,a,b);//k=2 -> dividi per 4
		mpz_clear(c);
		mpz_clear(v);
		mpz_clear(inverse_a);
		mpz_clear(inverse2b);
		mpz_clear(rad);
		mpz_clear(r2);
		mpz_clear(r3);
		mpz_clear(r4);
		mpz_clear(ntemp);//ntemp=n
		mpz_clear(b_temp);
		mpz_clear(p_to_k);
		mpz_clear(j1t);
		mpz_clear(j2t);
		mpz_clear(j3t);
		mpz_clear(j4t);
		mpz_clear(index);
		mpz_clear(l_temp);
		mpz_clear(j_temp);
		mpz_clear(h_temp);
		mpz_clear(t_temp);
		mpz_clear(z_temp);
		return count;
	}
	else if(k==3){//radici quadrate di n modulo 8
		mpz_mod_ui(ntemp,ntemp,8);
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 8
			handle_error_with_exit("n!=1 mod 8\n");//non esistono radici quadrate mod 8
			mpz_clear(c);
			mpz_clear(v);
			mpz_clear(inverse_a);
			mpz_clear(inverse2b);
			mpz_clear(rad);
			mpz_clear(r2);
			mpz_clear(r3);
			mpz_clear(r4);
			mpz_clear(ntemp);//ntemp=n
			mpz_clear(b_temp);
			mpz_clear(p_to_k);
			mpz_clear(j1t);
			mpz_clear(j2t);
			mpz_clear(j3t);
			mpz_clear(j4t);
			mpz_clear(index);
			mpz_clear(l_temp);
			mpz_clear(j_temp);
			mpz_clear(h_temp);
			mpz_clear(t_temp);
			mpz_clear(z_temp);
			return 0;
		}
		//radici quadrate note modulo 8
		mpz_set_si(rad,1);//r=1;
		mpz_set_si(r2,3);//r2=3;
		mpz_set_si(r3,5);//r3=5;
		mpz_set_si(r4,7);//r4=7;
	}
	else{//k>=3
		mpz_mod_ui(ntemp,ntemp,8);
		if(mpz_get_si(ntemp)!=1){//n non è congruo a 1 mod 8
			handle_error_with_exit("n!=1 mod 8,k>=3\n");//non esistono radici quadrate mod 2^k con k>=3
			return 0;
		}
		//calcolo delle radici quadrate di n,vedi note
		mpz_set(rad,y);//r=y
		mpz_mul_ui(rad,rad,2);//r=2y
		mpz_add_ui(rad,rad,1);//r=2*y+1==a

		mpz_neg(r2,rad);//r2=-r==b
		mpz_set(r3,p_to_k);//r3=2^k
		if(mpz_divisible_ui_p(r3,2)==0){
				handle_error_with_exit("error in mpz_divisible all by 2 to k\n");
		}
		mpz_divexact_ui(r3,r3,2);//r3=2^(k-1)
		mpz_add_ui(r3,r3,1);//r3=(1+pow(2,k-1))
		mpz_mul(r3,r3,rad);//r3=(1+pow(2,k-1))*r==c
		mpz_neg(r4,r3);//r4=-r3==d
		//riduzione modulo n delle radici
		mpz_mod(rad,rad,p_to_k);
		mpz_mod(r2,r2,p_to_k);
		mpz_mod(r3,r3,p_to_k);
		mpz_mod(r4,r4,p_to_k);
	}
	if(mpz_divisible_ui_p(a,2)==0){//se 2 non divide a calcolo i vari j,sono 4
		mpz_invert(inverse_a,a,p_to_k);//a_temp=a^-1 mod p^k
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_neg(j3t,b);//j3=-b
		mpz_neg(j4t,b);//j4=-b
		mpz_add(j1t,j1t,rad);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
		mpz_add(j3t,j3t,r3);//j3=-b+r3
		mpz_add(j4t,j4t,r4);//j4=-b+r4
		mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
		mpz_mul(j3t,j3t,inverse_a);//j3=(-b+r)*inverse_a;
		mpz_mul(j4t,j4t,inverse_a);//j1=(-b+r)*inverse_a;
	}
	else{//2 divide a,esiste solo una soluzione,si usa solo j1
		handle_error_with_exit("2 is a factor of a\n");
		mpz_mul_ui(b_temp,b_temp,2);//b_temp=2*b
		mpz_invert(inverse2b,b_temp,p_to_k);//(2b)^-1 mod p^k
		
		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){
				handle_error_with_exit("error in mpz_divisible divide all by 2 to k\n");
		}
		mpz_divexact(c,v,a);//c=(b^2-n)/a
		
		mpz_neg(j1t,c);//j1=-c
		mpz_mul(j1t,j1t,inverse2b);//j1=-c*inverse2b;
		j=1;
	}
	//riduzione dei j modulo 2^k
	mpz_mod(j1t,j1t,p_to_k);
	mpz_mod(j2t,j2t,p_to_k);
	mpz_mod(j3t,j3t,p_to_k);
	mpz_mod(j4t,j4t,p_to_k);
	l=0;
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0){//soluzioni di r1 con l positivo
		mpz_set(index,l_temp);//index=l*p_to_k
		mpz_add(index,index,j1t);//index=j1+l*p_to_k
		mpz_add_ui(index,index,M);//index=j1+l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log=mat->row[indexv].sum_log+r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		//mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k	
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0){//soluzioni di r1 con l negativo
		mpz_set(index,l_temp);//index=-l*p_to_k
		mpz_add(index,index,j1t);//index=j1-l*p_to_k
		mpz_add_ui(index,index,M);//index=j1-l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		//mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	}
	h=0;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,h_temp);//index=h*p_to_k
		mpz_add(index,index,j2t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo htemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	}
	
	h=1;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,h_temp);//index=-h*p_to_k
		mpz_add(index,index,j2t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo jtemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	}
	t=0;
	mpz_set_si(t_temp,t);//h_temp=h
	mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j3t,t_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,t_temp);//index=h*p_to_k
		mpz_add(index,index,j3t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		t++;
		//calcolo del nuovo htemp
		mpz_set_si(t_temp,t);//h_temp=h
		mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j3t,t_temp);//j_temp=j2t+h*p_to_k
	}
	t=1;
	mpz_set_si(t_temp,t);//h_temp=h
	mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(t_temp,t_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j3t,t_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,t_temp);//index=-h*p_to_k
		mpz_add(index,index,j3t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente  
		array_divided=1;//una divisione è stata effettuata
		t++;
		//calcolo del nuovo jtemp
		mpz_set_si(t_temp,t);//h_temp=h
		mpz_mul(t_temp,t_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(t_temp,t_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j3t,t_temp);//j_temp=j2t-h*p_to_k
	}
	z=0;
	mpz_set_si(z_temp,z);//h_temp=h
	mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j4t,z_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,z_temp);//index=h*p_to_k
		mpz_add(index,index,j4t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente
		array_divided=1;//una divisione è stata effettuata
		z++;
		//calcolo del nuovo htemp
		mpz_set_si(z_temp,z);//h_temp=h
		mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j4t,z_temp);//j_temp=j2t+h*p_to_k
	}

	z=1;
	mpz_set_si(z_temp,z);//h_temp=h
	mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(z_temp,z_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j4t,z_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,z_temp);//index=-h*p_to_k
		mpz_add(index,index,j4t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division z negativo\n");
		}
		mat->row[indexv].sum_log+=r.log_prime[1];//somma i log
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;//memorizza la fattorizzazione aggiungendo 1 all'esponente
		array_divided=1;//una divisione è stata effettuata
		z++;
		//calcolo del nuovo jtemp
		mpz_set_si(z_temp,z);//h_temp=h
		mpz_mul(z_temp,z_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(z_temp,z_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j4t,z_temp);//j_temp=j2t-h*p_to_k
	}
	mpz_clear(c);
	mpz_clear(v);
	mpz_clear(inverse_a);
	mpz_clear(inverse2b);
	mpz_clear(rad);
	mpz_clear(r2);
	mpz_clear(r3);
	mpz_clear(r4);
	mpz_clear(ntemp);//ntemp=n
	mpz_clear(b_temp);
	mpz_clear(p_to_k);
	mpz_clear(j1t);
	mpz_clear(j2t);
	mpz_clear(j3t);
	mpz_clear(j4t);
	mpz_clear(index);
	mpz_clear(l_temp);
	mpz_clear(j_temp);
	mpz_clear(h_temp);
	mpz_clear(t_temp);
	mpz_clear(z_temp);
	if(array_divided==1){
		return 1;//ritorna 1 se c'è stata almeno 1 divisione
	}*/
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

/*char divide_all_by_p_to_k(const mpz_t r,long p,int index_of_prime,long k,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,const mpz_t n,const mpz_t a,const mpz_t b){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p 

	//j1=r-x0;
	//j2=-r-x0;
	if((p<=1 && p!=-1 )|| r==NULL || n==NULL || a==NULL || b==NULL || k<=0 || array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL){
		handle_error_with_exit("invalid parameter divide all by p to k\n");
	}
	char array_divided=0;//0 se nessuna divisione effettuata,1 se sono state effettutate divisioni per p^k
	int l,h;
	long indexv;
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

	mpz_ui_pow_ui(p_to_k,p,k);//p_to_k=p^k
	mpz_sub(r2,p_to_k,r);//r2=p^k-r seconda radice quadrata di n modulo p^k
	long j=0;//indici dell'array divisibile per p^k,se j!=0 j2 non esiste

	if(mpz_divisible_ui_p(a,p)==0){//se p non divide a
		mpz_invert(inverse_a,a,p_to_k);//a_temp=a^-1 mod p^k
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add(j1t,j1t,r);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
		mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
	}
	else{//p divide a,esiste solo una soluzione,si usa solo j1

		//calcolo di c
		mpz_mul(v,b,b);//v=b^2
		mpz_sub(v,v,n);//v=b^2-n
		if(mpz_divisible_p(v,a)==0){
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
			calculate_root_poly_second_degree_mod_p_to_k(j1t,p_to_k,p,k,a,b,n);
			j=1;//solo 1 soluzione
		}
	}
	mpz_mod(j1t,j1t,p_to_k);
	mpz_mod(j2t,j2t,p_to_k);
	l=0;
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0){//soluzioni di r1 con l positivo
		mpz_set(index,l_temp);//index=l*p_to_k
		mpz_add(index,index,j1t);//index=j1+l*p_to_k
		mpz_add_ui(index,index,M);//index=j1+l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k	
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0){//soluzioni di r1 con l negativo
		mpz_set(index,l_temp);//index=-l*p_to_k
		mpz_add(index,index,j1t);//index=j1-l*p_to_k
		mpz_add_ui(index,index,M);//index=j1-l*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division\n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
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
	h=0;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,h_temp);//index=h*p_to_k
		mpz_add(index,index,j2t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo htemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	}

	h=1;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,h_temp);//index=-h*p_to_k
		mpz_add(index,index,j2t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		if(mpz_divisible_ui_p(array_of_number[indexv],p)==0){//se non è possibile dividerlo -> errore
			handle_error_with_exit("error in division \n");
		}
		mpz_divexact_ui(array_of_number[indexv],array_of_number[indexv],p);
		mpz_add_ui(matrix_factorization[indexv][index_of_prime+1],matrix_factorization[indexv][index_of_prime+1],1);//memorizza la fattorizzazione aggiungendo 1 all'esponente 
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo jtemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
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
}*/

/*char divide_all_by_p_to_k_f(const mpz_t rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b){//r=square root di n mod p^k,ritorna >zero se c'è stata almeno 1 divisione
//a(j)=aj^2+2bj+c,se polinomio forma semplice a=1 e b=x0
//(x0+j)^2==n mod p^k r:r^2==n mod p^k -> r=x0+j ->j=r-x0
//una volta trovate le soluzioni r1 e r2 e una volta trovati j e t le altre soluzioni sono della forma j+l*p^k t+h*p^k,cioè a salti di p^k rispetto a t e j,dove l ed h sono interi
//n.b j+l*p^k è sempre diverso da t+h*p^k per ogni l ed h interi,quindi non bisogna memorizzare quali elementi dell'array sono stati divisi per p 

	//j1=r-x0;
	//j2=-r-x0;
	if((p<=1 && p!=-1 )|| rad==NULL || n==NULL || a==NULL || b==NULL || k<=0){
		handle_error_with_exit("invalid parameter divide all by p to k\n");
	}
	char array_divided=0;//0 se nessuna divisione effettuata,1 se sono state effettutate divisioni per p^k
	int l,h;
	long indexv;
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

	mpz_ui_pow_ui(p_to_k,p,k);//p_to_k=p^k
	mpz_sub(r2,p_to_k,rad);//r2=p^k-r seconda radice quadrata di n modulo p^k
	long j=0;//indici dell'array divisibile per p^k,se j!=0 j2 non esiste

	if(mpz_divisible_ui_p(a,p)==0){//se p non divide a
		mpz_invert(inverse_a,a,p_to_k);//a_temp=a^-1 mod p^k
		mpz_neg(j1t,b);//j1=-b
		mpz_neg(j2t,b);//j2=-b
		mpz_add(j1t,j1t,rad);//j1=-b+r
		mpz_add(j2t,j2t,r2);//j2=-b+r2
		mpz_mul(j1t,j1t,inverse_a);//j1=(-b+r)*inverse_a;
		mpz_mul(j2t,j2t,inverse_a);//j2=(-b+r2)*inverse_a;
	}
	else{//p divide a,esiste solo una soluzione,si usa solo j1

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
			calculate_root_poly_second_degree_mod_p_to_k(j1t,p_to_k,p,k,a,b,n);
			j=1;//solo 1 soluzione
		}
	}
	mpz_mod(j1t,j1t,p_to_k);
	mpz_mod(j2t,j2t,p_to_k);
	l=0;
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0){//soluzioni di r1 con l positivo
		mpz_set(index,l_temp);//index=l*p_to_k
		mpz_add(index,index,j1t);//index=j1+l*p_to_k
		mpz_add_ui(index,index,M);//index=j1+l*p_to_k+M
		indexv=mpz_get_si(index);
		//if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			//handle_error_with_exit("error in division\n");
		//}

        //somma il logaritmo del primo index_of_prime nel numero all'indice indexv
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
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t+l*p_to_k	
	}
	l=1;//l=1 permette di non calcolare 2 volte la soluzione j+0*l
	mpz_set_si(l_temp,l);//l_temp=l
	mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
	mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
	mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0){//soluzioni di r1 con l negativo
		mpz_set(index,l_temp);//index=-l*p_to_k
		mpz_add(index,index,j1t);//index=j1-l*p_to_k
		mpz_add_ui(index,index,M);//index=j1-l*p_to_k+M
		indexv=mpz_get_si(index);
		//if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			//handle_error_with_exit("error in division\n");
		//}
        //somma il logaritmo del primo index_of_prime nel numero all'indice indexv
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
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		array_divided=1;//una divisione è stata effettuata
		l++;
		//calcolo del nuovo jtemp
		mpz_set_si(l_temp,l);//l_temp=l
		mpz_mul(l_temp,l_temp,p_to_k);//l_temp=l*p_to_k
		mpz_neg(l_temp,l_temp);//l_temp=-l*p_to_k
		mpz_add(j_temp,j1t,l_temp);//j_temp=j1t-l*p_to_k
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
	h=0;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	while(mpz_cmp_si(j_temp,M)<=0 && j==0){//r2+p^k,soluzioni di r2 con h positivo
		mpz_set(index,h_temp);//index=h*p_to_k
		mpz_add(index,index,j2t);//index=j2+h*p_to_k
		mpz_add_ui(index,index,M);//index=j2+h*p_to_k+M
		indexv=mpz_get_si(index);
		//if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
			//handle_error_with_exit("error in division\n");
		//}
        //somma il logaritmo del primo index_of_prime nel numero all'indice indexv
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
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo htemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t+h*p_to_k
	}

	h=1;
	mpz_set_si(h_temp,h);//h_temp=h
	mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
	mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
	mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
	while(mpz_cmp_si(j_temp,-M)>=0 && j==0){//soluzioni di r2 con h negativo
		mpz_set(index,h_temp);//index=-h*p_to_k
		mpz_add(index,index,j2t);//index=j2-h*p_to_k
		mpz_add_ui(index,index,M);//index=j2-h*p_to_k+M
		indexv=mpz_get_si(index);
		//if(mpz_divisible_ui_p(mat->row[indexv].num,p)==0){//se non è possibile dividerlo -> errore
		//	handle_error_with_exit("error in division\n");
		//}
        //somma il logaritmo del primo index_of_prime nel numero all'indice indexv
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
		//mpz_divexact_ui(mat.row[indexv].num,mat.row[indexv].num,p);
		//mat.row[indexv].factorization[index_of_prime+1]+=1;
		array_divided=1;//una divisione è stata effettuata
		h++;
		//calcolo del nuovo jtemp
		mpz_set_si(h_temp,h);//h_temp=h
		mpz_mul(h_temp,h_temp,p_to_k);//h_temp=h*p_to_k
		mpz_neg(h_temp,h_temp);//h_temp=-h*p_to_k
		mpz_add(j_temp,j2t,h_temp);//j_temp=j2t-h*p_to_k
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
 */
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
/*void factor_matrix(const mpz_t n,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,struct node*head_f_base,int cardinality_factor_base,const mpz_t a,const mpz_t b){
	if(n==NULL || a==NULL || b==NULL ||  mpz_sgn(n)<=0 || M<=0 || array_of_number==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || head_f_base==NULL || cardinality_factor_base<=0){
		handle_error_with_exit("error in factor matrix\n");
	}
	long p;
	char count,t;
	mpz_t y,n_temp,rootpk,p_temp,r1;//r1 square root mod p,r==square root mod p^k
	int k;
	
	mpz_init(y);
	mpz_init(n_temp);
	mpz_init(rootpk);
	mpz_init(p_temp);
	mpz_init(r1);
	struct node*ptr_list=head_f_base;
	for (int i=0;i<cardinality_factor_base;i++){//per ogni elemento della factor base
		count=1;
		mpz_set(p_temp,ptr_list->prime);//p_temp=primo iesimo della factor base
		ptr_list=ptr_list->next;
		p=mpz_get_si(p_temp);//p=p_temp
		if(p==-1){
			divide_all_by_min1(array_of_number,M,matrix_factorization);//divide la matrice di fattorizzazione per meno 1
			continue;
		}
		if(p==2){
			mpz_set_si(y,0);//y1=0,necessario per calcolare le radici quadrate di n modulo 2^k
			k=1;
			while(count==1){
				if(k>=4){
					calculate_new_y(y,n,k);//k=4 calcolo-> y2 da y1,k=5 calcolo y3 da y2 ecc
				}
				count=divide_all_by_2_to_k(array_of_number,M,matrix_factorization,n,k,a,b,y);
				k++;
			}
		}
		else{//p>2 e dispari
			mpz_set(n_temp,n);//n_temp=n
			mpz_set_si(p_temp,p);//p_temp=p
			mpz_mod(n_temp,n_temp,p_temp);//n_temp = n mod p
			t=quadratic_residue(r1,n_temp,p_temp);//r1 radice quadrata di n mod p
			if(t==-1 || t==0){
				handle_error_with_exit("error in calculate quadratic_residue\n");
			}
			k=1;
			while(count==1){//finquando c'è stata almeno 1 divisione continua a dividere
				if(k!=1){
					square_root_mod_p_to_k(rootpk,r1,p,n,k);//rootpk ottieni la radice modulo p alla k,con k>1
				}
				else{
					mpz_set(rootpk,r1);//rootpk=r1
				}
				count=divide_all_by_p_to_k(rootpk,p,2*i,k,array_of_number,M,matrix_factorization,n,a,b);
				k++;//passa alla prossima potenza
			}	
		}
	}
	
	mpz_clear(y);
	mpz_clear(n_temp);
	mpz_clear(rootpk);
	mpz_clear(p_temp);
	mpz_clear(r1);
	return;
}*/

/*void factor_matrix_f(const mpz_t n,long M,struct thread_data thread_data,int cardinality_factor_base,const mpz_t a){
	if(n==NULL || a==NULL || mpz_sgn(n)<=0 || M<=0 || cardinality_factor_base<=0){
		handle_error_with_exit("error in factor matrix_f\n");
	}
	long p;
	char count,t;
	mpz_t y,n_temp,rootpk,p_temp,r1;//r1 square root mod p,r==square root mod p^k
	int k;
	
	mpz_init(y);
	mpz_init(n_temp);
	mpz_init(rootpk);
	mpz_init(p_temp);
	mpz_init(r1);
	for (int i=0;i<cardinality_factor_base;i++){//per ogni elemento della factor base
		count=1;
		p=r.prime[i];//primo iesimo della factor base
		if(p==-1){
			//divide_all_by_min1_f(M,thread_data);//divide la matrice di fattorizzazione per meno 1
			continue;
		}
		if(p==2){
			mpz_set_si(y,0);//y1=0,necessario per calcolare le radici quadrate di n modulo 2^k
			k=1;
			while(count==1){
				//if(k>=4){
				//	calculate_new_y(y,n,k);//k=4 calcolo-> y2 da y1,k=5 calcolo y3 da y2 ecc
				//}
				count=divide_all_by_2_to_k_f(M,thread_data,n,k,a,thread_data.b,y);
				//k++;
				break;
			}
		}
		else{//p>2 e dispari
			mpz_set(n_temp,n);//n_temp=n
			mpz_set_si(p_temp,p);//p_temp=p
			mpz_mod(n_temp,n_temp,p_temp);//n_temp = n mod p
			t=quadratic_residue(r1,n_temp,p_temp);//r1=radice quadrata di n mod p
			if(t==-1 || t==0){
				handle_error_with_exit("error in calculate quadratic_residue\n");
			}
			k=1;
			while(count==1){//finquando c'è stata almeno 1 divisione continua a dividere
				if(k!=1){
					square_root_mod_p_to_k(rootpk,r1,p,n,k);//rootpk ottieni la radice modulo p alla k,con k>1
				}
				else{
					mpz_set(rootpk,r1);//rootpk=r1
				}
				count=divide_all_by_p_to_k_f(rootpk,p,i,k,M,thread_data,n,a,thread_data.b);
				//k++;//passa alla prossima potenza
				break;
			}	
		}
	}
	mpz_clear(y);
	mpz_clear(n_temp);
	mpz_clear(rootpk);
	mpz_clear(p_temp);
	mpz_clear(r1);
	return;
}*/
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
/*void factor_matrix_f(const mpz_t n,long M,struct thread_data *thread_data,int cardinality_factor_base,const mpz_t a,
		struct a_struct*array_a_struct,int s){
    if(n==NULL || a==NULL || mpz_sgn(n)<=0 || M<=0 || cardinality_factor_base<=0 || array_a_struct==NULL || s<0){
        handle_error_with_exit("error in factor matrix_f\n");
    }
    pthread_t *array_tid=NULL;
    struct factorization_thread_data*factorization_thread_data=NULL;
    int num_thread;
    if(NUM_THREAD_FACTORIZATION>cardinality_factor_base){
    	num_thread=cardinality_factor_base;
    }
    else{
    	num_thread=NUM_THREAD_FACTORIZATION;
    }
    array_tid = alloc_array_tid(num_thread);//alloca memoria per contenere tutti i tid
    factorization_thread_data = create_factorization_threads(array_tid,*thread_data,a,num_thread);//crea tutti i thread
    join_all_threads(array_tid,num_thread);//aspetta tutti i thread
    if (array_tid != NULL) {//libera memoria allocata
        free(array_tid);
        array_tid = NULL;
    }
	memcpy((*thread_data).numbers,factorization_thread_data[0].thread_data.numbers,sizeof(struct number)*(2*M+1));
	mpz_clear(factorization_thread_data[0].thread_data.b);
	factorization_thread_data[0].thread_data.head=NULL;
	factorization_thread_data[0].thread_data.tail=NULL;
	free(factorization_thread_data[0].thread_data.numbers);

	for(int i=1;i<num_thread;i++){
    	for(int j=0;j<2*M+1;j++){
			(*thread_data).numbers[j].sum_log+=factorization_thread_data[i].thread_data.numbers[j].sum_log;
			if((*thread_data).numbers[j].first_index_f_base==-1 && factorization_thread_data[i].thread_data.numbers[j].first_index_f_base!=-1){
				(*thread_data).numbers[j].first_index_f_base=factorization_thread_data[i].thread_data.numbers[j].first_index_f_base;
			}
			if((*thread_data).numbers[j].first_index_f_base>factorization_thread_data[i].thread_data.numbers[j].first_index_f_base && factorization_thread_data[i].thread_data.numbers[j].first_index_f_base!=-1){
				(*thread_data).numbers[j].first_index_f_base=factorization_thread_data[i].thread_data.numbers[j].first_index_f_base;
			}
			if((*thread_data).numbers[j].last_index_f_base<factorization_thread_data[i].thread_data.numbers[j].last_index_f_base){
				(*thread_data).numbers[j].last_index_f_base=factorization_thread_data[i].thread_data.numbers[j].last_index_f_base;
			}
    	}
		mpz_clear(factorization_thread_data[i].thread_data.b);
		factorization_thread_data[i].thread_data.head=NULL;
		factorization_thread_data[i].thread_data.tail=NULL;
		free(factorization_thread_data[i].thread_data.numbers);
    }
    if (factorization_thread_data != NULL) {
	    free(factorization_thread_data[0].pointer);
        free(factorization_thread_data);
        factorization_thread_data = NULL;
    }
    return;
}*/

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
/*int count_number_B_smooth_in_matrix(mpz_t**potential_matrix_B_smooth,int card_factor_base,long size){
//conta il numero di numeri B_smooth nell'array dei numeri e inserisce i valori sia nella lista degli indici B_smooth,sia nella lista delle radici quadrate dei numeri B_smoot
	mpz_t temp,i_temp;

	mpz_init(temp);
	mpz_init(i_temp);
	if(potential_matrix_B_smooth==NULL  || *potential_matrix_B_smooth==NULL || card_factor_base<=0 || size<=0){
		handle_error_with_exit("error in parameter count number_B_smooth\n");
	}
	int count=0;
	int index_residuos=card_factor_base*2+1;
	int index_square=card_factor_base*2;
	if(mpz_cmp_si(potential_matrix_B_smooth[0][index_residuos],1)!=0){
			mpz_clear(temp);
			mpz_clear(i_temp);
			return count;//0 numeri b_smooth trovati
	}
	for(long i=0;i<size-1;i++){//controlla tutte le righe tranne l'ultima
		if(mpz_cmp_si(potential_matrix_B_smooth[i][index_residuos],1)==0){//il residuo è uguale a 1
			if(mpz_cmp(potential_matrix_B_smooth[i][index_square],potential_matrix_B_smooth[i+1][index_square])!=0){
				//i quadrati sono diversi
				count++;
			}
			else{//i quadrati sono uguali
				mpz_set_si(potential_matrix_B_smooth[i][index_residuos],0);//imposta residuo della prima riga considerata a 0
			}
		}
		else{
			handle_error_with_exit("error in count_number_B_smooth in matrix,residuos!=1 founded\n");
		}
	}
	count++;//l'ultimo quadrato è sicuramente diverso dagli altri
	mpz_clear(temp);
	mpz_clear(i_temp);
	return count;
}*/

/*mpz_t**save_element_B_smooth_in_matrix(mpz_t**matrix_factorization,long size,int cardinality_f_base,int* num_of_B_smooth){//alloca matrice con la fattorizzazione dei numeri b_smmoth
	if(size<0 || cardinality_f_base<=0 || (size>0 && matrix_factorization==NULL) || 
		(size>0 && *matrix_factorization==NULL) || num_of_B_smooth==NULL){
		handle_error_with_exit("invalid parameter\n");
	}
	mpz_t**matrix_B_smooth;
	mpz_t index_temp;
	long num_B_smooth_left;
	long index_residuos=cardinality_f_base*2+1;
	if(size==0){
		return NULL;
	}
	*num_of_B_smooth=count_number_B_smooth_in_matrix(matrix_factorization,cardinality_f_base,size);
	if(*num_of_B_smooth<=cardinality_f_base || *num_of_B_smooth-cardinality_f_base<THRESOLD_RELATION){
		return NULL;//ritorna NULL,questa matrice delle fattorizzazioni non va bene per fattorizzare n
	}
	mpz_init(index_temp);
	matrix_B_smooth=alloc_matrix_mpz(*num_of_B_smooth,cardinality_f_base*2+1);//l'1 in più serve per memorizzare i quadrati
	int j=0;
	num_B_smooth_left=*num_of_B_smooth;
	while(num_B_smooth_left>0){
		if(mpz_cmp_si(matrix_factorization[j][index_residuos],1)==0){//se residuo è 1
			copy_array_mpz(matrix_B_smooth[*num_of_B_smooth-num_B_smooth_left],matrix_factorization[j],cardinality_f_base*2+1);
			j++;//passa alla prossima riga
			num_B_smooth_left--;
		}
		else{//residuo è 0
			j++;//passa alla prossima riga
		}
	}
	mpz_clear(index_temp);
	return matrix_B_smooth;
}*/

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

int rand_int(int a,int b){
	if(a<0 || b<0){
		handle_error_with_exit("error in rand_int\n");
	}
	if(a>b){
		return ((rand()%(a-b+1))+b);
	}
	else if(b>a){
		return ((rand()%(b-a+1))+a);
	}
	else{
		return a;
	}
}
int rand_long(long a,long b){
	if(a<0 || b<0){
		handle_error_with_exit("error in rand_long\n");
	}
	if(a>b){
		return ((rand()%(a-b+1))+b);
	}
	else if(b>a){
		return ((rand()%(b-a+1))+a);
	}
	else{
		return a;
	}
}

long floor_sqrt(long x){//cerca con ricerca binaria la radice intera di n
	long start = 1, end = x/2, ans,mid;  
	// Base cases 
	if(x<0){
		handle_error_with_exit("impossibile trovare radice di un numero negativo\n");
	}
    if (x == 0 || x == 1){
       return x;
 	}
    while (start <= end) {        
        mid = (start + end) / 2;
        // If x is a perfect square
        if (mid*mid == x){
            return mid;
 	}
        // Since we need floor, we update answer when mid*mid is 
        // smaller than x, and move closer to sqrt(x)
        if (mid*mid < x) {
            start = mid + 1;
            ans = mid;
        } 
        else{// If mid*mid is greater than x
            end = mid - 1;        
	}
    }
    return ans;
}


