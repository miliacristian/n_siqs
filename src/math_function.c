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




