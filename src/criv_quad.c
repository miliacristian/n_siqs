#include "criv_quad.h"

extern long k;
extern struct timespec timer;
extern struct timespec time_start;
extern int num_increment_M_and_B;
extern struct row_factorization r;
extern double thresold_relation;
extern char combined;


void adjust_array_bi(mpz_t *array_bi,int s,const mpz_t a){//if 2*b>a ->b=a-b
	if(array_bi==NULL || s<0 || a==NULL){
		handle_error_with_exit("error in adjust array_bi\n");
	}
	mpz_t temp;
	mpz_init(temp);
	long length=pow(2,s-1);
	for(long i=0;i<length;i++){
		mpz_mul_ui(temp,array_bi[i],2);//temp=2*b[i]
		if(mpz_cmp(temp,a)>0){//if 2*b[1]>a
			mpz_sub(array_bi[i],a,array_bi[i]);//b[1]=a-b[1]
		}
	}
	mpz_clear(temp);
	return;
}

void calculate_thresold_a(mpfr_t thresold_a,const mpz_t n,long M){
	if(M<=0 || thresold_a==NULL || n==NULL){
		handle_error_with_exit("error in calculate thresold_a\n");
	}
	mpfr_t rad2,radnn,m;
	mpfr_t t;

	mpfr_init(t);
	mpfr_init(radnn);
	mpfr_init(rad2);
	mpfr_init(m);

	mpfr_set_si(m,M,MPFR_RNDN);//m=M
	mpfr_set_d(rad2,RAD2,MPFR_RNDN);//rad2=1.4142
	mpfr_set_z(radnn,n,MPFR_RNDN);//radn=n
    mpfr_sqrt(radnn,radnn,MPFR_RNDN);//radn=rad(n)
	mpfr_mul(t,radnn,rad2,MPFR_RNDN);//t=rad(n)*rad2
	mpfr_div(t,t,m,MPFR_RNDN);//t2=rad(n)*rad2/M

	mpfr_set(thresold_a,t,MPFR_RNDN);
	mpfr_clear(t);
	mpfr_clear(rad2);
	mpfr_clear(radnn);
	mpfr_clear(m);
	return;
}

void calculate_x0(mpz_t x0,const mpz_t n,int k,char *factorized){
	if(x0 ==NULL || n==NULL || (k!=1 && k!=3 && k!=5 && k!=7) || factorized==NULL){
		handle_error_with_exit("error in calculate x0\n");
	}
	mpz_t s,n_temp;
	mpz_init(s);
	mpz_init(n_temp);

	if(mpz_divisible_ui_p(n,k)==0){
		handle_error_with_exit("error in mpz_divisible calculate x0,n not divisible by k\n");
	}
	mpz_divexact_ui(n_temp,n,k);//n_temp=n/k
	mpz_sqrt(x0,n_temp);//x0=rad(n)
	mpz_mul(s,x0,x0);//s=x0^2
	if(mpz_cmp(s,n_temp)==0){//se s=x0^2=n_temp
		gmp_printf("factorization of n=%Zd*%Zd\n",x0,x0);
		/*fprintf(file_log,"p=");
		mpz_out_str(file_log,10,x0);
		fprintf(file_log," ");
		fprintf(file_log,"q=");
		mpz_out_str(file_log,10,x0);
		fprintf(file_log," ");*/
		*factorized=1;
		mpz_clear(s);
		mpz_clear(n_temp);
		return;
	}

	mpz_sqrt(x0,n);//x0=rad(n)
	mpz_mul(s,x0,x0);//s=x0^2
	if(mpz_cmp(s,n)==0){//se s=x0^2=n
		if(mpz_divisible_ui_p(n,k*k)==0){
			handle_error_with_exit("error in calculate x0,n not divisible by k^2\n");
		}
		mpz_divexact_ui(s,n,k*k);//s=/k
		mpz_sqrt(s,s);
		gmp_printf("factorization of n=%Zd*%Zd*%d\n",s,s,k);
		/*fprintf(file_log,"p=");
		mpz_out_str(file_log,10,s);
		fprintf(file_log," ");
		fprintf(file_log,"q=");
		mpz_out_str(file_log,10,s);
		fprintf(file_log," ");
		fprintf(file_log,"z=%d ",k);*/
		*factorized=1;
		mpz_clear(s);
		mpz_clear(n_temp);
		return;
	}
	
	mpz_clear(s);
	mpz_clear(n_temp);
	return;
}




void calculate_p_min_p_max_i_f(long*p_min_i,long*p_max_i,struct node_factor_base*head_f_base_f,long cardinality_factor_base){
	if(p_min_i==NULL || p_max_i==NULL || head_f_base_f==NULL || cardinality_factor_base<=0){
		handle_error_with_exit("error in calculate _p_min_p_max_i\n");
	}
	struct node_factor_base*list=head_f_base_f;
	mpz_t temp;
	mpz_init(temp);
	long i=0;
	while(list!=NULL){
		mpz_set_si(temp,list->prime);//valore del primo della factor base
		if(*p_min_i==0 && mpz_cmp_si(temp,SIQS_MIN_PRIME_POLYNOMIAL)>=0){
			*p_min_i=i;
		}
		if(*p_max_i==0 && mpz_cmp_si(temp,SIQS_MAX_PRIME_POLYNOMIAL)>0){
			*p_max_i=i-1;
			break;
		}
		list=list->next;//vai al prossimo elemento
		i++;//incrementa indice
	}
	if(*p_max_i==0){//primi della factor base<SIQS_MAX_PRIME_POLYNOMIAL
		*p_max_i=cardinality_factor_base-1;
	}
	if(*p_min_i<=1){//primi della factor base<SIQS_MIN_PRIME_POLYNOMIAL
		*p_min_i=2;//salta -1 e 2
	}
	if((*p_max_i-*p_min_i)<20){//se la distanza tra p_max e p_min è minore di 20(s_max=20) allora imposta p_min al minimo tra 2 e 5
		*p_min_i=min(*p_min_i,5);
	}
	mpz_clear(temp);
	return;
}

void calculate_target_a1_f(mpfr_t target_a1,const mpfr_t target_a,struct node_factor_base*head_f_base_f,
		long p_min_i,long p_max_i,int 	cardinality_factor_base){
	if(target_a==NULL || head_f_base_f==NULL || p_min_i<0 || p_max_i<=0 || p_min_i>p_max_i || target_a1==NULL || 
	cardinality_factor_base<=p_max_i){
		handle_error_with_exit("error in calculate_target_a1\n");
	}
	mpz_t prime_min;
	mpz_t prime_max;
	mpz_t prime_sum;
	mpfr_t prime_avg;
	mpfr_t temp;
	mpfr_t sqrt_prime_avg;

	mpz_init(prime_sum);
	mpz_init(prime_min);
	mpz_init(prime_max);
	mpfr_init(prime_avg);
	mpfr_init(sqrt_prime_avg);
	mpfr_init(temp);
	int count=0;
	struct node_factor_base*p=head_f_base_f;
	while(count<p_min_i){
		p=p->next;
		count++;
	}
	if(count!=p_min_i){
		handle_error_with_exit("error in calculate target_a1\n");
	}

	mpz_set_si(prime_min,p->prime);//ottieni prime_min
	while(count<p_max_i){
		p=p->next;
		count++;
	}
	if(count!=p_max_i){
		handle_error_with_exit("error in calculate target_a1\n");
	}
	mpz_set_si(prime_max,p->prime);//ottieni prime_max

	mpz_add(prime_sum,prime_max,prime_min);//prime_sum=prime_min+prime_max
	mpfr_set_z(prime_avg,prime_sum,MPFR_RNDN);//prime_avg=prime_sum
	mpfr_div_2ui(prime_avg,prime_avg,1,MPFR_RNDN);//prime_avg=(pmax+pmin)/2
	mpfr_sqrt(sqrt_prime_avg,prime_avg,MPFR_RNDN);//sqrt_prime=rad((pmax+pmin)/2)
	mpfr_set(temp,target_a,MPFR_RNDN);//temp=target=rad(2*n)/M
	mpfr_div(temp,temp,sqrt_prime_avg,MPFR_RNDN);//temp=target/rad((pmax+pmin)/2)
	mpfr_set(target_a1,temp,MPFR_RNDN);//target_a1=temp=target/rad((pmax+pmin)/2)
	mpz_clear(prime_min);
	mpz_clear(prime_max);
	mpz_clear(prime_sum);
	mpfr_clear(prime_avg);
	mpfr_clear(sqrt_prime_avg);
	mpfr_clear(temp);
	return;
}

void calculate_a_f2(mpz_t a,const mpfr_t target_a,int*s,struct node_factor_base*head_f_base_f,long cardinality_factor_base,int**best_q,int**best_q_number){
	if(s==NULL || target_a==NULL || mpfr_sgn(target_a)<0 || head_f_base_f==NULL || cardinality_factor_base<=0 || a==NULL || best_q==NULL || best_q_number==NULL){
		handle_error_with_exit("error in calculate_a\n");
	}
	long p_min_i=0;//indice del minimo primo da scegliere rispetto alla factor base
	long p_max_i=0;//indice del massimo primo da scegliere rispetto alla factor base
	//p_max_i>p_min_i
	char initialized=0;//se è zero la prima volta che si trova un a,anche se il rapporto a/target è pessimo si
	        //prende quell'a
	int value;
	int p_temp2;
	long p_i;//indici dei primi scelti nella factor base per rappresentare a
	long s_max;//valore massimo di s=pmax-pmin+1,numero massimo di primi da scegliere
	mpz_t v;
	int iter=0,iter2=0;
	mpz_init(v);
	if(cardinality_factor_base<3){
		mpz_set_si(a,0);//poni a=0
		*s=0;
		mpz_clear(v);
		*best_q=NULL;
		*best_q_number=NULL;
		return;//ritorna array_of_prime_chosen_for_a==NULL e a=0
	}
	if(cardinality_factor_base==3){
		get_element_linked_list_f(&value,head_f_base_f,2);//a=valore del 3 primo della factor base
		mpz_set_si(a,value);
		*s=1;
		//fprintf(file_log,"p_min=p_max=2 ");
		*best_q_number=alloc_array_int(1);//array che contiene i dei primi scelti
		*best_q=alloc_array_int(1);	
		(*best_q_number)[0]=value;
		(*best_q)[0]=2;
		mpz_clear(v);
		return;
	}
	mpz_t p_temp;
	mpfr_t target_a1;
	mpz_t best_a;
	mpfr_t best_ratio;
	double ratio_double=0,best_ratio_double=0;
	mpfr_t ratio;
	mpz_t a2_int;
	mpfr_t a2_double;
	mpfr_t p_rational;
	double dist_best,dist_ratio;
	mpfr_init(p_rational);
	mpz_init(a2_int);
	mpfr_init(a2_double);
	mpz_init(p_temp);
	mpfr_init(target_a1);

	mpz_init(best_a);
	mpfr_init(best_ratio);//double
	mpfr_init(ratio);//double


	mpz_set_si(best_a,0);//best_a=0
	mpfr_set_si(best_ratio,0,MPFR_RNDN);//best_target=0

	calculate_p_min_p_max_i_f(&p_min_i,&p_max_i,head_f_base_f,cardinality_factor_base);
	s_max=p_max_i-p_min_i+1;//nel caso peggiore s=p_max-p_min+1
	if(s_max>S_MAX){
		s_max=S_MAX;
	}
	calculate_target_a1_f(target_a1,target_a,head_f_base_f,p_min_i,p_max_i,cardinality_factor_base);
	if(mpfr_cmp_si(target_a1,1)<=0){//se target_a1<0 poni s=0 e ritorna
		mpz_set_si(a,0);//poni a=0
		*s=0;
		mpz_clear(v);
		mpfr_clear(p_rational);
		mpz_clear(a2_int);
		mpfr_clear(a2_double);
		mpz_clear(p_temp);
		mpfr_clear(target_a1);
		mpz_clear(best_a);
		mpfr_clear(best_ratio);//double
		mpfr_clear(ratio);//double
		return;
	}
	//printf("target_a1=");
	//mpfr_out_str(stdout,10,0,target_a1,MPFR_RNDN);
	//printf("\n");
    //print_time_elapsed("time to calculate target_a1");

	int count=0;//conta quante volte while(mpz_cmp(a,target_a1)<0) è verificata
	int*q=alloc_array_int(s_max);//array che contiene gli indici dei primi scelti
	int*q_number=alloc_array_int(s_max);//array che contiene i dei primi scelti
	*best_q_number=alloc_array_int(s_max);//array che contiene i numeri dei primi scelti
	*best_q=alloc_array_int(s_max);//array che contiene gli indici dei primi scelti
	int length_best_q=0;//inizialmente la lunghezza della lista migliore è zero
	for(int i=0;i<NUM_ITER_FOR_CALCULATE_A;i++){//iterazioni per cercare di migliorare a
		mpz_set_si(a2_int,1);//a2=1 azzera a2
		memset(q,0,sizeof(int)*s_max);//azzera q
		memset(q_number,0,sizeof(int)*s_max);//azzera q
		count=0;//azzera le append a q
		iter=0;
		mpfr_set_z(a2_double,a2_int,MPFR_RNDN);
		while(mpfr_cmp(a2_double,target_a1)<0 && iter<=MAX_ITER){//finquando a2<target_a1 oppure sono state raggiunte tot iterazioni
			p_i=2;
			iter2=0;
			while((p_i==2 || is_in_array_int(q,s_max,p_i)) && iter2<=MAX_ITER2){//scegli un p_i che non è stato ancora scelto
				p_i=rand_long(p_min_i,p_max_i);
				if(p_i>p_max_i || p_i<p_min_i){
					handle_error_with_exit("error random long generator\n");
				}
				iter2++;
			}
			//numero diverso da 2 e che non è nell'array di q
			if(iter2>=MAX_ITER2){
				break;//ricomincia da capo e setta a2=1
			}
			get_element_linked_list_f(&p_temp2,head_f_base_f,p_i);//prendi l'elemento p_i-esimo dalla factor base
			mpz_set_si(p_temp,p_temp2);
			mpfr_set_z(p_rational,p_temp,MPFR_RNDN);//p_rational=p_temp
			mpfr_mul(a2_double,a2_double,p_rational,MPFR_RNDN);//a2=a2*p
			mpz_mul(a2_int,a2_int,p_temp);//a2=a2*p

			if(count>=s_max){
				handle_error_with_exit("error in calculate a,index out of bounds\n");
			}
			q_number[count]=mpz_get_si(p_temp);
			q[count]=p_i;
			count++;//aumenta il conto delle append fatte
			iter++;
			if(count>=s_max){
				break;
			}
		}
		//una sequenza di numeri è stata scelta vediamo quanto è buona(vedendo i rapporti della migliore attualmente),se è buona la 			aggiorniamo
		mpfr_div(ratio,a2_double,target_a,MPFR_RNDN);//ratio=a2/target_a

		ratio_double=mpfr_get_d(ratio,MPFR_RNDN);
		if(initialized==0){
		            initialized=1;
                    mpz_set(best_a,a2_int);//best_a=a2
                    mpfr_set(best_ratio,ratio,MPFR_RNDN);//best_ratio=ratio
                    best_ratio_double=ratio_double;
                    memcpy(*best_q,q,sizeof(int)*s_max);//best_q=q
                    memcpy(*best_q_number,q_number,sizeof(int)*s_max);//best_q=q
                    length_best_q=count;//la lunghezza di best_q coincide con il numero di append fatte
		}
		dist_best=fabs(RATIO_A-best_ratio_double);//distanza di best_ratio da 1
        	dist_ratio=fabs(RATIO_A-ratio_double);//distanza di ratio_da 1
		if(dist_ratio<dist_best){
		  mpz_set(best_a,a2_int);//best_a=a2
		  mpfr_set(best_ratio,ratio,MPFR_RNDN);//best_ratio=ratio
		  best_ratio_double=ratio_double;
		  memcpy(*best_q,q,sizeof(int)*s_max);//best_q=q
		  memcpy(*best_q_number,q_number,sizeof(int)*s_max);//best_q=q
		  length_best_q=count;//la lunghezza di best_q coincide con il numero di append fatte
		}
	}

	mpz_set_si(a,1);//a=1
	for(int i=0;i<length_best_q;i++){//moltiplica tutti i fattori di a
		mpz_mul_si(a,a,(*best_q_number)[i]);
	}
	if(mpz_cmp(best_a,a)!=0){
		gmp_printf("a=%Zd\n",a);
		gmp_printf("best_a=%Zd\n",best_a);
		handle_error_with_exit("error in calculate a\n");
	}
	*s=length_best_q;//imposta il valore di s
	if(*s>s_max){
		handle_error_with_exit("error in calculate a,s must be minor or equal to s_max\n");
	}
	if(*s==0 && mpz_cmp(a,0)!=0){
		handle_error_with_exit("error in function calculate_a s or a\n");
	}
	free(q);
	free(q_number);
	mpfr_clear(p_rational);
	mpz_clear(p_temp);
	mpfr_clear(target_a1);
	mpfr_clear(a2_double);
	mpz_clear(a2_int);
	mpz_clear(best_a);
	mpfr_clear(best_ratio);
	mpfr_clear(ratio);
	mpz_clear(v);
	return;
}

mpz_t*calculate_array_Bk_f(int*number_prime_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1){
	mpz_t*array_Bk=NULL;
	long vpk;
	char t;

	mpz_t ak,root,ak_inverse;
	mpz_t pk,n_temp;

	mpz_init(ak);
	mpz_init(root);
	mpz_init(ak_inverse);
	mpz_init(n_temp);
	mpz_init(pk);
	if(n==NULL || a==NULL || b1==NULL || card_factor_base<=0){
		handle_error_with_exit("error in calculate_array_Bk\n");
	}
	if(mpz_get_si(a)==0){//se a=0
		mpz_clear(ak);
		mpz_clear(root); 
		mpz_clear(ak_inverse);
		mpz_clear(n_temp);
		mpz_clear(pk);
		return NULL;//ritorna NULL
	}
	if(mpz_sgn(n)<=0 || s<0 || mpz_sgn(a)<=0){
		handle_error_with_exit("error in calculate_array_Bk\n");
	}
	long index=0;
	mpz_set_si(b1,0);//b1==0
	array_Bk=alloc_array_mpz(s);
	for(int k=0;k<s;k++){
		vpk=number_prime_a[index];
		index++;
		mpz_set_si(pk,vpk);//pk=vpk,pk=primo iesimo scelto per a
		if(mpz_divisible_p(a,pk)==0){
			handle_error_with_exit("error in mpz_divisible_calculate array_BK\n");
		}
		mpz_divexact(ak,a,pk);//ak=a/pk
		if(mpz_invert(ak_inverse,ak,pk)==0){
		    handle_error_with_exit("error in mpz_invert calculate array_Bk\n");
		}//ak_inverse=(ak)^-1 mod pk
		mpz_set(n_temp,n);//n_temp=n
		mpz_mod(n_temp,n,pk);//n_temp=n mod pk
		if(vpk==2){//la radice di n mod 2 è la riduzione di n mod 2,in realtà i numeri primi devono essere dispari
			mpz_set(root,n_temp);//root =n_temp
		}
		else{
			t=quadratic_residue(root,n_temp,pk);//root=radice quadrata di n modulo pk
			if(t==-1 || t==0){
				handle_error_with_exit("error in calculate quadratic_residue\n");
			}
		}
		mpz_mul(root,root,ak);//root=root*ak
		mpz_mul(root,root,ak_inverse);//tk*ak*(ak)^-1=root*ak*(ak)^-1
		mpz_mod(root,root,a);//Bk può essere ridotto modulo a
		mpz_set(array_Bk[k],root);//array_Bk[k]=root*ak*ak_inverse;
		mpz_add(b1,b1,array_Bk[k]);//b1=(b1)+array_Bk[k],alla fine del for b1=somma di tutti i Bk
	}

	mpz_clear(pk);
	mpz_clear(n_temp);
	mpz_clear(ak);
	mpz_clear(root);
	mpz_clear(ak_inverse);
	return array_Bk;
}
void multiply_n_for_k(mpz_t n,int *k,char*factorized){//n=n*k in modo tale che n è un quadrato modulo 8,cioè n==1 mod 8
	if(n==NULL || mpz_sgn(n)<=0 || n==NULL || k==NULL || factorized==NULL){
		handle_error_with_exit("error in multiply_n_for_k\n");
	}
	if(mpz_divisible_2exp_p(n,1)){
		mpz_divexact_ui(n,n,2);
		gmp_printf("factorization of n=2*%Zd\n",n);
		/*fprintf(file_log,"p=2 ");
		fprintf(file_log,"q=");
		mpz_out_str(file_log,10,n);
		fprintf(file_log," ");*/
		*factorized=1;
		return;
	}
	return;
	/*mpz_t v;//valore temporaneo
	mpz_init_set_si(v,8);
	mpz_mod(v,n,v);
	*k=mpz_get_si(v);//k congruo a n mod 8,se n congruo a 1 mod 8 allora k=1
	mpz_mul(n,n,v);//n=n*k
	mpz_clear(v);
	return;*/
}
int calculate_v(int i){
	int v=1;
	while(i%2==0){//finquando è divisibile per 2
		i=i/2;//dividilo per 2
		v++;
	}
	return v;
}

mpz_t *calculate_bi(mpz_t *array_Bk,const mpz_t b1,int s){
	if(array_Bk==NULL){
		return NULL;
	}
	if(mpz_sgn(b1)<=0 || s<=0){//se array_Bk diverso da NULL s deve essere almeno 1
		handle_error_with_exit("error in calculate_bi\n");
	}
	long length=pow(2,s-1);//length=2^(s-1)
	mpz_t*array_bi=alloc_array_mpz(length);
	int v;
	int exp;
	mpz_t value;
	mpz_init(value);
	mpz_set(array_bi[0],b1);//prima posizione dell'array è b1
	for(int i=1;i<length;i++){
		v=calculate_v(i);
		exp=ceil((double)i/pow(2,v));// parte intera superiore i/2^v
		mpz_set_si(value,2);//value=2
		if((exp & 1)==0){//exp è pari
			//2*(-1)^[i/2^v]==2
		}
		else{
			mpz_set_si(value,-2);//2*(-1)^[i/2^v]==-2
		}
		mpz_mul(value,value,array_Bk[v-1]);//2*(-1)^[i/2^v]*Bv
		//value=2*(pow(-1,exp))*array_Bk[v-1];//2*(-1)^[i/2^v]*Bv
		mpz_add(array_bi[i],array_bi[i-1],value);//bi=bi-1+2*(-1)^[i/2^v]*Bv
	}
	mpz_clear(value);
	return array_bi;
}

void calculate_a_and_b(int*solution,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base,mpz_t a,mpz_t b,const mpz_t n){
	if(solution==NULL || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0 || card_f_base<=0 || num_B_smooth<card_f_base || a==NULL || b==NULL || n==NULL){
		handle_error_with_exit("error in parameter calculate a and b\n");
	}
	//il vettore solution ci dice quali relazioni vanno moltiplicate,ogni elemento di solution che è 1 è una relazione da moltiplicare con le altre
	mpz_t temp,v_temp;
	mpz_t num_square;
	
	mpz_init(num_square);
	mpz_init(temp);
	mpz_init(v_temp);
	
	mpz_set_si(a,1);//a=1 temporaneo
	mpz_set_si(b,1);//b=1 temporaneo
	mpz_t*v=create_array_temp_factorization(card_f_base,matrix_B_smooth[0]);//-1 0 2 0 3 0 ecc primi della factor base con esponente zero
	for(int i=0;i<num_B_smooth;i++){
		if(solution[i]==1){//se solution==1 allora bisogna moltiplicare le relazioni
			mpz_set(num_square,matrix_B_smooth[i][card_f_base*2]);//moltiplica le radici quadrate dei numeri B-smooth mod n
			mpz_mul(a,a,num_square);//a=a*square
			mpz_mod(a,a,n);//a=a*square mod n
			sum_elem_multiple_of_2_mpz(v,matrix_B_smooth[i],card_f_base*2,card_f_base*2);//moltiplica le relazioni sommando opportunamente gli esponenti tra i 2 numeri
		}
	}
	divide_elem_multiple_of_2_by_x(v,card_f_base*2,2);//divide gli esponenti per 2 per calcolare successivamente b
	for(int j=2;j<card_f_base*2;j=j+2){//per ogni primo della factor base...
		mpz_set(v_temp,v[j]);//v_temp=v[j]
		mpz_powm(v_temp,v_temp,v[j+1],n);//v_temp=pow(v[j],v[j+1]) mod n
		mpz_mul(b,b,v_temp);//b=primo factor base*esponente del primo della factor base
		mpz_mod(b,b,n);//b mod n
	}

	mpz_clear(temp);
	mpz_clear(num_square);
	mpz_clear(v_temp);
	free_memory_array_mpz(v,card_f_base*2);
	return;
}
void calculate_a_and_b_siqs(const int*solution,struct node_square_relation*head,int num_B_smooth,int card_f_base,mpz_t a,mpz_t b,const mpz_t n){
    if(solution==NULL || head==NULL || num_B_smooth<=0 || card_f_base<=0 || num_B_smooth<card_f_base*ENOUGH_RELATION || a==NULL || b==NULL || n==NULL){
        handle_error_with_exit("error in parameter calculate a and b\n");
    }
    //il vettore solution ci dice quali relazioni vanno moltiplicate,ogni elemento di solution che è 1 è una relazione da moltiplicare con le altre
    mpz_t v_temp;
    mpz_t square;

    mpz_init(square);
    mpz_init(v_temp);

    int exponent,index;
    mpz_set_si(a,1);//a=1 temporaneo
    mpz_set_si(b,1);//b=1 temporaneo
    struct node_square_relation*p=head;
    int*sum_exponent_relation=alloc_array_int(card_f_base);//contiene la somma degli esponenti dei numeri usati per le relazioni che forniscono la soluzione del sistema lineare
    //mpz_t*v=create_array_temp_factorization(card_f_base,matrix_B_smooth[0]);//-1 0 2 0 3 0 ecc primi della factor base con esponente zero
    for(int i=0;i<num_B_smooth;i++){//scansiona tutto l'array solution
        if(solution[i]==1){//se solution==1 allora bisogna moltiplicare le relazioni
            mpz_set(square,p->square_relation.square);//moltiplica le radici quadrate dei numeri B-smooth mod n
            mpz_mul(a,a,square);//a=a*square
            mpz_mod(a,a,n);//a=a*square mod n
            struct node_factorization*q=p->square_relation.head_factorization;
            while(q!=NULL){//scansiona fattorizzazione del numero
                exponent=q->exp_of_number;//ottieni esponente
                index=q->index;//ottieni indice
                q=q->next;
                sum_exponent_relation[index]+=exponent;//metti nell'array in posizione indice l'esponente+la somma preceente degli esponenti
            }
        }
        p=p->next;//passa alla prossima relazione,se il prossimo elemento di solution[i] è 0 passa al prossimo
    }
    divide_vector_multiple_of_2_by_2(sum_exponent_relation,card_f_base);
    for(int j=0;j<card_f_base;j++){//per ogni primo della factor base
        if(sum_exponent_relation[j]==0){
            continue;
        }
        mpz_set_si(v_temp,r.prime[j]);//v_temp=primo
        mpz_powm_ui(v_temp,v_temp,sum_exponent_relation[j],n);//v_temp=primo^exp mod n
        mpz_mul(b,b,v_temp);//b=primo factor base*esponente del primo della factor base
        mpz_mod(b,b,n);//b mod n
    }
    mpz_clear(square);
    mpz_clear(v_temp);
	free(sum_exponent_relation);
    return;
}
int try_to_factor(const mpz_t a,const mpz_t b,const mpz_t n,mpz_t factor1,mpz_t factor2){//dati a,b ed n in input prova a fattorizzare n con il crivello quadratico
	//ritorna il numero di fattorizzazioni trovate usando a e b(0,1,2)
	//se c'è almeno una fattorizzazione non banale imposta factor1 e factor2 come fattori
	if(a==NULL || b==NULL || n==NULL || mpz_sgn(a)<0 || mpz_sgn(b)<0 || mpz_sgn(n)<=0 || factor1==NULL || factor2==NULL){
		handle_error_with_exit("error in parameter try to factor\n");
	}
	mpz_t n_temp;
	char not_trivial_factor=0;
	mpz_t sum;//a+b;
	mpz_t diff;//a-b
	mpz_t gcd_diff;
	mpz_t gcd_sum;

	mpz_init(sum);
	mpz_init(diff);
	mpz_init(gcd_diff);
	mpz_init(gcd_sum);
	mpz_init(n_temp);

	mpz_add(sum,a,b);//sum=a+b;
	mpz_sub(diff,a,b);//diff=a-b;
	mpz_set_si(gcd_sum,1);
	mpz_set_si(gcd_diff,1);
	mpz_mod(sum,sum,n);
	mpz_mod(diff,diff,n);
	if(mpz_cmp(sum,n)!=0){//sum!=n
		mpz_gcd(gcd_sum,sum,n);
	}
	if(mpz_cmp_si(diff,0)!=0){//diff!=0
		mpz_gcd(gcd_diff,diff,n);
	}
	if(mpz_cmp_si(gcd_diff,1)!=0 && mpz_cmp(gcd_diff,n)!=0){//gcd_diff !=1 && gcd_diff!=n
		mpz_set(n_temp,n);//n_temp==n
		not_trivial_factor++;
		if(mpz_divisible_p(n_temp,gcd_diff)==0){
			handle_error_with_exit("error in mpz_divisible try to factor\n");
		}
		mpz_divexact(n_temp,n_temp,gcd_diff);//n=n/gcd_diff
		gmp_printf("factorization of n=%Zd*%Zd\n",gcd_diff,n_temp);
		mpz_set(factor1,gcd_diff);
		mpz_set(factor2,n_temp);
		mpz_mul(n_temp,factor1,factor2);
		if(mpz_cmp(n_temp,n)!=0){
			handle_error_with_exit("invalid factorization of n gcd diff\n");
		}
	}
	if(mpz_cmp_si(gcd_sum,1)!=0 && mpz_cmp(gcd_sum,n)!=0){//gcd_sum !=1 && gcd_sum!=n
		mpz_set(n_temp,n);//n_temp==n
		not_trivial_factor++;
		if(mpz_divisible_p(n_temp,gcd_sum)==0){
			handle_error_with_exit("error in mpz_divisible try to factor\n");
		}
		mpz_divexact(n_temp,n_temp,gcd_sum);//n=n/gcd_sum
		gmp_printf("factorization of n=%Zd*%Zd\n",gcd_sum,n_temp);
		mpz_set(factor1,gcd_sum);
		mpz_set(factor2,n_temp);
		mpz_mul(n_temp,factor1,factor2);
		if(mpz_cmp(n_temp,n)!=0){
			handle_error_with_exit("invalid factorization of n gcd sum\n");
		}
	}
	mpz_clear(sum);
	mpz_clear(diff);
	mpz_clear(gcd_diff);
	mpz_clear(gcd_sum);
	mpz_clear(n_temp);
	return not_trivial_factor;
}

int find_factor_of_n_from_base_matrix_char(int **base_matrix,int num_row,int* num_column,char*matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,struct node_square_relation*head,int num_B_smooth,int card_f_base){
//calcola tutte le soluzioni una per una e per ognuna prova a vedere se trova una coppia a,b che fattorizza n se riesce allora ritorna al chiamante
    if(base_matrix==NULL || *base_matrix==NULL || num_row<=0 || *num_column<=0 || matrix_linear_system==NULL
       || num_row_matrix<=0 || num_col_matrix<=0 || head==NULL || card_f_base<=0){
        handle_error_with_exit("error in parameter find_factor_of_n_from_base_matrix\n");
    }
    mpz_t n_copy,a,b,factor1,factor2;
    int*solution,num_col=*num_column;
    char*combination;
    int j;
    int factor_founded_from_a_and_b=0;
    int count_combination=0;
    int*col_h=NULL,*vector_sum=NULL;
    mpz_init(n_copy);
    mpz_init(a);
    mpz_init(b);
    mpz_init(factor1);
    mpz_init(factor2);
    mpz_set(n_copy,n);
    if(num_col>MAX_DIM_SOL){//troppe soluzioni portano ad una allocazione troppo grande della memoria considerare solo i primi 10 vettori della base
        //ciò equivale ad allocare 2^MAX_DIM_SOL soluzioni invece di 2^num_col,ciò in termini di soluzione corrisponde a considerare tutti i vettori
        //non allocati ==0,quindi trova soluzioni riga del tipo=v1,v2,...v10,0,0,0,0,0...0
        num_col=MAX_DIM_SOL;
        *num_column=num_col;
    }
    solution=alloc_array_int(num_row);//alloca soluzione
    combination=alloc_array_char(num_col);//alloca array per mantenere traccia delle combinazioni ottenute,è un array temporaneo
    count_combination=0;//cont le combinazioni scansionate finora
    while(!array_is_fill_of_value(combination,num_col,1)){//while combination !=11111111111
        //somma binaria nell'array temporaneo
        j=num_col-1;//indice dell'ultimo elemento di combination,si parte dalla fine dell'array e si somma sempre 1 in modo binario finquando non si arriva a 1111111111111...1111
        combination[j]=combination[j]+1;
        while(combination[j]==2){//trasporta il riporto dagli ultimi elementi ai primi
            combination[j]=0;
            combination[j-1]=combination[j-1]+1;
            j--;//passa all'elemento più significativo dell'array
        }

        //creazione della combinazione
        for(int h=0;h<num_col;h++){//scansiona l'array combination dall'inizio alla fine per generare la soluzione
            if(combination[h]==0){
                continue;
            }
            col_h=get_coli(base_matrix,num_row,num_col,h);
            vector_sum=sum_vector(solution,col_h,num_row,num_row);//somma tutti gli array che hanno 					indice uguale a 1 nell'array combination
            memcpy(solution,vector_sum,sizeof(int)*num_row);
            free(vector_sum);
            free(col_h);
        }
        count_combination++;
        reduce_array_mod_n(solution,num_row,2);//riduce la matrice soluzione mod 2
        /*if(check_if_array_is_reduce_mod_n(solution,num_row,2)==0){
            handle_error_with_exit("error in calculate solution\n");
        }
        if(verify_solution_char(matrix_linear_system,num_row_matrix,num_col_matrix,solution)==0){
            handle_error_with_exit("invalid solution in find_factor_of_n_from_base_matrix\n");
        }*/
        calculate_a_and_b_siqs(solution,head,num_B_smooth,card_f_base,a,b,n_copy);//calcola un a e un b
        factor_founded_from_a_and_b=try_to_factor(a,b,n_copy,factor1,factor2);
        if(factor_founded_from_a_and_b>0){
            free(combination);
            free(solution);
            mpz_clear(n_copy);
            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(factor1);
            mpz_clear(factor2);
            return factor_founded_from_a_and_b;
        }
    }
    free(combination);
    free(solution);
    mpz_clear(n_copy);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(factor1);
    mpz_clear(factor2);
    return 0;

}
