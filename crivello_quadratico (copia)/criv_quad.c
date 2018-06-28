#include "criv_quad.h"
#include "math_function.h"
#include <gmp.h>
#include <mpfr.h>
#include "list_factor_base.h"
#include "list_square_relation.h"
#include "list_factorization.h"

extern long k;
extern struct timespec timer;
extern struct timespec time_start;
//extern FILE*file_log;
extern int num_increment_M_and_B;
extern struct row_factorization r;
void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B){
	//calcola valori di M e B cercando di minimizzare il numero di volte che i numeri devono essere aumentati
	if(n==NULL || digit_n<=0 || M==NULL || B==NULL){
		handle_error_with_exit("error in calculate best M and B\n");
	}
	//*B=100*1000*1000;
	//*M=5000000
	//*M=400000;
	//*B=4000000;

	if(digit_n<7){
		*M=25;
		*B=20;
		return;
	}
	if(digit_n<9){
		*M=68;
		*B=38;
		return;
	}
	if(digit_n<11){
		*M=125;
		*B=59;
		return;
	}
	if(digit_n<13){
		*M=303;
		*B=112;
		return;
	}
	if(digit_n<15){
		*M=492;
		*B=156;
		return;
	}
	if(digit_n<17){
		*M=1177;
		*B=277;
		return;
	}
	if(digit_n<19){
		*M=2835;
		*B=502;
		return;
	}
	if(digit_n<21){
		*M=9922;
		*B=1668;
		return;
	}
	if(digit_n<23){
		*M=14958;
		*B=2577;
		return;
	}
	if(digit_n<25){
		*M=22512;
		*B=3540;
		return;
	}
	if(digit_n<40){//da 37 cifre a 39,4 secondi di media
		*M=22512;
		*B=3540;
		return;
	}
	if(digit_n<45){//da 42 a 44 cifre 4 secondi di media
		handle_error_with_exit("error criv_quad <45\n");
		*M=20000;
		*B=8000;
		return;
	}
	if(digit_n<50){//da 47 a 49 cifre, 6,5 secondi di media
		handle_error_with_exit("error criv_quad <50\n");
		*M=20000;
		*B=19000;
		return;
	}
	if(digit_n<55){//da 53 a 54 cifre,
		*M=20000;
		*B=19000;
	}
	if(digit_n<60){
		handle_error_with_exit("error criv_quad <60\n");
		*M=35000;
		*B=15000;
	}
	if(digit_n<65){
		handle_error_with_exit("error criv_quad <65\n");
		*M=35000;
		*B=15000;
		return;
	}
	if(digit_n<70){
		handle_error_with_exit("error criv_quad <70\n");
		*M=35000;
		*B=15000;
		return;
	}
	if(digit_n<75){
		handle_error_with_exit("error criv_quad <75\n");
		*M=35000;
		*B=15000;
		return;
	}
	if(digit_n<80){
		handle_error_with_exit("error criv_quad <80\n");
		*M=35000;
		*B=15000;
		return;
	}
	if(digit_n>=80){
		handle_error_with_exit("error criv_quad >=80\n");
		*M=35000;
		*B=15000;
		return;
	}
	handle_error_with_exit("num_too long\n");
	return;
}
void calculate_news_M_and_B(long*M,long*B){
	//calcola i valori nuovi di M e B secondo una formula scelta opportunamente
	//new_m=(m+perc_m)+(m+perc_m)*perc_m/100
	//new_b=(b+perc_b)+(b+perc_b)*perc_b/100
	return;
	if(M==NULL || B==NULL || *M<=0 || *B<=0){
		handle_error_with_exit("error in calculate_news_M_and_B\n");
	}
	double temp;
	//calculate new M
	temp=*M+PERC_INCREMENT_M;
	temp=temp+(temp*PERC_INCREMENT_M)/100;
	*M=temp;
	
	//calculate new B
	temp=*B+PERC_INCREMENT_B;
	temp=temp+(temp*PERC_INCREMENT_B)/100;
	*B=temp;
	num_increment_M_and_B++;
	return;
}
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

/*void calculate_b_mpqs(mpz_t b_mpqs,const mpz_t q,const mpz_t n,const mpz_t a_mpqs){
	//da thresold_q calcolare una radice di N mod thresold_q 
	// e con la radice trovata calcolare la radice di N mod q^2
	//thresold_q è primo
	if(b_mpqs==NULL || q==NULL || n==NULL || a_mpqs==NULL){
		handle_error_with_exit("error in calculate b_mpqs\n");
	}
	long p;
	char t;
	mpz_t root;
	mpz_init(root);
	t=quadratic_residue(root,n,q);//root^2 congruo a n mod q
	if(mpz_cmp_si(a_mpqs,0)==0){
		mpz_set_si(b_mpqs,0);
		return;
	}
	if(t==-1 || t==0){
		handle_error_with_exit("error in calculate b_mpqs\n");
	}
	p=mpz_get_si(q);//p=q
	square_root_mod_p_to_k(b_mpqs,root,p,n,2);//b_mpqs=radice quadrata di n modulo q^2
	mpz_clear(root);
	return;
}*/

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



/*void calculate_best_M(const mpz_t n,long*M){
	if(M==NULL || n==NULL){
		handle_error_with_exit("error in calculate best M\n");
	}
	mpz_t temp,temp1,ten;
	mpz_init(temp);
	mpz_init(temp1);
	mpz_init(ten);

	mpz_set_si(ten,10);//ten=10
	mpz_pow_ui(temp,ten,0);//temp=10^0
	mpz_pow_ui(temp1,ten,10);//temp=10^10
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^0<=n<10^10
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,10);//temp=10^10
	mpz_pow_ui(temp1,ten,20);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^10<=n<10^20
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,20);//temp=10^10
	mpz_pow_ui(temp1,ten,30);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^20<=n<10^30
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,30);//temp=10^10
	mpz_pow_ui(temp1,ten,40);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^30<=n<10^40
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,40);//temp=10^10
	mpz_pow_ui(temp1,ten,50);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^40<=n<10^50
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,50);//temp=10^10
	mpz_pow_ui(temp1,ten,60);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^50<=n<10^60
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,60);//temp=10^10
	mpz_pow_ui(temp1,ten,70);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^60<=n<10^70
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,70);//temp=10^10
	mpz_pow_ui(temp1,ten,80);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^70<=n<10^80
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,80);//temp=10^10
	mpz_pow_ui(temp1,ten,90);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^80<=n<10^90
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,90);//temp=10^10
	mpz_pow_ui(temp1,ten,100);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^90<=n<10^100
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	
	mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
	return;
}
void calculate_best_B(const mpz_t n,long*B){
	if(B==NULL || n==NULL){
		handle_error_with_exit("error in calculate best M\n");
	}
	mpz_t temp,temp1,ten;
	mpz_init(temp);
	mpz_init(temp1);
	mpz_init(ten);

	mpz_set_si(ten,10);//ten=10
	mpz_pow_ui(temp,ten,0);//temp=10^0
	mpz_pow_ui(temp1,ten,10);//temp=10^10
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^0<=n<10^10
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,10);//temp=10^10
	mpz_pow_ui(temp1,ten,20);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^10<=n<10^20
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,20);//temp=10^10
	mpz_pow_ui(temp1,ten,30);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^20<=n<10^30
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,30);//temp=10^10
	mpz_pow_ui(temp1,ten,40);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^30<=n<10^40
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,40);//temp=10^10
	mpz_pow_ui(temp1,ten,50);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^40<=n<10^50
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,50);//temp=10^10
	mpz_pow_ui(temp1,ten,60);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^50<=n<10^60
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,60);//temp=10^10
	mpz_pow_ui(temp1,ten,70);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^60<=n<10^70
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,70);//temp=10^10
	mpz_pow_ui(temp1,ten,80);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^70<=n<10^80
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,80);//temp=10^10
	mpz_pow_ui(temp1,ten,90);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^80<=n<10^90
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	mpz_pow_ui(temp,ten,90);//temp=10^10
	mpz_pow_ui(temp1,ten,100);//temp=10^20
	if(mpz_cmp(n,temp)>=0 && mpz_cmp(n,temp1)<0){ //10^90<=n<10^100
		mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
		return;
	}
	
	mpz_clear(temp);
	mpz_clear(temp1);
	mpz_clear(ten);
	return;
}*/
/*struct node_factor_base* create_factor_base_f(int*cardinality_factor_base,long B,struct node_factor_base**tail,const mpz_t n,int *start_prime_factor_base){//crea la factor base aggiungendo i numeri primi minori o uguali a B
	if(B<2 || tail==NULL || cardinality_factor_base==NULL || n==NULL){
		handle_error_with_exit("error in parameter\n");
	}
	struct node_factor_base *head=NULL;
	struct node_factor_base*node;
	long v,m;
	mpz_t p,pp,value;
	

	mpz_init(p);
	mpz_init(pp);
	mpz_init(value);

	mpz_set_si(p,2);//p=2
	*cardinality_factor_base=1;
	insert_ordered_f(2,n,&head,tail);//inserisce 2 nella factor base
	for(long i=2;i<=B;){
		mpz_set_si(p,i);//p=i,inizialmente p=2,ad ogni inizio ciclo deve essere pari
		mpz_nextprime(p,p);//p=next_prime,ritorna il prossimo primo maggiore strettamente di p,è bene che p sia pari
		if(mpz_cmp_si(p,B)<=0){//primo<=B
			v=mpz_get_si(p);
			v=(v-1)>>1;//shift a destra divide per 2
			mpz_powm_ui(value,n,(unsigned long int)v,p);//value=n^v mod p
			m=mpz_get_si(value);//m=n^((p-1)/2) mod p
			if(m==1){//n è quadrato modulo p
				node=get_new_node_f(mpz_get_si(p),n);
				insert_at_tail_f(node,&head,tail);
				(*cardinality_factor_base)++;
			}
			i=mpz_get_si(p)+1;//il numero diventa pari e il prossimo numero primo sarà dispari
		}
		else{//B è stato superato
			break;
		}
	}
	mpz_set_si(p,-1);//temp=-1
	insert_ordered_f(-1,n,&head,tail);//inserisci -1
	(*cardinality_factor_base)++;
	mpz_clear(p);
	mpz_clear(pp);
	mpz_clear(value);
	return head;
}*/
struct node_factor_base*initialize_factor_base(int*cardinality_factor_base,long B,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base){
	if(B<2 || tail==NULL || cardinality_factor_base==NULL || n==NULL || last_prime_factor_base==NULL){
		handle_error_with_exit("error in parameter\n");
	}
	struct node_factor_base *head=NULL;
	mpz_t p;
	mpz_init(p);
	mpz_set_si(p,-1);//temp=-1
	insert_ordered_f(-1,n,&head,tail);//inserisci -1 nella factor base
	(*cardinality_factor_base)++;
	mpz_set_si(p,2);//p=2
	(*cardinality_factor_base)++;
	insert_ordered_f(2,n,&head,tail);//inserisce 2 nella factor base
	*last_prime_factor_base=2;
	mpz_clear(p);
	return head;
}

void create_factor_base_f(int*cardinality_factor_base,long B,struct node_factor_base**head,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base){//crea la factor base aggiungendo i numeri primi minori o uguali a B
    if(B<2 || tail==NULL || head==NULL || cardinality_factor_base==NULL || n==NULL || last_prime_factor_base==NULL ){
        handle_error_with_exit("error in parameter\n");
    }
    struct node_factor_base*node;
    long v,v_temp,m;
    mpz_t p,pp,value;


    mpz_init(p);
    mpz_init(pp);
    mpz_init(value);
	long start;
	if(*last_prime_factor_base==1) {
		*head=initialize_factor_base(cardinality_factor_base, B,tail, n,last_prime_factor_base);
		//initialize factor base pone last_prime_factor_base a 1
	}
	if((*last_prime_factor_base & 1)==0){//start=last_prime_factor_base è pari==2
		start=*last_prime_factor_base;
	}
	else{//last_prime_factor_base è dispari
		start=*last_prime_factor_base+1;//start è pari
	}
    for(long i=start;i<=B;){
        mpz_set_si(p,i);//p=i,inizialmente p=2,ad ogni inizio ciclo deve essere pari
        mpz_nextprime(p,p);//p=next_prime,ritorna il prossimo primo maggiore strettamente di p,è bene che p sia pari
        if(mpz_cmp_si(p,B)<=0){//primo<=B
            v=mpz_get_si(p);
            v_temp=v;
            v=(v-1)>>1;//shift a destra divide per 2,v=(v-1/2)
            mpz_powm_ui(value,n,(unsigned long int)v,p);//value=n^v mod p
            m=mpz_get_si(value);//m=n^((p-1)/2) mod p
            if(m==1){//n è quadrato modulo p
                node=get_new_node_f(mpz_get_si(p),n);
                insert_at_tail_f(node,head,tail);
                (*cardinality_factor_base)++;
                *last_prime_factor_base=v_temp;
            }
            i=mpz_get_si(p)+1;//il numero diventa pari e il prossimo numero primo sarà dispari
        }
        else{//B è stato superato
            break;
        }
    }
    mpz_clear(p);
    mpz_clear(pp);
    mpz_clear(value);
    return;
}
/*struct node* create_factor_base(int*cardinality_factor_base,long B,struct node**tail,const mpz_t n,mpz_t q,const mpz_t thresold_q){//crea la factor base aggiungendo i numeri primi minori o uguali a B
	if(B<2 || tail==NULL || cardinality_factor_base==NULL || n==NULL || q==NULL || thresold_q==NULL){
		handle_error_with_exit("error in parameter\n");
	}
	struct node *head=NULL;
	struct node*node;
	char found=0;
	long v,m;
	mpz_t p,pp,value;
	mpz_t delta_min,delta_max,min,max;
	

	mpz_init(p);
	mpz_init(pp);
	mpz_init(value);
	mpz_init(delta_min);
	mpz_init(delta_max);
	mpz_init(min);
	mpz_init(max);

	mpz_set_si(p,2);//p=2
	*cardinality_factor_base=1;
	insert_ordered(p,&head,tail);//inserisce 2 nella factor base
	mpz_set_si(min,2);
	for(long i=3;i<=B;){
		mpz_set_si(p,i);//p=i
		if(mpz_divisible_ui_p(p,2)!=0){//divisibile per 2,p è pari
		}
		else{
			mpz_sub_ui(p,p,1);//p è dispari,p da dispari diventa pari
		}
		mpz_nextprime(p,p);//p=next_prime,ritorna il prossimo primo maggiore di p
		if(mpz_cmp_si(p,B)<=0){//primo<=B
			mpz_set(pp,p);//pp=p
			mpz_sub_ui(pp,pp,1);//pp=(p-1)
			if(mpz_divisible_ui_p(pp,2)==0){
				handle_error_with_exit("error in mpz_divisible,create factor base\n");
			}
			mpz_divexact_ui(pp,pp,2);//p=(p-1)/2
			v=mpz_get_si(pp);//v=(p-1)/2
			mpz_powm_ui(value,n,(unsigned long int)v,p);
			m=mpz_get_si(value);//m=n^((p-1)/2) mod p
			if(m==1){//n è quadrato modulo p
				node=get_new_node(p);
				insert_at_tail(node,&head,tail);
				(*cardinality_factor_base)++;
				if(mpz_cmp(p,thresold_q)<0){//p minore della soglia
					mpz_set(min,p);//min=p
					mpz_sub(delta_min,thresold_q,p);//delta_min=thresold_q-p
				}
				else if(found==0){//p maggiore della soglia
					mpz_set(max,p);//max=p
					mpz_sub(delta_max,p,thresold_q);//delta_max=p-thresold_q
					found=1;
				}
			}
			mpz_add_ui(p,p,1);
			i=mpz_get_si(p);//il numero diventa pari e il prossimo numero primo sarà dispari
		}
		else{//B è stato superato
			break;
		}
	}
	mpz_set_si(p,-1);//temp=-1
	insert_ordered(p,&head,tail);//inserisci -1
	(*cardinality_factor_base)++;
	if(mpz_cmp(delta_min,delta_max)<0 && mpz_cmp_si(min,2)!=0){//delta min è più piccolo di delta max
		mpz_set(q,min);
	}
	else{//se la distanza tra q e max è minore o se min era uguale a 2 scegli max
		mpz_set(q,max);
	}
	mpz_clear(p);
	mpz_clear(pp);
	mpz_clear(value);
	mpz_clear(delta_min);
	mpz_clear(delta_max);
	mpz_clear(min);
	mpz_clear(max);
	return head;
}*/

/*void calculate_p_min_p_max_i(long*p_min_i,long*p_max_i,struct node*head_f_base,long cardinality_factor_base){
	if(p_min_i==NULL || p_max_i==NULL || head_f_base==NULL || cardinality_factor_base<=0){
		handle_error_with_exit("error in calculate _p_min_p_max_i\n");
	}
	struct node*list=head_f_base;
	mpz_t temp;
	mpz_init(temp);
	long i=0;
	while(list!=NULL){
		mpz_set(temp,list->prime);//valore del primo della factor base
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
}*/
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
/*void calculate_target_a1(mpfr_t target_a1,const mpfr_t target_a,struct node*head_f_base,long p_min_i,long p_max_i,int cardinality_factor_base){
	if(target_a==NULL || head_f_base==NULL || p_min_i<0 || p_max_i<=0 || p_min_i>p_max_i || target_a1==NULL || 
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
	struct node*p=head_f_base;
	while(count<p_min_i){
		p=p->next;
		count++;
	}
	if(count!=p_min_i){
		handle_error_with_exit("error in calculate target_a1\n");
	}
	mpz_set(prime_min,p->prime);
	while(count<p_max_i){
		p=p->next;
		count++;
	}
	if(count!=p_max_i){
		handle_error_with_exit("error in calculate target_a1\n");
	}
	mpz_set(prime_max,p->prime);
	mpz_add(prime_sum,prime_max,prime_min);//prime_sum=prime_min+prime_max
	mpfr_set_z(prime_avg,prime_sum,MPFR_RNDN);//prime_avg=prime_sum
	mpfr_div_2ui(prime_avg,prime_avg,1,MPFR_RNDN);//prime_avg=(pmax+pmin)/2
	mpfr_sqrt(sqrt_prime_avg,prime_avg,MPFR_RNDN);//rad((pmax+pmin)/2)
	mpfr_set(temp,target_a,MPFR_RNDN);//temp=target
	mpfr_div(temp,temp,sqrt_prime_avg,MPFR_RNDN);//temp=target/rad((pmax+pmin)/2)
	mpfr_set(target_a1,temp,MPFR_RNDN);//target_a1=temp=target/rad((pmax+pmin)/2)

	mpz_clear(prime_min);
	mpz_clear(prime_max);
	mpz_clear(prime_sum);
	mpfr_clear(prime_avg);
	mpfr_clear(sqrt_prime_avg);
	mpfr_clear(temp);
	print_time_elapsed("time to calculate target_a1");
	return;
}*/
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
	printf("prime_avg=");
	mpfr_out_str(stdout,10,0,prime_avg,MPFR_RNDN);
	printf("\n");
	mpfr_sqrt(sqrt_prime_avg,prime_avg,MPFR_RNDN);//sqrt_prime=rad((pmax+pmin)/2)
	mpfr_set(temp,target_a,MPFR_RNDN);//temp=target=rad(2*n)/M
	mpfr_div(temp,temp,sqrt_prime_avg,MPFR_RNDN);//temp=target/rad((pmax+pmin)/2)
	mpfr_set(target_a1,temp,MPFR_RNDN);//target_a1=temp=target/rad((pmax+pmin)/2)
	printf("target_a1=");
	mpfr_out_str(stdout,10,0,target_a1,MPFR_RNDN);
	printf("\n");
	mpz_clear(prime_min);
	mpz_clear(prime_max);
	mpz_clear(prime_sum);
	mpfr_clear(prime_avg);
	mpfr_clear(sqrt_prime_avg);
	mpfr_clear(temp);
	return;
}


/*void calculate_a(mpz_t a,const mpfr_t target_a,int*s,struct node*head_f_base,long cardinality_factor_base,long**array_of_prime_chosen_for_a){
	if(s==NULL || target_a==NULL || mpfr_sgn(target_a)<0 || head_f_base==NULL || cardinality_factor_base<=0 || *array_of_prime_chosen_for_a!=NULL || a==NULL){
		handle_error_with_exit("error in calculate_a\n");
	}
	long p_min_i=0;//indice del massimo primo da scegliere rispetto alla factor base
	long p_max_i=0;//indice del minimo primo da scegliere rispetto alla factor base
	//p_max_i>p_min_i

	long p_i;//indici dei primi scelti nella factor base per rappresentare a
	long s_max;//valore massimo di s=pmax-pmin+1
	mpz_t v;
	int iter=0,iter2=0;
	mpz_init(v);
	if(cardinality_factor_base<3){
		mpz_set_si(a,0);//poni a=0
		*s=0;
		mpz_clear(v);
		*array_of_prime_chosen_for_a=NULL;
		fprintf(file_log,"pmin=0,p_max=0 ");
		return;//ritorna array_of_prime_chosen_for_a==NULL e a=0
	}
	if(cardinality_factor_base==3){
		get_element_linked_list(a,head_f_base,2);//a=valore del 3 primo della factor base
		*s=1;
		*array_of_prime_chosen_for_a=alloc_array_long(cardinality_factor_base*2);//deve essere riempito con 1 nelle posizioni in cui c'è il 			primo scelto
		for(long j=0;j<cardinality_factor_base;j++){
			get_element_linked_list(v,head_f_base,j);//riempie gli elementi di posto dispari(indice pari) delle righe con i primi della 				factor base 
			(*array_of_prime_chosen_for_a)[2*j]=mpz_get_si(v);//scrive il valore del primo della factor base in posizione dispari
		}
		fprintf(file_log,"p_min=p_max=2 ");
		(*array_of_prime_chosen_for_a)[5]=1;
		mpz_clear(v);
		return;
	}
	mpz_t p_temp;
	mpfr_t target_a1;
	mpfr_t best_a;
	mpfr_t best_ratio;
	mpfr_t ratio;
	mpfr_t a2;
	mpfr_t p_rational;
	
	mpfr_init(p_rational);
	mpfr_init(a2);
	mpz_init(p_temp);
	mpfr_init(target_a1);

	mpfr_init(best_a);
	mpfr_init(best_ratio);//double
	mpfr_init(ratio);//double
	
	mpfr_set_si(best_a,0,MPFR_RNDN);//best_a=0
	mpfr_set_si(best_ratio,0,MPFR_RNDN);//best_target=0

	calculate_p_min_p_max_i(&p_min_i,&p_max_i,head_f_base,cardinality_factor_base);
	printf("p_min=%ld p_max=%ld\n",p_min_i,p_max_i);
	fprintf(file_log,"p_min=%ld p_max=%ld ",p_min_i,p_max_i);
	s_max=p_max_i-p_min_i+1;//nel caso peggiore s=p_max-p_min+1
	printf("s_max=%ld\n",s_max);
	calculate_target_a1(target_a1,target_a,head_f_base,p_min_i,p_max_i,cardinality_factor_base);
	if(mpfr_cmp_si(target_a1,1)<=0){
		mpz_set_si(a,0);//poni a=0
		*s=0;
		mpz_clear(v);
		*array_of_prime_chosen_for_a=NULL;
		mpfr_clear(p_rational);
		mpfr_clear(a2);
		mpz_clear(p_temp);
		mpfr_clear(target_a1);
		mpfr_clear(best_a);
		mpfr_clear(best_ratio);//double
		mpfr_clear(ratio);//double
		return;
	}
	printf("target_a1=");
	mpfr_out_str(stdout,10,0,target_a1,MPFR_RNDN);
	printf("\n");
	char f,b,c,d;
	int count=0;//conta quante volte while(mpz_cmp(a,target_a1)<0) è verificata
	long*q=alloc_array_long(s_max);//array che contiene gli indici dei primi scelti
	long*q_number=alloc_array_long(s_max);//array che contiene i dei primi scelti
	long*best_q_number=alloc_array_long(s_max);//array che contiene i dei primi scelti
	long*best_q=alloc_array_long(s_max);
	int length_best_q=0;
	for(int i=0;i<NUM_ITER_FOR_CALCULATE_A;i++){//iterazioni per cercare di migliorare a
		mpfr_set_si(a2,1,MPFR_RNDN);//a2=1 azzera a2
		memset(q,0,sizeof(long)*s_max);//azzera q
		memset(q_number,0,sizeof(long)*s_max);//azzera q
		count=0;//azzera le append a q
		iter=0;
		while(mpfr_cmp(a2,target_a1)<0 && iter<=MAX_ITER){//finquando a2<target_a1 oppure sono state raggiunte tot iterazioni
			p_i=2;
			iter2=0;
			while((p_i==2 || is_in_array_long(q,s_max,p_i)) && iter2<=MAX_ITER2){//scegli un p_i che non è stato ancora scelto
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
			get_element_linked_list(p_temp,head_f_base,p_i);//prendi l'elemento p_i-esimo dalla factor base
			mpfr_set_z(p_rational,p_temp,MPFR_RNDN);//p_rational=p_temp
			mpfr_mul(a2,a2,p_rational,MPFR_RNDN);//a2=a2*p
			if(count>=s_max){
				handle_error_with_exit("error in calculate a,index out of bounds\n");
			}
			q_number[count]=mpz_get_si(p_temp);
			q[count]=p_i;
			count++;//aumenta il conto delle append fatte
			iter++;
		}
		//una sequenza di numeri è stata scelta vediamo quanto è buona(vedendo i rapporti della migliore attualmente),se è buona la 			aggiorniamo
		mpfr_div(ratio,a2,target_a,MPFR_RNDN);//ratio=a2/target_a
		f=(mpfr_cmp_si(best_ratio,0)==0);
		b=mpfr_cmp_d(best_ratio,0.9)<0;
		c=mpfr_cmp(ratio,best_ratio)>0;
		d=( (mpfr_cmp_d(ratio,0.9)>=0) && (mpfr_cmp(ratio,best_ratio)<0) );
		if((f || d || b)  && c  ){
			mpfr_set(best_a,a2,MPFR_RNDN);//best_a=a2
			mpfr_set(best_ratio,ratio,MPFR_RNDN);//best_ratio=ratio
			memcpy(best_q,q,sizeof(long)*s_max);//best_q=q
			memcpy(best_q_number,q_number,sizeof(long)*s_max);//best_q=q
			length_best_q=count;//la lunghezza di best_q coincide con il numero di append fatte
		}
	}
	print_array_long(best_q,length_best_q);
	print_array_long(best_q_number,length_best_q);
	mpz_set_si(a,1);
	for(int i=0;i<length_best_q;i++){
		mpz_mul_si(a,a,best_q_number[i]);
	}
	*s=length_best_q;//imposta il valore di s
	printf("s=%d,s_max=%ld\n",*s,s_max);
	print_time_elapsed("time to choose factor for a");
	if(*s>s_max){
		handle_error_with_exit("error in calculate a,s must be minor or equal to s_max\n");
	}
	//creazione array dei primi scelti
	*array_of_prime_chosen_for_a=alloc_array_long(cardinality_factor_base*2);//deve essere riempito con 1 nelle posizioni in cui c'è il primo 		scelto
	struct node*p=head_f_base;
	for(long j=0;j<cardinality_factor_base;j++){
			//get_element_linked_list(v,head_f_base,j);//riempie gli elementi di posto dispari(indice pari) delle righe con i primi 				della factor base 
			(*array_of_prime_chosen_for_a)[2*j]=mpz_get_si(p->prime);//scrive il valore del primo della factor base in posizione 				dispari
			p=p->next;

	}
	adjust_array_of_prime(*array_of_prime_chosen_for_a,cardinality_factor_base*2,best_q,*s);//aggiusta array dei primi scelti rispetto a best_q
	free(q);
	free(best_q);
	mpfr_clear(p_rational);
	mpz_clear(p_temp);
	mpfr_clear(target_a1);
	mpfr_clear(a2);
	mpfr_clear(best_a);
	mpfr_clear(best_ratio);
	mpfr_clear(ratio);
	mpz_clear(v);
	return;
}*/

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
		//fprintf(file_log,"pmin=0,p_max=0 ");
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
	printf("p_min=%ld p_max=%ld\n",p_min_i,p_max_i);
	//fprintf(file_log,"p_min=%ld p_max=%ld ",p_min_i,p_max_i);
	s_max=p_max_i-p_min_i+1;//nel caso peggiore s=p_max-p_min+1
	if(s_max>S_MAX){
		s_max=S_MAX;
	}
	printf("s_max=%ld\n",s_max);
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
	printf("target_a1=");
	mpfr_out_str(stdout,10,0,target_a1,MPFR_RNDN);
	printf("\n");
    print_time_elapsed("time to calculate target_a1");

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
	printf("s=%d,s_max=%ld\n",*s,s_max);
	if(*s>s_max){
		handle_error_with_exit("error in calculate a,s must be minor or equal to s_max\n");
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
/*void add_remainder_to_matrix_factorization(mpz_t **matrix_factorization,mpz_t*array_number,long M,int cardinality_factor_base){
	//aggiunge il fattore rimanente dall'elemento dell'array di numeri una volta che viene diviso per tutti i numeri della factor base
	if(*matrix_factorization==NULL || matrix_factorization==NULL || array_number==NULL || M<=0 || cardinality_factor_base<=0){
		handle_error_with_exit("error in add_remainder_to_matrix_factorization\n");
	}
	for(long i=0;i<2*M+1;i++){
		mpz_set(matrix_factorization[i][cardinality_factor_base*2+1],array_number[i]);
	}
	return;
}*/

/*mpz_t*calculate_array_Bk(long*array_of_prime_chosen_for_a,int card_factor_base,const mpz_t n,long s,const mpz_t a,mpz_t b1){
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
	long length_array=2*card_factor_base;
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
		while(array_of_prime_chosen_for_a[2*index+1]!=1){//ricerca i valori diversi da 0 nell'array dei numeri scelti
			index++;
		}
		if(index>=length_array || array_of_prime_chosen_for_a[2*index+1]!=1){
			handle_error_with_exit("error in calculate array_Bk index too big\n");
		}
		vpk=array_of_prime_chosen_for_a[2*index];//pk è "piccolo",vpk è valore del primo k-esimo di a
		index++;
		mpz_set_si(pk,vpk);//pk=vpk
		if(mpz_divisible_p(a,pk)==0){
			handle_error_with_exit("error in mpz_divisible_calculate array_BK\n");
		}
		mpz_divexact(ak,a,pk);//ak=a/pk
		mpz_invert(ak_inverse,ak,pk);//ak_inverse=(ak)^-1 mod pk
		mpz_set(n_temp,n);//n_temp=n
		mpz_mod(n_temp,n,pk);//n_temp=n mod pk
		if(vpk==2){//la radice di n mod 2 è la riduzione di n mod 2
			mpz_set(root,n_temp);//root =n_temp
		}
		else{
			t=quadratic_residue(root,n_temp,pk);//root=radice quadrata di n modulo pk
			if(t==-1 || t==0){
				handle_error_with_exit("error in calculate quadratic_residue\n");
			}
		}
		mpz_mul(root,root,ak);//root=root*ak
		mpz_mul(root,root,ak_inverse);//tk*ak*(ak)^-1
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
}*/
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
/*mpz_t *create_array_of_number(const mpz_t a,const mpz_t b,long M,const mpz_t n){//crea array di numeri da setacciare per vedere quali sono B smooth
	if(a==NULL || b==NULL || n==NULL || mpz_sgn(n)<=0 || M<=0 || mpz_sgn(a)<=0 || mpz_sgn(b)<=0){
		handle_error_with_exit("invalid parameter create_array_of_number\n");
	}
	mpz_t *array=alloc_array_mpz(2*M+1);
	long j=-M;
	mpz_t c;
	mpz_t v,v1;

	mpz_init(c);
	mpz_init(v);
	mpz_init(v1);

	//calcolo di c
	mpz_mul(v,b,b);//v=b^2
	mpz_sub(v,v,n);//v=b^2-n
	if(mpz_divisible_p(v,a)==0){
		handle_error_with_exit("error in mpz_divisible create array_of_number\n");
	}
	mpz_divexact(c,v,a);//c=(b^2-n)/a

	for(int i=0;i<2*M+1;i++){
		//array[i]=(x0+j)*(x0+j)-n;//se a=1 b=x0
		//array[i]=a*j^2+2*j*b+c;//a(j)=a*(j^2)+2*bj+c se a=1 b=x0 a(j)=(x0+j)^2-n

		mpz_mul_si(v,a,j*j);//v=a*j^2
		mpz_set_si(v1,2*j);//v1=2*j
		mpz_mul(v1,v1,b);//v1=2*j*b
		mpz_add(v,v,v1);//v=a*j^2+2*j*b
		mpz_add(v,v,c);//v=a*j^2+2*j*b+c
		//array[i]=a*j^2+2*j*b+c;//a(j)=a*(j^2)+2*bj+c se a=1 b=x0 a(j)=(x0+j)^2-n
		mpz_set(array[i],v);//array[i]=a*j^2+2*j*b+c;
		j++;
	}

	mpz_clear(c);
	mpz_clear(v);
	mpz_clear(v1);
	return array;
}*/
/*mpz_t *create_array_of_number(const mpz_t a,const mpz_t b,long M,const mpz_t n){//crea array di numeri da setacciare per vedere quali sono B smooth
	if(a==NULL || b==NULL || n==NULL || mpz_sgn(n)<=0 || M<=0 || mpz_sgn(a)<=0 || mpz_sgn(b)<=0){
		handle_error_with_exit("invalid parameter create_array_of_number\n");
	}
	mpz_t *array=alloc_array_mpz(2*M+1);
	long j=-M;
	long h=1+2*j;
	mpz_t c;
	mpz_t v,v1;
	mpz_t double_a,double_b,t_j;

	mpz_init(c);
	mpz_init(v);
	mpz_init(v1);
	mpz_init(double_a);
	mpz_init(double_b);
	mpz_init(t_j);

	//calcolo di c
	mpz_mul(v,b,b);//v=b^2
	mpz_sub(v,v,n);//v=b^2-n
	if(mpz_divisible_p(v,a)==0){
		handle_error_with_exit("error in mpz_divisible create array_of_number\n");
	}
	mpz_divexact(c,v,a);//c=(b^2-n)/a

	//calcolo di double_a
	mpz_set(double_a,a);//double_a=a
	mpz_mul_ui(double_a,double_a,2);//double_a=2*a
	
	//calcolo di double_b
	mpz_set(double_b,b);//double_b=b
	mpz_mul_ui(double_b,double_b,2);//double_b=2*b
	
	//imposta primo valore dell'array
	mpz_mul_si(v,a,j*j);//v=a*j^2
	mpz_set_si(v1,2*j);//v1=2*j
	mpz_mul(v1,v1,b);//v1=2*j*b
	mpz_add(v,v,v1);//v=a*j^2+2*j*b
	mpz_add(v,v,c);//v=a*j^2+2*j*b+c
	//array[i]=a*j^2+2*j*b+c;//a(j)=a*(j^2)+2*bj+c se a=1 b=x0 a(j)=(x0+j)^2-n
	mpz_set(array[0],v);//array[i]=a*j^2+2*j*b+c;
	
	//calcolo di t_j=a*(1+2*j)
	mpz_mul_si(t_j,a,h);//t_j=a*(1+2*j)
	
	//imposta valori successivi
	for(int i=1;i<2*M+1;i++){
		mpz_add(v,v,t_j);//v=v+t_j
		mpz_add(v,v,double_b);//v=v+t_j+2b
		mpz_set(array[i],v);//imposta il valore dell'array
		mpz_add(t_j,t_j,double_a);//calcola il nuovo t_j
	}

	mpz_clear(c);
	mpz_clear(v);
	mpz_clear(v1);
	mpz_clear(double_a);
	mpz_clear(double_b);
	mpz_clear(t_j);
	return array;
}*/
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
    if(solution==NULL || head==NULL || num_B_smooth<=0 || card_f_base<=0 || num_B_smooth<card_f_base || a==NULL || b==NULL || n==NULL){
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
    for(int i=0;i<num_B_smooth;i++){
        if(solution[i]==1){//se solution==1 allora bisogna moltiplicare le relazioni
            mpz_set(square,p->square_relation.square);//moltiplica le radici quadrate dei numeri B-smooth mod n
            mpz_mul(a,a,square);//a=a*square
            mpz_mod(a,a,n);//a=a*square mod n
            struct node_factorization*q=p->square_relation.head_factorization;
            while(q!=NULL){
                exponent=q->exp_of_number;
                index=q->index;
                q=q->next;
                sum_exponent_relation[index]+=exponent;
            }
        }
        p=p->next;//passa alla prossima relazione
    }
    divide_vector_multiple_of_2_by_2(sum_exponent_relation,card_f_base);
    for(int j=0;j<card_f_base;j++){//per ogni primo della factor base
        if(sum_exponent_relation[j]==0){
            continue;
        }
        mpz_set_si(v_temp,r.prime[j]);
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
/*int find_factor_of_n_from_base_matrix(int **base_matrix,int num_row,int* num_column,int **matrix_linear_system,int num_row_matrix,int num_col_matrix,mpz_t n,mpz_t**matrix_B_smooth,int num_B_smooth,int card_f_base){
//calcola tutte le soluzioni una per una e per ognuna prova a vedere se trova una coppia a,b che fattorizza n se riesce allora ritorna al chiamante
	if(base_matrix==NULL || *base_matrix==NULL || num_row<=0 || *num_column<=0 || matrix_linear_system==NULL || *matrix_linear_system==NULL
		|| num_row_matrix<=0 || num_col_matrix<=0 || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || card_f_base<=0){
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
		//ciò equivale ad allocare 2^10 soluzioni invece di 2^num_col,ciò in termini di soluzione corrisponde a considerare tutti i vettori
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
		if(check_if_array_is_reduce_mod_n(solution,num_row,2)==0){
			handle_error_with_exit("error in calculate solution\n");
		}
		if(verify_solution(matrix_linear_system,num_row_matrix,num_col_matrix,solution)==0){
			handle_error_with_exit("invalid solution in find_factor_of_n_from_base_matrix\n");
		}
		calculate_a_and_b(solution,matrix_B_smooth,num_B_smooth,card_f_base,a,b,n_copy);//calcola un a e un b
		factor_founded_from_a_and_b=try_to_factor(a,b,n_copy,factor1,factor2);
		if(factor_founded_from_a_and_b>0){
			fprintf(file_log,"p=");
			mpz_out_str(file_log,10,factor1);
			fprintf(file_log," ");
			fprintf(file_log,"q=");
			mpz_out_str(file_log,10,factor2);
			fprintf(file_log," ");
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
	
}*/
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
        if(check_if_array_is_reduce_mod_n(solution,num_row,2)==0){
            handle_error_with_exit("error in calculate solution\n");
        }
        if(verify_solution_char(matrix_linear_system,num_row_matrix,num_col_matrix,solution)==0){
            handle_error_with_exit("invalid solution in find_factor_of_n_from_base_matrix\n");
        }
        calculate_a_and_b_siqs(solution,head,num_B_smooth,card_f_base,a,b,n_copy);//calcola un a e un b
        factor_founded_from_a_and_b=try_to_factor(a,b,n_copy,factor1,factor2);
        if(factor_founded_from_a_and_b>0){
            /*fprintf(file_log,"p=");
            mpz_out_str(file_log,10,factor1);
            fprintf(file_log," ");
            fprintf(file_log,"q=");
            mpz_out_str(file_log,10,factor2);
            fprintf(file_log," ");*/
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
/*int**generate_all_solution(int **base_matrix,int num_row,int* num_column,int **matrix_linear_system,int num_row_matrix,int num_col_matrix){//row della base matrix=num_b_smooth,num col della base matrix=num_b_smooth-rango_matrice
	//base_matrix=v1,v2,vn vettori colonna

	//matrix solution= sol1  vettori riga
	//		sol2
	//		...
	//		sol n
	
	if(*base_matrix==NULL || base_matrix==NULL || num_row<=0 || *num_column<=0){
		handle_error_with_exit("error in parameter generate_all_solution\n");
	}
	int**solution,num_col=*num_column;
	char*combination;
	int j;
	int count_combination=0;
	int*col_h=NULL,*vector_sum=NULL;
	if(num_col>MAX_DIM_SOL){//troppe soluzioni portano ad una allocazione troppo grande della memoria considerare solo i primi 10 vettori della base
		//ciò equivale ad allocare 2^10 soluzioni invece di 2^num_col,ciò in termini di soluzione corrisponde a considerare tutti i vettori
		//non allocati ==0,quindi trova soluzioni riga del tipo=v1,v2,...v10,0,0,0,0,0...0
		num_col=MAX_DIM_SOL;
		*num_column=num_col;
	}
	solution=alloc_matrix_int((pow(2,num_col)-1),num_row);//alloca tutta la memoria per contenere tutte le soluzioni
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
			vector_sum=sum_vector(solution[count_combination],col_h,num_row,num_row);//somma tutti gli array che hanno 				indice uguale a 1 nell'array combination
			memcpy(solution[count_combination],vector_sum,sizeof(int)*num_row);
			free(vector_sum);
			free(col_h);
		}
		count_combination++;
	}
	free(combination);
	reduce_matrix_mod_n(solution,(pow(2,num_col)-1),num_row,2);//riduce la matrice soluzione mod 2
	if(check_if_matrix_is_reduce_mod_n(solution,(pow(2,num_col)-1),num_row,2)==0){
			handle_error_with_exit("error in main,calculate matrix_solution\n");
	}
	for(int i=0;i<(pow(2,num_col)-1);i++){
		if(verify_solution(matrix_linear_system,num_row_matrix,num_col_matrix,solution[i])==0){
			handle_error_with_exit("invalid solution in generate all solution\n");
		}
	}
	return solution;
}

int factor_n(int **matrix_solution,int dim_sol,int num_B_smooth,mpz_t**matrix_B_smooth,int card_f_base,const mpz_t n){//numero righe=numero sol numero colonne=numero b_smooth
	//matrice soluzione=matrice i cui vettori sono soluzioni non banali del sistema lineare(combinazioni lineari dei vettori della base
	mpz_t a,b,factor1,factor2;
	int factorizations_founded=0;//conta il numero di fattorizzazioni trovate con tutti gli a e tutti i b
	char factor_founded_from_a_and_b=0,writed=0;//conta il numero di fattorizzazioni trovate con un a ed un b
	if(n==NULL || mpz_sgn(n)<0 || dim_sol<=0 || matrix_solution==NULL || 
	*matrix_solution==NULL || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0 || 
	card_f_base<=0 || num_B_smooth<card_f_base){
		handle_error_with_exit("error in parameter factor n\n");
	}
	mpz_t n_copy;
	mpz_init(n_copy);
	mpz_init(a);
	mpz_init(b);
	mpz_init(factor1);
	mpz_init(factor2);
	mpz_set(n_copy,n);//n_copy=n
	long num_sol=pow(2,dim_sol)-1;//numero soluzioni==2^n-1 dove n è il numero di vettori della base linearmente indipendenti
	for(long i=0;i<num_sol;i++){
		calculate_a_and_b(matrix_solution[i],matrix_B_smooth,num_B_smooth,card_f_base,a,b,n_copy);//calcola 1 a e un b
		factor_founded_from_a_and_b=try_to_factor(a,b,n_copy,factor1,factor2);
		if(factor_founded_from_a_and_b>0 && writed==0){
			writed=1;
			fprintf(file_log,"p=");
			mpz_out_str(file_log,10,factor1);
			fprintf(file_log," ");
			fprintf(file_log,"q=");
			mpz_out_str(file_log,10,factor2);
			fprintf(file_log," ");
		}
		factorizations_founded+=factor_founded_from_a_and_b;
	}
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(factor1);
	mpz_clear(factor2);
	mpz_clear(n_copy);
	return factorizations_founded;
}*/










