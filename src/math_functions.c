#include "math_functions.h"

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
