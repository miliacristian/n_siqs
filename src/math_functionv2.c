#include "math_functionv2.h"
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
