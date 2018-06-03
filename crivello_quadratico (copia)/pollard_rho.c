#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic.h"
#include "math_function.h"
long f(long x,long n){
	if(n<=0 || x<0 || x>=n){
		handle_error_with_exit("error in parameter\n");
	}
	long t;
	t=pow(x,2)+1;
	reduce_mod_n(&t,n);
	return t;
}


void pollard_rho(long n){
	long x=2;
	long y=2;
	long d=1;
	while (d==1){
		x=f(x,n);
		y=f(f(y,n),n);
		d=gcd(abs(x-y),n);
	}
	if(d==n){
		handle_error_with_exit("no factorization founded\n");
	}
	else{
		printf("n= %ld*%ld\n",d,n/d);
	}
	return;
}
int main(){
	pollard_rho(47*53);
	return 0;
}