#include "M_and_B_functions.h"
extern double thresold_relation;
extern int num_increment_M_and_B;
extern char combined;
void calculate_best_M_and_B(const mpz_t n,int digit_n,long*M,long*B){
	//calcola valori di M e B cercando di minimizzare il numero di volte che i numeri devono essere aumentati
	if(n==NULL || digit_n<=0 || M==NULL || B==NULL){
		handle_error_with_exit("error in calculate best M and B\n");
	}
	if(digit_n<7){
		*M=25;
		*B=20;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<9){
		*M=68;
		*B=38;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<11){
		*M=125;
		*B=59;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<13){
		*M=303;
		*B=112;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<15){
		*M=492;
		*B=156;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<17){
		*M=1177;
		*B=277;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<19){
		*M=2835;
		*B=502;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<21){
		*M=9922;
		*B=1668;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<23){
		*M=14958;
		*B=2577;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<25){
		*M=22512;
		*B=3540;
        thresold_relation=0.0;
        combined=1;
		return;
	}
	if(digit_n<30){
		*M=22512;
		*B=3540;
		thresold_relation=0.3;
		return;
	}
	if(digit_n<35){
		*M=22512;
		*B=3540;
		thresold_relation=0.3;
		return;
	}
	if(digit_n<40){//fino a 40 cifre meno di 1 secondo
		*M=22512;
		*B=3540;
        thresold_relation=0.3;
		return;
	}
	if(digit_n<45){//1,5 sec
		*M=20000;
		*B=8000;
        thresold_relation=0.4;
		return;
	}
	if(digit_n<50){//2-3-4 sec
		*M=20000;
		*B=19000;
        thresold_relation=0.4;
		return;
	}
	if(digit_n<55){//20 sec
		*M=60000;
		*B=97000;
        thresold_relation=0.6;
		return;
	}
    if(digit_n<60){//50 secondi
        *M=40000;
        *B=130000;
        thresold_relation=0.6;
        return;
    }
	if(digit_n<65){//2 minuti 30 secondi
		*M=35000;
		*B=140000;
		thresold_relation=0.6;
		return;
	}
    if(digit_n<70){//6,5 minuti
        *M=40000;
        *B=180000;
        thresold_relation=0.6;
        return;
    }
    if(digit_n<75){//35 minuti
        *M=45000;
        *B=300000;
        thresold_relation=0.4;
        return;
    }
	if(digit_n<80){
		*M=60000;
		*B=750000;
		thresold_relation=0.4;
		return;
	}
	if(digit_n>=80){
		handle_error_with_exit("error criv_quad >=80\n");
		*M=35000;
		*B=15000;
		return;
	}
	handle_error_with_exit("num not handle\n");
	return;
}
void increment_M_and_B(long*M,long*B){
    //calcola i valori nuovi di M e B secondo una formula scelta opportunamente
    //new_m=(m+perc_m)+(m+perc_m)*perc_m/100
    //new_b=(b+perc_b)+(b+perc_b)*perc_b/100
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