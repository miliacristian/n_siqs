#include "matrix_function.h"

extern struct row_factorization r;
extern struct timespec timer;
extern struct timespec time_start;








































char value_is_in_sorted_array(int index_of_prime,struct a_struct*array_a_struct,int length){
    if(array_a_struct==NULL || length<0){
        handle_error_with_exit("error in value is in sorted_array\n");
    }
    //se è minore del primo o maggiore dell'ultimo allora non è nell'array
    if(index_of_prime<array_a_struct[0].index_prime_a || index_of_prime>array_a_struct[length-1].index_prime_a){
        return 0;
    }
    //se sono uguali ritorna 1 altrimenti 0
    for(int i=0;i<length;i++){
        if(index_of_prime==array_a_struct[i].index_prime_a){
            return 1;
        }
    }
    return 0;
}

//scriverla in modo parallelizzato
void create_row_factorization(struct node_factor_base*head_f_base_f,int card_f_base,const mpz_t a,struct a_struct*array_a_struct,int s){
	if(head_f_base_f==NULL || card_f_base<=0 || (array_a_struct==NULL && s>0) || s<0){
		handle_error_with_exit("error in create row factorization\n");
	}
	mpz_t temp;
	int index=0;
	mpz_init(temp);
	struct node_factor_base*p=head_f_base_f;
	for(int i=0;i<card_f_base;i++){
        r.prime[i]=p->prime;//metti il primo della factor base in posizione pari
		r.root_n_mod_p[i]=p->root_n_mod_prime;
		r.root2_n_mod_p[i]= r.prime[i]-r.root_n_mod_p[i];//rad2=p-rad1 mod p
		if((i==0 && s>0) || (i==1 && s>0) ||
		   (index<s && i==array_a_struct[index].index_prime_a) ) {//p divide a non esiste inverso modulo p
			r.inverse_a_mod_p[i]=-1;
            if(i!=0 && i!=1) {
				index++;
			}
        }
        else if(s==0){//a=1
			r.inverse_a_mod_p[i]=1;
		}
        else{//p non divide a ->esiste l'inverso
		    mpz_set_si(temp,r.prime[i]);//temp=p
            if(mpz_invert(temp,a,temp)==0){//temp=inverse of a mod p
				if(i==array_a_struct[index].index_prime_a){
					handle_error_with_exit("error in inverse\n");
				}
				printf("i=%d,s=%d\n",i,s);
            	handle_error_with_exit("error in mpz_invert create_row_factorization\n");
            }
			r.inverse_a_mod_p[i]=mpz_get_si(temp);
		}
		if(p->prime!=-1){
			r.log_prime[i]=(int)round(log2f((float)p->prime));
		}
		p=p->next;
	}
	if(index!=s){
	    handle_error_with_exit("error in initialize inverse mod p create_row_factorization\n");
	}
	mpz_clear(temp);
	return;
}

mpz_t*create_array_temp_factorization(int card_f_base,mpz_t*row_matrix_b_smooth){// crea array -1 0 2 0 30...primo factor base 0...di 2*cardinalità 		della factor base elementi
	if(card_f_base<=0 || row_matrix_b_smooth==NULL){
		handle_error_with_exit("error in parameter create_array_temp_factorization\n");
	}
	mpz_t*v;
	v=alloc_array_mpz(card_f_base*2);//alloca memoria
	for(int i=0;i<card_f_base*2;i=i+2){//ciclo a salti di 2
		mpz_set(v[i],row_matrix_b_smooth[i]);//sui posti dispari (indice pari) viene messo il primo della factor base
	}
	return v;
}

int**create_linear_system(int cardinality_f_base,mpz_t** matrix_B_smooth,int num_B_smooth){//alloca matrice con la fattorizzazione dei numeri b_smmoth
	if(cardinality_f_base<=0 || matrix_B_smooth==NULL || *matrix_B_smooth==NULL || num_B_smooth<=0){
		handle_error_with_exit("error in create_linear_system\n");
	}
	int**linear_system=alloc_linear_system(cardinality_f_base,num_B_smooth);
	for(int i=0;i<cardinality_f_base;i++){
		for(int j=0;j<num_B_smooth;j++){
			if(mpz_divisible_2exp_p(matrix_B_smooth[j][2*i+1],1)==0){//non è divisibile per 2
				linear_system[i][j]=1;
			}
			else{//è divisibile per 2
				linear_system[i][j]=0;
			}
			//linear system[i][j]=esponente del'iesimo elemento della factor base preso dalla fattorizzazione j-esimo numero B-smooth
		}	
	}
	return linear_system;
}



char*create_linear_system_f(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth){
	if(head==NULL || cardinality_factor_base<=0 || num_B_smooth<=0){
		handle_error_with_exit("error in create_linear_system\n");
	}
	char*linear_system=alloc_array_char(cardinality_factor_base*num_B_smooth);
	struct node_square_relation*p=head;
	int col_index=0;
	long index;
	while(p!=NULL){//cicla su tutte le relazioni quadratiche
		struct node_factorization*sq=p->square_relation.head_factorization;
		while(sq!=NULL){//cicla su tutti i fattori della singola relazione quadratica
			if((sq->exp_of_number & 1)!=0){//se l'esponente non è divisibile per 2
				//metti 1 in posizione indice nella colonna iesima
				index=get_index(sq->index,col_index,num_B_smooth);
				linear_system[index]=1;
			}
			sq=sq->next;
		}
		p=p->next;//passa alla prossima relazione quadratica
		col_index++;
	}
	if(col_index!=num_B_smooth){
		handle_error_with_exit("error initialize_linear_system\n");
	}

	return linear_system;
}
unsigned long**create_binary_linear_system(struct node_square_relation*head,int cardinality_factor_base,int num_B_smooth,int*num_col_binary_matrix){
    if(head==NULL || cardinality_factor_base<=0 || num_B_smooth<=0 || num_col_binary_matrix==NULL){
        handle_error_with_exit("error in create_linear_system\n");
    }
    int remainder=num_B_smooth%BIT_OF_UNSIGNED_LONG;
    if(remainder==0){
        *num_col_binary_matrix=num_B_smooth/BIT_OF_UNSIGNED_LONG;
    }
    else{
        *num_col_binary_matrix=((num_B_smooth-remainder)/BIT_OF_UNSIGNED_LONG)+1;
    }
    unsigned long**binary_linear_system=alloc_matrix_unsigned_long(cardinality_factor_base,*num_col_binary_matrix);
    struct node_square_relation*p=head;
    int col_index=0;
    while(p!=NULL){//cicla su tutte le relazioni quadratiche
        struct node_factorization*sq=p->square_relation.head_factorization;
        while(sq!=NULL){//cicla su tutti i fattori della singola relazione quadratica
            if((sq->exp_of_number & 1)!=0){//se l'esponente non è divisibile per 2
                //metti 1 in posizione indice nella colonna iesima
                int remainder=col_index%BIT_OF_UNSIGNED_LONG;//va da 0 a bit_unsigned_long->ci dice qual'è l'indice del bit da mettere a 1
                int index_element=col_index/BIT_OF_UNSIGNED_LONG;
                binary_linear_system[sq->index][index_element]^= 1UL << (BIT_OF_UNSIGNED_LONG-remainder-1);//toggle bit
            }
            sq=sq->next;
        }
        p=p->next;//passa alla prossima relazione quadratica
        col_index++;
    }
    if(col_index!=num_B_smooth){
        handle_error_with_exit("error initialize_linear_system\n");
    }

    return binary_linear_system;
}








