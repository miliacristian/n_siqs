#include "mpz_functions.h"
//print mpz
void print_array_mpz(mpz_t*array,int length){
	if(length<=0 || array==NULL){
		handle_error_with_exit("error in print_array_mpz\n");
	}
	if(not_print_array(length)==1){
		return;
	}
	for(int i=0;i<length;i++){
		gmp_printf("%Zd ",array[i]);
	}
	printf("\n");
	return;
}

void print_list_mpz(struct node*head){
	if (head==NULL){
		printf("impossible print list head is NULL\n");
		return;
	}
	struct node *p=head;
	while(p!=NULL){
		gmp_printf("%Zd ",p->prime);
		p=p->next;
	}
	printf("\n");
	return;
}
void print_matrix_mpz(mpz_t**matrix,int num_row,int num_col){
	if(matrix==NULL){
		handle_error_with_exit("pointer null\n");
	}
	if(*matrix==NULL){
		handle_error_with_exit("empty matrix\n");
	}
	if(num_row <=0 || num_col <=0){
		handle_error_with_exit("error in print_matrix_mpz\n");
	}
	if(not_print_matrix(num_row,num_col)==1){
		return;
	}
	for(int i=0;i<num_row;i++){
			print_array_mpz(matrix[i],num_col);
	}
	return;
}
//funzioni per allocare array o matrici 
void copy_array_mpz(mpz_t*array1,mpz_t*array2,int length){
	if(array1==NULL || array2==NULL || length<=0){
		handle_error_with_exit("error in copy_array_mpz\n");
	}
	for(int i=0;i<length;i++){
		mpz_set(array1[i],array2[i]);
	}
	return;
}
void copy_matrix_mpz(mpz_t**matrix1,mpz_t**matrix2,int num_row,int num_col){
	if(matrix1==NULL || *matrix1==NULL || matrix2==NULL || *matrix2==NULL
	|| num_row<=0 || num_col<=0){
		handle_error_with_exit("error in copy_array_mpz\n");
	}
	for(int i=0;i<num_row;i++){
		copy_array_mpz(matrix1[i],matrix2[i],num_col);//copia riga
	}
	return;
}

mpz_t*alloc_array_mpz(int length){
	if(length<=0){
		handle_error_with_exit("error in parameter alloc_array_mpz\n");
	}
	mpz_t *array=malloc(sizeof(mpz_t)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc\n");
	}
	for(int i=0;i<length;i++){
		mpz_init(array[i]);
	}
	return array;
}

mpz_t***alloc_array_matrix_mpz(int length_array_matrix){
	mpz_t***array_matrix=NULL;
	if(length_array_matrix<=0){
		handle_error_with_exit("error in parameter alloc aray_matrix_mpz\n");
	}
	array_matrix=malloc(sizeof(mpz_t**)*length_array_matrix);
	if(array_matrix==NULL){
		handle_error_with_exit("error in malloc alloc array_matrix_mpz\n");
	}
	memset(array_matrix,0,sizeof(mpz_t**)*length_array_matrix);
	return array_matrix;
}

mpz_t**realloc_array_pointer_to_mpz(mpz_t**t,int new_length){
	if(new_length<=0 || t==NULL || *t==NULL){
		handle_error_with_exit("error in parameteralloc array pointer to mpz\n");
	}
	mpz_t **array=realloc(t,sizeof(mpz_t*)*(new_length));
	if(array==NULL){
		handle_error_with_exit("error in realloc alloc array pointer to long\n");
	}
	return array;
}

mpz_t**alloc_array_pointer_to_mpz(int length){
	if(length<=0){
		handle_error_with_exit("error in parameteralloc array pointer to mpz\n");
	}
	mpz_t **array=malloc(sizeof(mpz_t*)*(length));
	if(array==NULL){
		handle_error_with_exit("error in malloc alloc array pointer to long\n");
	}
	memset(array,0,sizeof(mpz_t*)*(length));
	return array;
}

mpz_t **alloc_matrix_mpz(int num_row,int num_col){
	if(num_row<=0 || num_col <=0){
		handle_error_with_exit("error in alloc_matrix_long\n");
	}
	mpz_t**matrix;
	matrix=alloc_array_pointer_to_mpz(num_row);
	for(int i=0;i<num_row;i++){
		matrix[i]=alloc_array_mpz(num_col);
	}
	return matrix;
}
void copy_column_matrix_mpz(mpz_t *v,mpz_t**matrix_factorization,long num_row,int index_col){
	if(v==NULL || matrix_factorization==NULL || *matrix_factorization==NULL || num_row<=0 || index_col<0){
		handle_error_with_exit("error in copy_column matrix mpz\n");
	}
	for(long i=0;i<num_row;i++){
		mpz_set(v[i],matrix_factorization[i][index_col]);
	}
	return;
}

void sum_elem_multiple_of_2_mpz(mpz_t*vector1,mpz_t*vector2,int length1,int length2){//vettore risultato=vector1=copia elementi di posto dispari(indice pari) e somma elementi di posto pari(indice dispari)
//vector1 e vector2 devono avere gli elementi di posto dispari uguali e la loro lunghezza deve essere pari
	if(vector1==NULL || vector2==NULL || length1<=0 || length2<=0 || length1!=length2 || length1%2!=0 ){
		handle_error_with_exit("invalid parameter sum_elem_multiple_of_2\n");
	}
	mpz_t t;
	mpz_init(t);
	for(int i=0;i<length1;i=i+2){
		mpz_add(t,vector1[i+1],vector2[i+1]);
		mpz_set(vector1[i+1],t);
	}
	mpz_clear(t);
	return;
}

void find_max_array_mpz(mpz_t max,mpz_t*array,long length){
	if(array==NULL || max==NULL || length<=0){
		handle_error_with_exit("error in sort_array_number\n");
	}
	mpz_set(max,array[0]);
	if(length==1){
		return;
	}
	for(long i=1;i<length;i++){
		if(mpz_cmp(max,array[i])<0){//se il valore nell'array è maggiore prendi quello
			mpz_set(max,array[i]);
		}
	}
	return;
}


void divide_elem_multiple_of_2_by_x(mpz_t*vector,int length,double x){//vector = num1 exp1 num2 exp2 ecc,divide elementi di posto pari per 2
	//il vettore deve avere lunghezza pari
	if(vector==NULL || length<=0 || length%2!=0 || x==0){
		handle_error_with_exit("divide_elem_multiple_of_2_by_2\n");
	}
	if(x!=2){
		handle_error_with_exit("error in parameter\n");
	}
	for(int i=1;i<length;i=i+2){
		if(mpz_divisible_ui_p(vector[i],(unsigned long)x)==0){
			handle_error_with_exit("error in mpz_divisible,divide_elem_multiple_of_2_by_2\n");
		}
		mpz_divexact_ui(vector[i],vector[i],x);
	}
	return;
}
void free_memory_array_mpz(mpz_t*array,long length){
	if(length<0){
		handle_error_with_exit("error in free memory array mpz\n");
	}
	if(length==0){
 		if(array==NULL){
			return;
		}
		else{
			handle_error_with_exit("error in parameter array free memory array mpz\n");
		}
	}
	for(long i=0;i<length;i++){
		mpz_clear(array[i]);
	}
	free(array);
}

void free_memory_matrix_mpz(mpz_t**matrix,int num_row,int num_col){
	if(matrix==NULL || *matrix==NULL || num_row<=0 || num_col<=0){
		handle_error_with_exit("error in free memory matrix mpz\n");
	}
	for(int i=0;i<num_row;i++){
		free_memory_array_mpz(matrix[i],num_col);
	}
	free(matrix);
	return;
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
void print_array_matrix_same_dimension(mpz_t***array_matrix_mpz,int length_array,int num_row,int num_col){
    if(array_matrix_mpz==NULL || *array_matrix_mpz==NULL || **array_matrix_mpz==NULL || length_array<=0
       || num_row<=0 ||  num_col<=0){
        handle_error_with_exit("error in print_array_mpz_same_dimension\n");
    }
    for(int i=0;i<length_array;i++){
        print_matrix_mpz(array_matrix_mpz[i],num_row,num_col);
    }
    return;
}

