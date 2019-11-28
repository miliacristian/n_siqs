#include "square_relations_functions.h"

extern struct row_factorization r;

void free_memory_list_square_relation(struct node_square_relation*head){
    if(head==NULL){
        return;
    }
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL){
        q=p->next;
        mpz_clear(p->square_relation.square);
        mpz_clear(p->square_relation.num);
        mpz_clear(p->square_relation.residuos);
        free_list_factorization(p->square_relation.head_factorization);
        free(p);
        p=q;
    }
    return;
}

void print_struct_square_relation(struct square_relation square_relation){
    gmp_printf("square=%Zd,residuos=%Zd\n",square_relation.square,square_relation.residuos);
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
void print_list_square_relation(struct node_square_relation*head,int length){
    if (head==NULL){
        printf("impossible print list head is NULL\n");
        return;
    }
    if(length<0){
        handle_error_with_exit("error in print_list_square_relation\n");
    }
    if(length==0){
        printf("list is empty\n");
    }
    if(not_print_list(length)==1){
        return;
    }
    struct node_square_relation*p=head;
    while(p!=NULL){
        gmp_printf("square=%Zd,residuos=%Zd,",p->square_relation.square,p->square_relation.residuos);
        print_factorization(p->square_relation.num,p->square_relation.head_factorization);
        p=p->next;
    }
    printf("\n");
    return;
}
char verify_square_relation(struct square_relation square_relation,const mpz_t n){
    if(TEST==0){
        return 1;
    }
    struct node_factorization*head_factor=square_relation.head_factorization;
    mpz_t temp,num_temp,square;
    mpz_init(num_temp);
    mpz_init(temp);
    mpz_init(square);

    if(calculate_num_from_factorization(num_temp,head_factor)==0){
        handle_error_with_exit("error in calculate_num_from_factorization\n");
    }
    mpz_mod(num_temp,num_temp,n);//num modulo n
    mpz_mul(square,square_relation.square,square_relation.square);//square=x^2
    mpz_mod(square,square,n);
    if(mpz_cmp(num_temp,square)!=0){//non sono congrui modulo n
        handle_error_with_exit("error in verify square relation\n");
    }
    mpz_clear(temp);
    mpz_clear(num_temp);
    mpz_clear(square);
    return 1;
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