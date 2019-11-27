#include "square_relations_functions.h"

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
