#include "basic.h"
#ifndef CRIVELLO_QUADRATICO_COPIA_1_LIST_SQUARE_RELATION_H
#define CRIVELLO_QUADRATICO_COPIA_1_LIST_SQUARE_RELATION_H
struct node_square_relation {
    struct square_relation square_relation;
    struct node_square_relation*next;
    struct node_square_relation*prev;
};
#endif //CRIVELLO_QUADRATICO_COPIA_1_LIST_SQUARE_RELATION_H
void insert_ordered_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail);
void union_list_square(struct node_square_relation**head,struct node_square_relation**tail,struct node_square_relation*head_square,struct node_square_relation*tail_square);
void remove_same_num(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth);
void insert_ordered_residuos_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail);
void combine_relation_B_smooth_and_semi_B_smooth(struct node_square_relation*head,struct node_square_relation*tail,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,mpz_t n,int*num_B_smooth);
void remove_same_square(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth);
char verify_sorted_num_square_rel_list(struct node_square_relation*head);
void insert_ordered_num_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail);
void remove_same_square_and_sort_by_num(struct node_square_relation*head,struct node_square_relation*tail,struct node_square_relation**head_final,struct node_square_relation**tail_final,int*num_B_smooth,int*num_semi_B_smooth);
char verify_sorted_square_rel_list(struct node_square_relation*head);
char verify_sorted_residuos_square_rel_list(struct node_square_relation*head);