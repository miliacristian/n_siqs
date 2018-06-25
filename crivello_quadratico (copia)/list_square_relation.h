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