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
char combine_relation_B_smooth_and_semi_B_smooth(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,struct node_square_relation*head_sort_residuos,mpz_t n,int*num_B_smooth,int*num_semi_B_smooth,int*combined_relations);
void remove_same_square(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth);
char verify_sorted_num_square_rel_list(struct node_square_relation*head);
void insert_ordered_num_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail);
void remove_same_square_and_sort_by_num(struct node_square_relation*head,struct node_square_relation*tail,struct node_square_relation**head_final,struct node_square_relation**tail_final,int*num_B_smooth,int*num_semi_B_smooth);
char verify_sorted_square_rel_list(struct node_square_relation*head);
char verify_sorted_residuos_square_rel_list(struct node_square_relation*head);
void sort_relation_by_residuos(struct node_square_relation*head,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos);
void sort_relation_by_num(struct node_square_relation*head,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos);
void union_list_residuos(struct node_square_relation**head_residuos,struct node_square_relation**tail_residuos,struct node_square_relation*phead_residuos,struct node_square_relation*ptail_residuos);
void add_square_relation_to_list_sorted(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,struct node_square_relation*head_square);
void add_relation_semi_B_smooth_to_list(struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,struct node_square_relation*head_residuos);
void insert_at_tail_square_rel(struct node_square_relation *new_node,struct node_square_relation**head,struct node_square_relation** tail);
void insert_at_tail_square_relation(struct square_relation relation,struct node_square_relation**head,struct node_square_relation** tail);
char verify_cardinality_list_square_relation(struct node_square_relation*head,int length);
void union_list_residuos_v2(struct node_square_relation**head_residuos1,struct node_square_relation**tai_residuos1,struct node_square_relation**head_residuos2,struct node_square_relation**tail_residuos2);
void union_list_square_v2(struct node_square_relation**head_square1,struct node_square_relation**tail_square1,struct node_square_relation**head_square2,struct node_square_relation**tail_square2);
char combine_relation_B_smooth_and_semi_B_smooth_v2(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,
                                                    struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,mpz_t n,int*num_B_smooth,int*num_semi_B_smooth,int*combined_relations);
void swap_square_relation(struct square_relation *pnode1,struct square_relation *pnode2);
struct node_square_relation *lastNode(struct node_square_relation *root);
struct node_square_relation* partition(struct node_square_relation *l,struct node_square_relation *h);
void _quickSort(struct node_square_relation* l,struct node_square_relation *h);
void quickSort(struct node_square_relation *head);