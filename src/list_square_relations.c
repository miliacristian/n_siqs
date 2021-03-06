#include "list_square_relations.h"

// A utility function to find last node of linked list
void swap_square_relation(struct square_relation *pnode1,struct square_relation *pnode2){
    struct square_relation temp = *pnode1;
    *pnode1 = *pnode2;
    *pnode2 = temp;
    return;
}

/* Considers last element as pivot, places the pivot element at its
   correct position in sorted array, and places all smaller (smaller than
   pivot) to left of pivot and all greater elements to right of pivot */
struct node_square_relation *lastNode(struct node_square_relation *root)
{
    while (root && root->next)
        root = root->next;
    return root;
}

struct node_square_relation* partition_residuos(struct node_square_relation *l,struct node_square_relation *h)
{
    // set pivot as h element
    struct square_relation x  = h->square_relation;

    // similar to i = l-1 for array implementation
    struct node_square_relation *i = l->prev;

    // Similar to "for (int j = l; j <= h- 1; j++)"
    for (struct node_square_relation *j = l; j != h; j = j->next)
    {
        if (mpz_cmp(j->square_relation.residuos,x.residuos)<=0)
        {
            // Similar to i++ for array
            i = (i == NULL)? l : i->next;
            swap_square_relation(&(i->square_relation), &(j->square_relation));
        }
    }
    i = (i == NULL)? l : i->next; // Similar to i++
    swap_square_relation(&(i->square_relation), &(h->square_relation));
    return i;
}
struct node_square_relation* partition_square(struct node_square_relation *l,struct node_square_relation *h)
{
    // set pivot as h element
    struct square_relation x  = h->square_relation;

    // similar to i = l-1 for array implementation
    struct node_square_relation *i = l->prev;

    // Similar to "for (int j = l; j <= h- 1; j++)"
    for (struct node_square_relation *j = l; j != h; j = j->next)
    {
        if (mpz_cmp(j->square_relation.square,x.square)<=0)
        {
            // Similar to i++ for array
            i = (i == NULL)? l : i->next;
            swap_square_relation(&(i->square_relation), &(j->square_relation));
        }
    }
    i = (i == NULL)? l : i->next; // Similar to i++
    swap_square_relation(&(i->square_relation), &(h->square_relation));
    return i;
}
/* A recursive implementation of quicksort for linked list */
void _quickSort_square(struct node_square_relation* l,struct node_square_relation *h)
{
    if (h != NULL && l != h && l != h->next)
    {
        struct node_square_relation *p = partition_square(l, h);
        _quickSort_square(l, p->prev);
        _quickSort_square(p->next, h);
    }
}
void _quickSort_residuos(struct node_square_relation* l,struct node_square_relation *h)
{
    if (h != NULL && l != h && l != h->next)
    {
        struct node_square_relation *p = partition_residuos(l, h);
        _quickSort_residuos(l, p->prev);
        _quickSort_residuos(p->next, h);
    }
}
// The main function to sort a linked list. It mainly calls _quickSort()
void quickSort_square(struct node_square_relation *head)
{
    // Find last node
    struct node_square_relation *h = lastNode(head);

    // Call the recursive QuickSort
    _quickSort_square(head, h);
}
void quickSort_residuos(struct node_square_relation *head)
{
    // Find last node
    struct node_square_relation *h = lastNode(head);

    // Call the recursive QuickSort
    _quickSort_residuos(head, h);
}
#if DEBUG==1
char verify_sorted_num_square_rel_list(struct node_square_relation*head){
    if(head==NULL){
        return 1;
    }
    struct node_square_relation*p=head;
    while(p!=NULL && p->next!=NULL) {
        if(mpz_cmp(p->square_relation.num,p->next->square_relation.num)>0){
            return 0;
        }
        p = p->next;
    }
    return 1;
}
char verify_cardinality_list_square_relation(struct node_square_relation*head,int length){
    if(length<0){
        printf("length minore di zero\n");
        return 0;
    }
    int counter=0;
    if(head==NULL && length==0){
        return 1;
    }
    struct node_square_relation*p=head;
    while(p!=NULL){
        counter++;
        p=p->next;
    }
    if(counter!=length){
        printf("elementi contati=%d\n",counter);
        printf("lunghezza errata\n");
        return 0;
    }
    return 1;
}
char verify_sorted_residuos_square_rel_list(struct node_square_relation*head){
    if(head==NULL){
        return 1;
    }
    struct node_square_relation*p=head;
    while(p!=NULL && p->next!=NULL) {
        if(mpz_cmp(p->square_relation.residuos,p->next->square_relation.residuos)>0){
            return 0;
        }
        p = p->next;
    }
    return 1;
}
char verify_sorted_square_rel_list(struct node_square_relation*head){
    if(head==NULL){
        return 1;
    }
    struct node_square_relation*p=head;
    while(p!=NULL && p->next!=NULL) {
        if(mpz_cmp(p->square_relation.square,p->next->square_relation.square)>0){
            return 0;
        }
        p = p->next;
    }
    return 1;
}
#endif
int count_element_linked_list_square_rel(struct node_square_relation*head){
    int count=0;
    if (head==NULL){
        handle_error_with_exit("head is NULL\n");
    }
    struct node_square_relation *p=head;
    while(p!=NULL){
        count++;
        p=p->next;
    }
    return count;
}

void remove_after_node_square_rel(struct node_square_relation**ppos,struct node_square_relation**tail){
    if(ppos==NULL || tail==NULL ){
        handle_error_with_exit("error in parameter remove_after_node\n");
    }
    if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
    }
    struct node_square_relation *r = *ppos;
    struct node_square_relation*q=r->next;
    if(q==NULL){//fine della lista,bisogna aggiornare la coda
        *ppos=r->next;
        *tail=r->prev;
        mpz_clear(r->square_relation.num);
        mpz_clear(r->square_relation.square);
        mpz_clear(r->square_relation.residuos);
        free_memory_list_factor((r)->square_relation.head_factorization);
        free(r);
        r=NULL;
    }
    else {
        q->prev = r->prev;
        *ppos = r->next;
        mpz_clear(r->square_relation.num);
        mpz_clear(r->square_relation.square);
        mpz_clear(r->square_relation.residuos);
        free_memory_list_factor((r)->square_relation.head_factorization);
        free(r);
        r=NULL;
    }
    return;
}

int delete_head_square_rel(struct node_square_relation** head){//non è importante il valore iniziale di oldhead
    //initializza oldhead con il primo nodo della lista e distrugge il primo nodo della lista
    if(head==NULL){
        handle_error_with_exit("error in delete head\n");
    }
    if(*head == NULL){//nessuna testa da eliminare
        return -1;
    }
    if ((*head)-> next == NULL){//c'è un solo nodo in lista
        free(*head);
        *head = NULL;
    }else{
        struct node_square_relation*temp=(*head)->next;
        free_memory_list_factor((*head)->square_relation.head_factorization);
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first_square_rel(struct node_square_relation *new_node, struct node_square_relation **head, struct node_square_relation **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail_square_rel(struct node_square_relation *new_node,struct node_square_relation**head,struct node_square_relation** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}


//alloca e inizializza un nodo della lista dinamica ordinata
struct node_square_relation* get_new_node_square_rel(struct square_relation square_relation) {
    if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in parameter get_new_node\n");
    }
    struct node_square_relation* new_node = (struct node_square_relation*)malloc(sizeof(struct node_square_relation));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    new_node->square_relation=square_relation;//copy field element by element
    new_node->prev = NULL;
    new_node->next = NULL;
    return new_node;
}
void insert_at_tail_square_relation(struct square_relation relation,struct node_square_relation**head,struct node_square_relation** tail){//inserisce un nodo in coda
    #if DEBUG==1
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    #endif
    struct node_square_relation*new_node=get_new_node_square_rel(relation);
    if(*head == NULL) {
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}

void insert_at_head_square_rel(struct node_square_relation* new_node,struct node_square_relation** head,struct node_square_relation** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller_num_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.num,node2.square_relation.num)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}
char first_is_smaller_sort_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.square,node2.square_relation.square)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}
char first_is_smaller_residuos_square_rel(struct node_square_relation node1, struct node_square_relation node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.square_relation.residuos,node2.square_relation.residuos)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}

void sort_relation_by_residuos(struct node_square_relation*head,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos){
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL) {
        q = p->next;
        insert_ordered_residuos_square_rel(p->square_relation, head_sort_residuos, tail_sort_residuos);
        free(p);
        p = q;
    }
    return;
}
void sort_relation_by_num(struct node_square_relation*head,struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos){
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL) {
        q = p->next;
        insert_ordered_num_square_rel(p->square_relation, head_sort_residuos, tail_sort_residuos);
        free(p);
        p = q;
    }
    return;
}
struct square_relation create_relation_large_prime(struct square_relation rel1,struct square_relation rel2,mpz_t n,char*factorization_founded){
    struct square_relation new_relation;
    int number,exp_of_number,index;
    mpz_t temp,temp2;
    struct node_factorization*tail=NULL;
    if(factorization_founded==NULL){
        handle_error_with_exit("error in create relation_large prime factorization_founded is NULL\n");
    }
    if(mpz_cmp(rel1.residuos,rel2.residuos)!=0){
        handle_error_with_exit("error in create_relation_large_prime residuos are different\n");
    }
    else if(mpz_cmp_si(rel1.residuos,1)==0){
        handle_error_with_exit("error in create_relation_large_prime residuos are equal to one\n");
    }
    mpz_init(temp);
    mpz_init(temp2);

    mpz_init(new_relation.square);
    mpz_init(new_relation.residuos);
    mpz_init(new_relation.num);

    mpz_set_si(new_relation.residuos,1);//residuo=1
    mpz_set_si(new_relation.num,0);//num=0
    new_relation.head_factorization=NULL;

    //new_square=square1*square2*(residuos)^-1
    //square=residuos^-1 mod n
    if(mpz_invert(new_relation.square,rel1.residuos,n)==0){//non esiste l'inverso modulo n,gcd(residuos,n)!=1
        mpz_gcd(temp,rel1.residuos,n);//temp=gcd(residuos,n)
        if(mpz_cmp(temp,n)!=0){//se il gcd(residuos,n)!=n
            mpz_divexact(temp2,n,temp);
            gmp_printf("factorization of n=%Zd*%Zd\n",temp2,temp);
            *factorization_founded=1;
            mpz_clear(temp);
            mpz_clear(temp2);
            mpz_clear(new_relation.square);
            mpz_clear(new_relation.residuos);
            mpz_clear(new_relation.num);
            return new_relation;
        }
        else{//gcd(residuos,n)!=n
            *factorization_founded=-1;
            mpz_clear(temp);
            mpz_clear(temp2);
            mpz_clear(new_relation.square);
            mpz_clear(new_relation.residuos);
            mpz_clear(new_relation.num);
            return new_relation;
        }

    }
    mpz_mul(new_relation.square,new_relation.square,rel1.square);
    mpz_mod(new_relation.square,new_relation.square,n);//square=residuos^-1*square1 mod n
    mpz_mul(new_relation.square,new_relation.square,rel2.square);
    mpz_mod(new_relation.square,new_relation.square,n);//square=residuos^-1*square1*square2 mod n

    //print_struct_square_relation(rel1);
    //print_struct_square_relation(rel2);

    struct node_factorization*p1=rel1.head_factorization;
    struct node_factorization*p2=rel2.head_factorization;
    while(p1!=NULL || p2!=NULL){//finisci quando sei arrivato alla fine di entrambe le liste
        if(p1!=NULL && p2!=NULL && p1->index<p2->index){
            //se le liste non sono finite e indice minore nella lista 1 aggiungi il nodo della lista 1
            // e vai avanti di un elemento
            index=p1->index;
            exp_of_number=p1->exp_of_number;
            number=p1->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p1=p1->next;//vai avanti al nodo successivo nella lista 1
        }
        else if(p1!=NULL && p2!=NULL && p1->index==p2->index){
            //se le liste non sono finite e indici uguali somma gli esponenti e vai avanti in entrambe le liste
            if(p1->number!=-1) {
                index = p1->index;//l'indice è uguale
                exp_of_number = p1->exp_of_number + p2->exp_of_number;
                number = p1->number;//il numero è uguale
                insert_ordered_factor(number, exp_of_number, index, &(new_relation.head_factorization), &tail);
            }
            p1=p1->next;//vai avanti in entrambe le liste
            p2=p2->next;//vai avanti in entrambe le liste
        }
        else if(p1!=NULL && p2!=NULL && p1->index>p2->index){
            //se le liste non sono finite e indice della lista 2 è minore,aggiungi nodo della lista 2 e vai avanti di un nodo nella lista 2
            index=p2->index;
            exp_of_number=p2->exp_of_number;
            number=p2->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p2=p2->next;//vai avanti nella lista 2
        }
        else if(p1==NULL && p2!=NULL){//se la prima lista è finita aggiungi nodo della seconda lista e vai avanti di un nodo nella lista 2
            index=p2->index;
            exp_of_number=p2->exp_of_number;
            number=p2->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p2=p2->next;//vai avanti di un nodo nella lista 2
        }
        else if(p1!=NULL && p2==NULL){//se la seconda lista è finita aggiungi nodo della prima lista
            index=p1->index;
            exp_of_number=p1->exp_of_number;
            number=p1->number;
            insert_ordered_factor(number,exp_of_number,index,&(new_relation.head_factorization),&tail);
            p1=p1->next;//vai avanti di un nodo nella lista 1 vai avanti di un nodo nella lista 1
        }
        else{
            handle_error_with_exit("error in create relation large prime,caso non gestito\n");
        }
    }
    //printf("new relation\n");//
    //print_struct_square_relation(new_relation);//
    mpz_clear(temp);
    mpz_clear(temp2);
    return new_relation;
}
void insert_ordered_sort_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in insert_ordered\n");
    }
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    struct node_square_relation* temp = *tail;
    struct node_square_relation* next_node = NULL;
    struct node_square_relation* new_node = get_new_node_square_rel(square_relation);

    if(*head == NULL){
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    if(first_is_smaller_sort_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_sort_square_rel(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_square_rel(new_node, head, tail);
                return;
            }
        }
        next_node = temp->next;
        new_node->prev = temp;
        new_node->next = next_node;
        temp->next = new_node;
        next_node->prev = new_node;
    }
    return;
}

char combine_relation_B_smooth_and_semi_B_smooth(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,
        struct node_square_relation*head_sort_residuos,mpz_t n,int*num_B_smooth,int*num_semi_B_smooth,int*combined_relations){
    if(head_sort_square==NULL || tail_sort_square==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL || combined_relations==NULL){
        handle_error_with_exit("error in combine_relation_B_smooth and semi_B_smooth\n");
    }
    char factorization_founded=0;
    if(head_sort_residuos==NULL){//nothing to do
        return factorization_founded;
    }
    struct node_square_relation*p=head_sort_residuos;//p=nodo
    struct node_square_relation*q;
    while(p!=NULL) {
        q=p->next;
        //ciclo sulla lista per trovare residui uguali e creare nuove relazioni B_smooth e ordinale per square
        while (p!=NULL && q!=NULL && mpz_cmp(p->square_relation.residuos, q->square_relation.residuos) == 0) {//residui uguali
            if (mpz_cmp(p->square_relation.square, q->square_relation.square) != 0) {
                (*num_B_smooth)++;
                (*combined_relations)++;
                struct square_relation new_square_relation = create_relation_large_prime(p->square_relation,
                                                                                         q->square_relation, n,
                                                                                         &factorization_founded);
                if (factorization_founded == 1) {
                    free_memory_list_square_relation(p);
                    head_sort_residuos = NULL;
                    return factorization_founded;
                } else if (factorization_founded != -1) {
                    insert_ordered_sort_square_rel(new_square_relation, head_sort_square, tail_sort_square);
                    q = q->next;
                    if(verify_square_relation(new_square_relation,n)==0){
                        handle_error_with_exit("error in create relation with combine function\n");
                    }
                }
            }
            else{//relazioni con lo stesso residuo ma con lo stesso quadrato
                q=q->next;
            }
        }
        //una volta che ho creato tutte le relazioni quadratiche sfruttando un numero semi_B_smooth lo tolgo dalla lista
        if(p!=NULL) {
            q = p->next;
            mpz_clear(p->square_relation.square);
            mpz_clear(p->square_relation.num);
            mpz_clear(p->square_relation.residuos);
            free_list_factorization(p->square_relation.head_factorization);
            free(p);
            p = q;
        }
        if(q==NULL){
            return factorization_founded;
        }
    }
    return factorization_founded;
}
char combine_relation_B_smooth_and_semi_B_smooth_v3(struct node_square_relation**head_square,struct node_square_relation**tail_square,
                                                    struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,mpz_t n,int*num_B_smooth,int*num_semi_B_smooth,int*combined_relations){
    if(head_square==NULL || tail_square==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL || combined_relations==NULL || head_sort_residuos==NULL || tail_sort_residuos==NULL){
        handle_error_with_exit("error in combine_relation_B_smooth and semi_B_smooth\n");
    }
    char factorization_founded=0;
    if(*head_sort_residuos==NULL){//nothing to do
        return factorization_founded;
    }
    struct node_square_relation*p=*head_sort_residuos;//p=nodo
    struct node_square_relation*q;
    while(p!=NULL) {
        //ciclo sulla lista per trovare residui uguali e creare nuove relazioni B_smooth e ordinale per square
        while (p!=NULL && p->next!=NULL && mpz_cmp(p->square_relation.residuos, p->next->square_relation.residuos) == 0) {//residui uguali
            if (mpz_cmp(p->square_relation.square, p->next->square_relation.square) != 0) {
                (*num_B_smooth)++;
                (*combined_relations)++;
                struct square_relation new_square_relation = create_relation_large_prime(p->square_relation,
                                                                                         p->next->square_relation, n,
                                                                                         &factorization_founded);
                if (factorization_founded == 1) {
                    free_memory_list_square_relation(p);
                    head_sort_residuos = NULL;
                    return factorization_founded;
                } else if (factorization_founded != -1) {
                    insert_at_tail_square_relation(new_square_relation, head_square, tail_square);
                    remove_after_node_square_rel(&(p->next),tail_sort_residuos);
                    if(verify_square_relation(new_square_relation,n)==0){
                        handle_error_with_exit("error in create relation with combine function\n");
                    }
                }
            }
            else{//relazioni con lo stesso residuo ma con lo stesso quadrato
                remove_after_node_square_rel(&(p->next),tail_sort_residuos);
            }
        }
        //una volta che ho creato tutte le relazioni quadratiche sfruttando un numero semi_B_smooth lo tolgo dalla lista
        q=p->next;
        remove_after_node_square_rel(&p,tail_sort_residuos);
        p=q;
        if(q==NULL){
            return factorization_founded;
        }
    }
    return factorization_founded;
}
char combine_relation_B_smooth_and_semi_B_smooth_v2(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,
                                                 struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,mpz_t n,int*num_B_smooth,int*num_semi_B_smooth,int*combined_relations){
    if(head_sort_square==NULL || tail_sort_square==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL || combined_relations==NULL || head_sort_residuos==NULL || tail_sort_residuos==NULL){
        handle_error_with_exit("error in combine_relation_B_smooth and semi_B_smooth\n");
    }
    char factorization_founded=0;
    if(*head_sort_residuos==NULL){//nothing to do
        return factorization_founded;
    }
    struct node_square_relation*p=*head_sort_residuos;//p=nodo
    struct node_square_relation*q;
    while(p!=NULL) {
        //ciclo sulla lista per trovare residui uguali e creare nuove relazioni B_smooth e ordinale per square
        while (p!=NULL && p->next!=NULL && mpz_cmp(p->square_relation.residuos, p->next->square_relation.residuos) == 0) {//residui uguali
            if (mpz_cmp(p->square_relation.square, p->next->square_relation.square) != 0) {
                (*num_B_smooth)++;
                (*combined_relations)++;
                struct square_relation new_square_relation = create_relation_large_prime(p->square_relation,
                                                                                         p->next->square_relation, n,
                                                                                         &factorization_founded);
                if (factorization_founded == 1) {
                    free_memory_list_square_relation(p);
                    head_sort_residuos = NULL;
                    return factorization_founded;
                } else if (factorization_founded != -1) {
                    insert_ordered_sort_square_rel(new_square_relation, head_sort_square, tail_sort_square);
                    remove_after_node_square_rel(&(p->next),tail_sort_residuos);
                    if(verify_square_relation(new_square_relation,n)==0){
                        handle_error_with_exit("error in create relation with combine function\n");
                    }
                }
            }
            else{//relazioni con lo stesso residuo ma con lo stesso quadrato
                remove_after_node_square_rel(&(p->next),tail_sort_residuos);
            }
        }
        //una volta che ho creato tutte le relazioni quadratiche sfruttando un numero semi_B_smooth lo tolgo dalla lista
        q=p->next;
        remove_after_node_square_rel(&p,tail_sort_residuos);
        p=q;
        if(q==NULL){
            return factorization_founded;
        }
    }
    return factorization_founded;
}
void remove_same_num(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth){
    if(head==NULL || tail==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL){
        handle_error_with_exit("error in remove_same_num\n");
    }
    struct node_square_relation*l=*head;
    while(l!=NULL){
        while(l->next!=NULL && (mpz_cmp(l->square_relation.num,(l->next)->square_relation.num)==0) && (mpz_cmp_si(l->square_relation.num,0)!=0 )) {
            if (mpz_cmp_si(l->square_relation.residuos, 1) == 0){//trovati due numeri uguali con residuo uguale a 1
                (*num_B_smooth)--;
            }
            else{//trovati due numeri uguali con residuo diverso da 1
                (*num_semi_B_smooth)--;
            }
            remove_after_node_square_rel(&(l->next),tail);
        }
        if(l->next==NULL){
            return;
        }
        l=l->next;
    }
    return;
}
int remove_same_square(struct node_square_relation**head,struct node_square_relation**tail,int*num_B_smooth,int*num_semi_B_smooth){
    if(head==NULL || tail==NULL || num_B_smooth==NULL || num_semi_B_smooth==NULL){
        handle_error_with_exit("error in remove_same_num\n");
    }
    int removed=0;
    struct node_square_relation*l=*head;
    while(l!=NULL){
        while(l->next!=NULL && (mpz_cmp(l->square_relation.square,(l->next)->square_relation.square)==0) ) {
            if (mpz_cmp_si(l->square_relation.residuos, 1) == 0){//trovati due numeri uguali con residuo uguale a 1
                (*num_B_smooth)--;
                removed++;
            }
            else{//trovati due numeri uguali con residuo diverso da 1
                handle_error_with_exit("error in remove same square\n");
                (*num_semi_B_smooth)--;
            }
            remove_after_node_square_rel(&(l->next),tail);
        }
        if(l->next==NULL){
            return removed;
        }
        l=l->next;
    }
    return removed;
}
void union_list_square_v2(struct node_square_relation**head_square1,struct node_square_relation**tail_square1,struct node_square_relation**head_square2,struct node_square_relation**tail_square2){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    if(head_square1==NULL || tail_square1==NULL || head_square2==NULL || tail_square2==NULL){
        handle_error_with_exit("error in union_list_square_v2\n");
    }
    if(*head_square2==NULL){//nulla da fare
        return;
    }
    if(*tail_square1==NULL){//lista vuota,prendi l'altra lista
        *head_square1=*head_square2;
        *tail_square1=*tail_square2;
        return;
    }
    (*head_square2)->prev=(*tail_square1);//il primo nodo della seconda lista punta all'ultimo nodo della prima lista
    (*tail_square1)->next=*head_square2;//l'ultimo nodo della prima lista punta al primo nodo dell'altra lista
    *tail_square1=*tail_square2;//la coda della prima lista punta alla coda della seconda lista
    *tail_square2=NULL;
    *head_square2=NULL;
    return;
}
void union_list_residuos_v2(struct node_square_relation**head_residuos1,struct node_square_relation**tai_residuos1,struct node_square_relation**head_residuos2,struct node_square_relation**tail_residuos2){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    if(head_residuos1==NULL || tai_residuos1==NULL || head_residuos2==NULL || tail_residuos2==NULL){
        handle_error_with_exit("error in union_list_square_v2\n");
    }
    if(*head_residuos2==NULL){//nulla da fare
        return;
    }
    if(*tai_residuos1==NULL){//lista vuota,prendi l'altra lista
        *head_residuos1=*head_residuos2;
        *tai_residuos1=*tail_residuos2;
        return;
    }
    (*head_residuos2)->prev=(*tai_residuos1);//il primo nodo della seconda lista punta all'ultimo nodo della prima lista
    (*tai_residuos1)->next=*head_residuos2;//l'ultimo nodo della prima lista punta al primo nodo dell'altra lista
    *tai_residuos1=*tail_residuos2;//la coda della prima lista punta alla coda della seconda lista
    *tail_residuos2=NULL;
    *head_residuos2=NULL;
    return;
}

void union_list_square(struct node_square_relation**head_square,struct node_square_relation**tail_square,struct node_square_relation*phead_square,struct node_square_relation*ptail_square){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    #if DEBUG==1
    if(head_square==NULL || tail_square==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    #endif
    if(phead_square==NULL){//nulla da fare
        return;
    }
    if(*tail_square==NULL){//lista vuota,prendi l'altra lista
        *head_square=phead_square;
        *tail_square=ptail_square;
        return;
    }
    phead_square->prev=(*tail_square);//il primo nodo della seconda lista punta all'ultimo nodo della prima lista
    (*tail_square)->next=phead_square;//l'ultimo nodo della prima lista punta al primo nodo dell'altra lista
    *tail_square=ptail_square;//la coda della prima lista punta alla coda della seconda lista
    return;
}
void union_list_residuos(struct node_square_relation**head_residuos,struct node_square_relation**tail_residuos,struct node_square_relation*phead_residuos,struct node_square_relation*ptail_residuos){
    //concatena la prima lista e la seconda lista =L1 unito L2=L1,L2
    #if DEBUG==1
    if(head_residuos==NULL || tail_residuos==NULL){
        handle_error_with_exit("error in union_list_square\n");
    }
    #endif
    if(phead_residuos==NULL){
        return;
    }
    if(*tail_residuos==NULL){//lista vuota,prendi l'altra lista
        *head_residuos=phead_residuos;
        *tail_residuos=ptail_residuos;
        return;
    }
    phead_residuos->prev=(*tail_residuos);//il primo nodo della seconda lista punta all'ultimo nodo della prima lista
    (*tail_residuos)->next=phead_residuos;//l'ultimo nodo della prima lista punta al primo nodo dell'altra lista
    *tail_residuos=ptail_residuos;//la coda della prima lista punta alla coda della seconda lista
    return;
}
void add_relation_semi_B_smooth_to_list(struct node_square_relation**head_sort_residuos,struct node_square_relation**tail_sort_residuos,struct node_square_relation*head_residuos){
    if(head_sort_residuos==NULL || tail_sort_residuos==NULL){
        handle_error_with_exit("error in add_relation_semi_B_smooth_to_list\n");
    }
    if(head_residuos==NULL){//nessuan relazione da aggiungere
        return;
    }
    struct node_square_relation*p=head_residuos;
    while(p!=NULL){
        struct node_square_relation*q=p->next;
        insert_ordered_residuos_square_rel(p->square_relation,head_sort_residuos,tail_sort_residuos);
        free(p);
        p=q;
    }
    return;
}

void add_square_relation_to_list_sorted(struct node_square_relation**head_sort_square,struct node_square_relation**tail_sort_square,struct node_square_relation*head_square){
    if(head_sort_square==NULL || tail_sort_square==NULL){
        handle_error_with_exit("error in add_square_relation_to_list_sorted\n");
    }
    if(head_square==NULL){//nessuan relazione da aggiungere
        return;
    }
    struct node_square_relation*p=head_square;
    while(p!=NULL){
        struct node_square_relation*q=p->next;
        insert_ordered_sort_square_rel(p->square_relation,head_sort_square,tail_sort_square);
        free(p);
        p=q;
    }
    return;
}

void insert_ordered_residuos_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in insert_ordered\n");
    }
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    struct node_square_relation* temp = *tail;
    struct node_square_relation* next_node = NULL;
    struct node_square_relation* new_node = get_new_node_square_rel(square_relation);

    if(*head == NULL){
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    if(first_is_smaller_residuos_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_residuos_square_rel(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_square_rel(new_node, head, tail);
                return;
            }
        }
        next_node = temp->next;
        new_node->prev = temp;
        new_node->next = next_node;
        temp->next = new_node;
        next_node->prev = new_node;
    }
    return;
}
void insert_ordered_num_square_rel(struct square_relation square_relation, struct node_square_relation** head, struct node_square_relation** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
     if(square_relation.head_factorization==NULL){
        handle_error_with_exit("error in insert_ordered\n");
     }
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    struct node_square_relation* temp = *tail;
    struct node_square_relation* next_node = NULL;
    struct node_square_relation* new_node = get_new_node_square_rel(square_relation);

    if(*head == NULL){
        insert_first_square_rel(new_node, head, tail);
        return;
    }
    if(first_is_smaller_num_square_rel((**tail),*new_node)){
        insert_at_tail_square_rel(new_node,head, tail);
    }
    else{
        while(!first_is_smaller_num_square_rel(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_square_rel(new_node, head, tail);
                return;
            }
        }
        next_node = temp->next;
        new_node->prev = temp;
        new_node->next = next_node;
        temp->next = new_node;
        next_node->prev = new_node;
    }
    return;
}

void free_memory_list_square_rel(struct node_square_relation*head){
    if(head==NULL){
        return;
    }
    struct node_square_relation*p=head;
    struct node_square_relation*q;
    while(p!=NULL){
        q=p->next;
        free_memory_list_factor(p->square_relation.head_factorization);
        free(p);
        p=q;
    }
    return;
}
