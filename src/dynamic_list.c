#include "dynamic_list.h"

int count_element_linked_list(struct node*head){
	int count=0;
	if (head==NULL){
		handle_error_with_exit("head is NULL\n");
	}
	struct node *p=head;
	while(p!=NULL){
		count++;
		p=p->next;
	}
	return count;
}
void get_element_linked_list(mpz_t elem,struct node*head,int index){//index start at 0
	if(index<0 || head==NULL || elem==NULL){
		handle_error_with_exit("index must be zero or positive\n");
	}
	struct node*p=head;
	int count=0;
	while(p!=NULL){
		if(count==index){
			mpz_set(elem,p->prime);
			return;
		}
		count++;
		p=p->next;
	}
	return;
}
void remove_after_node(struct node**ppos,struct node**tail){
	if(ppos==NULL || tail==NULL ){
		handle_error_with_exit("error in parameter remove_after_node\n");
	}
	 if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
  }
	struct node *r = *ppos;
	struct node*q=r->next;
	if(q==NULL){//fine della lista,bisogna aggiornare la coda
        *ppos=r->next;
        *tail=r->prev;
        free(r);
	r=NULL;
    }
    else {
        q->prev = r->prev;
        *ppos = r->next;
        free(r);
	r=NULL;
    }
    return;
}

int delete_head(struct node** head){//non è importante il valore iniziale di oldhead
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
	struct node*temp=(*head)->next;
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first(struct node *new_node, struct node **head, struct node **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail(struct node *new_node,struct node**head,struct node** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}


//alloca e inizializza un nodo della lista dinamica ordinata
struct node* get_new_node(mpz_t num) {
	if(num==NULL){
		handle_error_with_exit("error in parameter get_new_node\n");
	}
    struct node* new_node = (struct node*)malloc(sizeof(struct node));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    mpz_init(new_node->prime);
    mpz_set(new_node->prime,num);
    new_node->prev = NULL;
    new_node->next = NULL;
    return new_node;
}


void insert_at_head(struct node* new_node,struct node** head,struct node** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller(struct node node1, struct node node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
    if(mpz_cmp(node1.prime,node2.prime)>=0){//node1 è più grande di node2
        return 0;
    }
    return 1;//node1 è più piccolo di node 2
}

void insert_ordered(mpz_t num, struct node** head, struct node** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    struct node* temp = *tail;
    struct node* next_node = NULL;
    struct node* new_node = get_new_node(num);
    if(head==NULL || tail==NULL || num==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    if(*head == NULL){
        insert_first(new_node, head, tail);
        return;
    }
    if(first_is_smaller((**tail),*new_node)){
        insert_at_tail(new_node,head, tail);
    }
		else{
        while(!first_is_smaller(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head(new_node, head, tail);
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
void free_memory_list(struct node*head){
	if(head==NULL){
		return;
	}
	struct node*p=head;
	struct node*q;
	while(p!=NULL){
		q=p->next;
		mpz_clear(p->prime);
		free(p);
		p=q;
	}
	return;
}
