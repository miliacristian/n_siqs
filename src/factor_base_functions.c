#include "factor_base_functions.h"
#include "timing.h"

extern long B;
extern int cardinality_factor_base;
extern struct timespec timer_test;
extern double thresold_relation;

#if DEBUG==1
char verify_sorted_list(struct node_factor_base*head,int length){
    if(head==NULL || length<=0){
        handle_error_with_exit("error in verify sorted list\n");
    }
    int count_elem=0;
    struct node_factor_base*p=head;
    while(p!=NULL) {
        if (p->prime > p->next->prime) {
            return 0;
        }
        count_elem++;//aumenta numero elementi
        p = p->next;
    }
    if(count_elem!=length){
       return 0;
    }
    return 1;
}
char verify_factor_base(struct node_factor_base*head,int cardinality_factor_base,int last_prime_factor_base){
    if(head==NULL || cardinality_factor_base<=0 || last_prime_factor_base<=0){
      handle_error_with_exit("error in verify factor base");
    }
    int count_elem=0;
    struct node_factor_base*p=head;
    while(p!=NULL){
        if(p->next==NULL){//ultimo elemento
            if(p->prime!=last_prime_factor_base){
                return 0;
            }
            count_elem++;
            p=p->next;
            continue;
        }
        if(count_elem==0 && p->prime!=-1){
            return 0;
        }
        if(count_elem==1 && p->prime!=2){
            return 0;
        }
        if(p->prime>=p->next->prime){
            return 0;
        }
        count_elem++;//aumenta numero elementi
        p=p->next;
    }
    if(count_elem!=cardinality_factor_base){//se il numero di elementi è diverso dalla cardinalità errore
        return 0;
    }
    return 1;
}
char verify_cardinality_list_factor_base(struct node_factor_base*head,int length){
    if(length<0){
        printf("length minore di zero\n");
        return 0;
    }
    int counter=0;
    if(head==NULL && length==0){
        return 1;
    }
    struct node_factor_base*p=head;
    while(p!=NULL){
        counter++;
        p=p->next;
    }
    if(counter!=length){
        printf("lunghezza errata\n");
        return 0;
    }
    return 1;
}
#endif

void print_list_factor_base(struct node_factor_base*head,int length){
	if (head==NULL || length<0){
		printf("impossible print list head is NULL\n");
		return;
	}
	if(length==0){
		printf("list is empty\n");
	}
	if(not_print_list(length,THRESOLD_PRINT_LIST)==1){
		return;
	}
	struct node_factor_base *p=head;
	while(p!=NULL){
		printf("%d sq=%d,",p->prime,p->root_n_mod_prime);
		p=p->next;
	}
	printf("\n");
	return;
}

void print_factor_base(struct node*head_f_base){
	if(head_f_base==NULL){
		printf("factor_base is empty\n");
		return;
	}
	printf("factor_base=");
	print_list_mpz(head_f_base);
	return;
}

int count_element_linked_list_f(struct node_factor_base*head){
	int count=0;
	if (head==NULL){
		handle_error_with_exit("head is NULL\n");
	}
	struct node_factor_base *p=head;
	while(p!=NULL){
		count++;
		p=p->next;
	}
	return count;
}
void get_element_linked_list_f(int *elem,struct node_factor_base*head,int index){//index start at 0
	if(index<0 || head==NULL || elem==NULL){
		handle_error_with_exit("index must be zero or positive\n");
	}
	struct node_factor_base*p=head;
	int count=0;
	while(p!=NULL){
		if(count==index){
			*elem=p->prime;
			return;
		}
		count++;
		p=p->next;
	}
	return;
}
void remove_after_node_f(struct node_factor_base**ppos,struct node_factor_base**tail){
	if(ppos==NULL || tail==NULL ){
		handle_error_with_exit("error in parameter remove_after_node\n");
	}
	 if(*ppos==NULL || *tail==NULL){
        handle_error_with_exit("impossible remove node,empty list\n");
  }
	struct node_factor_base *r = *ppos;
	struct node_factor_base*q=r->next;
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

int delete_head_f(struct node_factor_base** head){//non è importante il valore iniziale di oldhead
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
	struct node_factor_base*temp=(*head)->next;
        free(*head);
        *head =temp;
        (*head)-> prev = NULL;
    }
    return 0;
}

void insert_first_f(struct node_factor_base *new_node, struct node_factor_base **head, struct node_factor_base **tail){//inserisce il primo nodo
    if(new_node==NULL || head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_first\n");
    }
    *head = new_node;
    *tail = new_node;
    return;
}

void insert_at_tail_f(struct node_factor_base *new_node,struct node_factor_base**head,struct node_factor_base** tail){//inserisce un nodo in coda
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_f(new_node, head, tail);
        return;
    }
    (*tail)->next = new_node;
    new_node->prev = *tail;
    *tail = new_node;
}
void union_list_factor_base(struct node_factor_base**head1,struct node_factor_base**tail1,int*cardinality_factor_base1,int*last_prime_factor_base1,
        struct node_factor_base*head2,struct node_factor_base*tail2,int cardinality_factor_base2,int last_prime_factor_base2){
    //concatena la prima lista e la seconda lista ----> L1 unito L2=L1,L2 suppone che la lista 1 contenga primi minori della lista 2
    if(head1==NULL || tail1==NULL || cardinality_factor_base1==NULL || (*cardinality_factor_base1)<0 || (*last_prime_factor_base1)<=0
            || cardinality_factor_base2<0 || last_prime_factor_base2<=0){
        handle_error_with_exit("error 1\n");
    }
    if(head2==NULL){
        return;
    }
    if(*tail1==NULL){//lista vuota,prendi l'altra lista
        *head1=head2;
        *tail1=tail2;
        *cardinality_factor_base1=cardinality_factor_base2;
        *last_prime_factor_base1=last_prime_factor_base2;
        return;
    }
    head2->prev=(*tail1);//il primo nodo della seconda lista punta al nodo ultimo della lista
    (*tail1)->next=head2;//l'ultimo nodo punta al primo nodo dell'altra lista
    *tail1=tail2;//la coda punta alla coda dell'altra lista
    *cardinality_factor_base1=*cardinality_factor_base1+cardinality_factor_base2;
    *last_prime_factor_base1=last_prime_factor_base2;
    return;
}

//alloca e inizializza un nodo della lista dinamica ordinata
struct node_factor_base* get_new_node_f(int num,const mpz_t n) {
	if(num<=0 && num!=-1){
		handle_error_with_exit("error in get_new_node\n");
	}
	int t;
	mpz_t n_temp,p_temp,root;
    struct node_factor_base* new_node = (struct node_factor_base*)malloc(sizeof(struct node_factor_base));
    if(new_node==NULL){
        handle_error_with_exit("error in malloc get_new_node\n");
    }
    mpz_init(n_temp);
    mpz_init(p_temp);
    mpz_init(root);

    new_node->prime=num;
    new_node->prev = NULL;
    new_node->next = NULL;

    if(num!=-1 && num!=2) {
        mpz_set(n_temp, n);//n_temp=n
        mpz_set_si(p_temp, num);//p_temp=p
        mpz_mod(n_temp, n_temp, p_temp);//n_temp = n mod p
        t = quadratic_residue(root, n_temp, p_temp);//r1=radice quadrata di n mod p
        if (t == -1 || t == 0) {
            handle_error_with_exit("error in calculate quadratic_residue\n");
        }
        new_node->root_n_mod_prime=mpz_get_si(root);
    }
    else{
        new_node->root_n_mod_prime=0;
    }

    mpz_clear(n_temp);
    mpz_clear(p_temp);
    mpz_clear(root);
    return new_node;
}


void insert_at_head_f(struct node_factor_base* new_node,struct node_factor_base** head,struct node_factor_base** tail) {//inserisce un nodo in testa alla lista
    if(head==NULL){
        handle_error_with_exit("error in insert_at_head **head is NULL\n");
    }
    if(tail==NULL){
        handle_error_with_exit("error in insert_at_head **tail is NULL\n");
    }
    if(*head == NULL) {
        insert_first_f(new_node, head, tail);
        return;
    }
    (*head)->prev = new_node;
    new_node->next = *head;
    *head = new_node;
    return;
}


char first_is_smaller_f(struct node_factor_base node1, struct node_factor_base node2){//verifica se il primo nodo contiene tempi più piccoli del secondo nodo
	if(node1.prime>=node2.prime){
		return 0;
	}
    return 1;//node1 è più piccolo di node 2
}

void insert_ordered_f(int num,const mpz_t n, struct node_factor_base** head, struct node_factor_base** tail){
    //inserisce ordinatamente un nodo nella lista ordinata per istanti temporali
    struct node_factor_base* temp = *tail;
    struct node_factor_base* next_node = NULL;
    struct node_factor_base* new_node = get_new_node_f(num,n);
    if(head==NULL || tail==NULL){
        handle_error_with_exit("error in insert_ordered,head or tail are NULL\n");
    }
    if(*head == NULL){
        insert_first_f(new_node, head, tail);
        return;
    }
    if(first_is_smaller_f((**tail),*new_node)){
        insert_at_tail_f(new_node,head, tail);
    }
		else{
        while(!first_is_smaller_f(*temp,*new_node)){
            if(temp->prev != NULL){
                temp = temp->prev;
            }else{
                insert_at_head_f(new_node, head, tail);
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
void free_memory_list_f(struct node_factor_base*head){
	if(head==NULL){
		return;
	}
	struct node_factor_base*p=head;
	struct node_factor_base*q;
	while(p!=NULL){
		q=p->next;
		free(p);
		p=q;
	}
	return;
}

void print_estimated_time(int cardinality_factor_base,int num_B_smooth){
    if(cardinality_factor_base<=0 || num_B_smooth<0){
        handle_error_with_exit("error in print_estimated_time");
    }
    long ns,ms,sec,min,hour,temp;
    struct timespec empty_struct,time_sub;
    empty_struct.tv_nsec=0;
    empty_struct.tv_sec=0;
    time_sub=diff_timespec(timer_test,empty_struct);
    ns=time_sub.tv_nsec%1000000;
    time_sub.tv_nsec-=ns;
    ms=time_sub.tv_nsec/1000000;
    sec=time_sub.tv_sec%60;//i secondi sono modulo 60
    min=(time_sub.tv_sec-sec)/60;
    temp=min%60;//temp min
    hour=(min-temp)/60;
    min=temp;

    long times=((double)cardinality_factor_base*thresold_relation)/(double)num_B_smooth;
    ms=ms*times;
    long remainder=ms%1000;
    int plus_sec=(ms-remainder)/1000;//secondi in più
    ms=remainder;
    sec=sec*times;
    remainder=sec%60;
    int plus_min=(sec-remainder)/60;
    sec=remainder;
    min+=plus_min;
    sec+=plus_sec;
    remainder=sec%60;
    plus_min=(sec-remainder)/60;
    sec=remainder;
    min+=plus_min;
    printf("hour=%ld min=%ld sec:%ld ms=%ld ns:%ld\n",hour,min,sec,ms,ns);
    return;
}

struct factor_base_data* alloc_array_factor_base_data(int length){
    if(length<=0){
        handle_error_with_exit("error in alloc_array_factor_base_data\n");
    }
    struct factor_base_data*array_factor_base=malloc(sizeof(struct factor_base_data)*length);
    if(array_factor_base==NULL){
        handle_error_with_exit("error in malloc array_factor_base\n");
    }
    memset(array_factor_base,0, sizeof(struct factor_base_data)*length);
    return array_factor_base;
}

struct node_factor_base* initialize_factor_base(int*cardinality_factor_base,long B,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base){
    if(B<2 || tail==NULL || cardinality_factor_base==NULL || n==NULL || last_prime_factor_base==NULL){
        handle_error_with_exit("error in parameter\n");
    }
    struct node_factor_base *head=NULL;
    mpz_t p;
    mpz_init(p);
    mpz_set_si(p,-1);//temp=-1
    insert_ordered_f(-1,n,&head,tail);//inserisci -1 nella factor base
    (*cardinality_factor_base)++;
    mpz_set_si(p,2);//p=2
    (*cardinality_factor_base)++;
    insert_ordered_f(2,n,&head,tail);//inserisce 2 nella factor base
    *last_prime_factor_base=2;
    mpz_clear(p);
    return head;
}

void create_factor_base_f(int*cardinality_factor_base,long B,struct node_factor_base**head,struct node_factor_base**tail,const mpz_t n,int *last_prime_factor_base){//crea la factor base aggiungendo i numeri primi minori o uguali a B
    if(B<2 || tail==NULL || head==NULL || cardinality_factor_base==NULL || n==NULL || last_prime_factor_base==NULL ){
        handle_error_with_exit("error in parameter\n");
    }
    struct node_factor_base*node;
    long v,v_temp,m;
    mpz_t p,pp,value;


    mpz_init(p);
    mpz_init(pp);
    mpz_init(value);
    long start;
    if(*last_prime_factor_base==1) {
        *head=initialize_factor_base(cardinality_factor_base, B,tail, n,last_prime_factor_base);
        //initialize factor base pone last_prime_factor_base a 1
    }
    if((*last_prime_factor_base & 1)==0){//start=last_prime_factor_base è pari==2
        start=*last_prime_factor_base;
    }
    else{//last_prime_factor_base è dispari
        start=*last_prime_factor_base+1;//start è pari
    }
    for(long i=start;i<=B;){
        mpz_set_si(p,i);//p=i,inizialmente p=2,ad ogni inizio ciclo deve essere pari
        mpz_nextprime(p,p);//p=next_prime,ritorna il prossimo primo maggiore strettamente di p,è bene che p sia pari
        if(mpz_cmp_si(p,B)<=0){//primo<=B
            v=mpz_get_si(p);
            v_temp=v;
            v=(v-1)>>1;//shift a destra divide per 2,v=(v-1/2)
            mpz_powm_ui(value,n,(unsigned long int)v,p);//value=n^v mod p
            m=mpz_get_si(value);//m=n^((p-1)/2) mod p
            if(m==1){//n è quadrato modulo p
                node=get_new_node_f(mpz_get_si(p),n);
                insert_at_tail_f(node,head,tail);
                (*cardinality_factor_base)++;
                *last_prime_factor_base=v_temp;
            }
            i=mpz_get_si(p)+1;//il numero diventa pari e il prossimo numero primo sarà dispari
        }
        else{//B è stato superato
            break;
        }
    }
    mpz_clear(p);
    mpz_clear(pp);
    mpz_clear(value);
    return;
}