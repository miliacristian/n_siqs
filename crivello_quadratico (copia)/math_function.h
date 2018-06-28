#include "dynamic_list.h"
#include "print.h"
#include <gmp.h>
#include "list_factor_base.h"
#include "basic.h"
void square_root_mod_p_to_k(mpz_t rootpk,const mpz_t x,long p,const mpz_t n,int k);
void reduce_mod_n(long *a,long n);//reduce a mod n
void reduce_mod_n_int(int *a,int n);//reduce a mod n

long power_mod_n(long base,long exponent,long n);// base^exponent mod n

void xgcd(long result [3], long x, long y);//ax+by=gcd(x,y)
	//result[0]=gcd(x,y),result[1]=a,result[2]=b,l'array result deve essere passato vuoto

long gcd(long a,long b);

long inverse_mod_n(long x,long n);

long floor_sqrt(long x);  //cerca con ricerca binaria la radice intera di n

void select_square_mod_p(struct node**p,struct node** tail,const mpz_t n);

long tonelli_shanks(long n, long p);
char is_prime(const mpz_t num,struct node*head);

void create_num(mpz_t num,const mpz_t a,const mpz_t b,const mpz_t n,long j);
char divide_all_by_p_to_k(const mpz_t r,long p,int index_of_prime,long k,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,const mpz_t n,const mpz_t a,const mpz_t b);//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione
	
char divide_all_by_min1(mpz_t*array_of_number,long M,mpz_t**matrix_factorization);//r=square root mod p^k,ritorna zero se c'è stata almeno 1 divisione	
void factor_matrix(const mpz_t n,mpz_t*array_of_number,long M,mpz_t**matrix_factorization,struct node*head_f_base,int cardinality_factor_base,const mpz_t a,const mpz_t b);

int count_number_B_smooth(const mpz_t*array_number,long M,struct node**head_index_B_smooth,struct node**tail_index_B_smooth,const mpz_t a,const mpz_t b,struct node**head_square_B_smooth,struct node**tail_square_B_smooth);
long**save_element_B_smooth(const mpz_t*array_number,mpz_t**matrix_factorization,long M,int cardinality_f_base,int* num_of_B_smooth,const mpz_t a,const mpz_t b,struct node**head_square_B_smooth,struct node**tail_square_B_smooth);//alloca matrice con la fattorizzazione dei numeri b_smmoth
long calculate_mod_n(long a,long n);
int max(int i,int j);
int min(int i,int j);
void factor_matrix_f(const mpz_t n,long M,struct thread_data thread_data,int cardinality_factor_base,const mpz_t a,
                     struct a_struct*array_a_struct,int s);
int quadratic_residue(mpz_t x,const mpz_t q,const mpz_t n);
void add_exponent_of_a(long**matrix_B_smooth,int num_B_smooth,int s,const mpz_t a,long* array_of_prime_chosen_for_a);
void adjust_n(mpz_t n,int *k);
void reduce_array_mpz_mod_n(mpz_t*array,int length,const mpz_t a);
void set_to_odd_array_mpz_mod_n(mpz_t *array,int length,const mpz_t a);
mpz_t**save_element_B_smooth_in_matrix(mpz_t**matrix_factorization,long size,int cardinality_f_base,int* num_of_B_smooth);
int rand_int(int a,int b);
int rand_long(long a,long b);
void reduce_int_mod_n(int *a,int n);
int reduce_mod_2(int a);
int reduce_int_mod_n_v2(int a,int n);//reduce a mod n,rendendo a positivo maggiore o uguale a zero e minore di n,n è positivo
float calculate_log_thresold(const mpz_t n,long M);
void find_list_square_relation(struct thread_data thread_data, int *num_B_smooth,int*num_semi_B_smooth, int *num_potential_B_smooth, long M,
                               struct node_square_relation **head_square, struct node_square_relation **tail_square,struct node_square_relation **head_residuos, struct node_square_relation **tail_residuos,
                               const mpz_t n,const mpz_t a,struct a_struct*array_a_struct,int s);
void calculate_index_min_max_a(int*number_prime_a,int*index_prime_a,int length,int*min_a,int*max_a);
struct a_struct*create_array_a_struct(int*number_prime_a,int*index_number_a,int length);
void*thread_factorization_job(void*arg);
char divide_all_by_p_to_k_with_thread(int rad,long p,int index_of_prime,long k,long M,struct thread_data thread_data,const mpz_t n,const mpz_t a,const mpz_t b,struct a_struct*array_a_struct);











