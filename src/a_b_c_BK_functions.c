//
// Created by cristian on 25/11/19.
//

#include "a_b_c_BK_functions.h"
int compare_a_struct( const void* a, const void* b)
{
    struct a_struct s_a = * ( (struct a_struct*) a );
    struct a_struct s_b = * ( (struct a_struct*) b );
    if(s_a.index_prime_a==s_b.index_prime_a){
        return 0;
    }
    else if(s_a.index_prime_a>s_b.index_prime_a){
        return 1;
    }
    return -1;
}