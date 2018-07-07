#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <mpfr.h>
int main(){
	system("make install");
	system("valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./criv_quad biprime2.txt");
	system("rm criv_quad");
	return 0;
}
