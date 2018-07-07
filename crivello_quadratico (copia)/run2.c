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
	system("./criv_quad biprime2.txt");
	system("rm criv_quad");
	return 0;
}
