CC = gcc
CFLAGS =-Wall -Wextra -funroll-loops -march=native -Ofast
CFILES=$(shell ls *.c)
PROGS=$(CFILES:%.c=%)

ifdef DEBUG
CFLAGS:= $(CFLAGS) -DDEBUG=$(DEBUG)
else
CFLAGS:= $(CFLAGS) -DDEBUG=0
endif

#wrap malloc and free
CFLAGS:= $(CFLAGS) -DCUSTOM_MALLOC=1	
CFLAGS:= $(CFLAGS) -Xlinker --wrap=malloc -Xlinker --wrap=free

CFLAGS:= $(CFLAGS) -DCACHE_LINE_SIZE=$(shell getconf LEVEL1_DCACHE_LINESIZE) #cache line size
CFLAGS:= $(CFLAGS) -DN_CPU=$(shell grep -c ^processor /proc/cpuinfo) #number of cpucores
CFLAGS:= $(CFLAGS) -DN_NUMA_NODES=$(shell numactl --hardware | grep -c cpus)

LIBS=-lmpfr -lm -lgmp -lrt -pthread

ifdef MAX_ALLOCABLE_GIGAS
CFLAGS:= $(CFLAGS) -DMAX_ALLOCABLE_GIGAS=$(MAX_ALLOCABLE_GIGAS)
else
CFLAGS:= $(CFLAGS) -DMAX_ALLOCABLE_GIGAS=4
endif

ifdef NUM_OF_N_TO_FACTORIZE
CFLAGS:= $(CFLAGS) -DNUM_OF_N_TO_FACTORIZE=$(NUM_OF_N_TO_FACTORIZE)
endif

ifdef LINEAR_PINNING
CFLAGS:= $(CFLAGS) -DLINEAR_PINNING=$(LINEAR_PINNING)
else
CFLAGS:= $(CFLAGS) -DLINEAR_PINNING=0
endif

install:criv_quad
	
#factorize one number,print initialize,print time total,no check(debug macro disabled),print factorization,check factorization
criv_quad:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize one number,print initialize,print all timing,no check(debug macro disabled),print factorization,check_factorization
criv_quad_with_timing:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize all numbers in file,print initialize,print all timing,print all things,all checks(debug macro enabled),print factorization,check_factorization
criv_quad_crash_test:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

criv_quad_debug:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize one number,print initialize,print all timing,print all things,all checks(debug macro enabled),print factorization,check_factorization,run with valgrind
criv_quad_valgrind:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)
	

criv_quad_clean:
	rm criv_quad

criv_quad_with_timing_clean:
	rm criv_quad_with_timing

criv_quad_crash_test_clean:
	rm criv_quad_crash_test

criv_quad_debug_clean:
	rm criv_quad_debug

criv_quad_valgrind_clean:
	rm criv_quad_valgrind

clean:
	make criv_quad_clean