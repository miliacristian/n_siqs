CC = gcc
CFLAGS =-Wall -Wextra -o3 -finline-functions -funroll-loops
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

ifdef LINEAR_PINNING
CFLAGS:= $(CFLAGS) -DLINEAR_PINNING=$(LINEAR_PINNING)
else
CFLAGS:= $(CFLAGS) -DLINEAR_PINNING=0
endif

install:criv_quad
	
#factorize one number,print initialize,print time total,no check,print factorization,check factorization
criv_quad:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize one number,print initialize,print all timing,no check,print factorization,check_factorization
criv_quad_with_timimg:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize all numbers in file,print initialize,print all timing,print all things,all checks,print factorization,check_factorization
criv_quad_crash_test:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)

#factorize one number,print initialize,print all timing,print all things,all checks,print factorization,check_factorization,run with valgrind
criv_quad_debug:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^ $(LIBS)
	

criv_quad_clean:
	rm criv_quad

criv_quad_with_timimg_clean:
	rm criv_quad_with_timimg

criv_quad_crash_test_clean:
	rm criv_quad_crash_test

criv_quad_debug_clean:
	rm criv_quad_debug

clean:
	make criv_quad_clean