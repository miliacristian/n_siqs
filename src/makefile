CC = gcc
CFLAGS =-Wall -Wextra -o3 -Xlinker --wrap=malloc -Xlinker --wrap=free
CFILES=$(shell ls *.c)
PROGS=$(CFILES:%.c=%)

ifdef DEBUG
CFLAGS:= $(CFLAGS) -DDEBUG=$(DEBUG)
else
CFLAGS:= $(CFLAGS) -DDEBUG=0
endif

install:criv_quad
	
#factorize one number,print initialize,print time total,no check,print factorization,check factorization
criv_quad:$(CFILES)
		$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread 

#factorize one number,print initialize,print all timing,no check,print factorization,check_factorization
criv_quad_with_timimg:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread

#factorize all numbers in file,print initialize,print all timing,print all things,all checks,print factorization,check_factorization
criv_quad_crash_test:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread

#factorize one number,print initialize,print all timing,print all things,all checks,print factorization,check_factorization,run with valgrind
criv_quad_debug:$(CFILES)
	$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread
	

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