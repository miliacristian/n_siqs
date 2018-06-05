CC = gcc
CFLAGS =-Wall -ggdb3 -Wextra -o2
CFILES=$(shell ls *.c)
PROGS=$(CFILES:%.c=%)

install:criv_quad
	
criv_quad:criv_quad.c criv_quad.h math_function.c math_function.h basic.c basic.h dynamic_list.c dynamic_list.h matrix_function.c matrix_function.h print.c print.h main.c miller_rabin.c miller_rabin.h list_factor_base.h list_factor_base.c list_factorization.c list_factorization.h list_square_relation.c list_square_relation.h
	$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread

client:Client.c Client.h basic.c basic.h file_lock.c file_lock.h io.c io.h dynamic_list.c dynamic_list.h parser.c  parser.h timer.c timer.h get_client.c get_client.h communication.h communication.c list_client.c list_client.h put_client.h put_client.c
	$(CC) $(CFLAGS) -pthread  -o  $@ $^ -lrt
