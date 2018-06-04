CC = gcc
CFLAGS =-Wall -ggdb3 -Wextra -o2
CFILES=$(shell ls *.c)
PROGS=$(CFILES:%.c=%)

install:criv_quad
	
criv_quad:criv_quad.c criv_quad.h math_function.c basic.c dynamic_list.c matrix_function.c print.c main.c miller_rabin.c miller_rabin.h list_factor_base.h list_factor_base.c
	$(CC) $(CFLAGS) -o  $@ $^  -lmpfr -lm -lgmp -lrt -pthread

client:Client.c Client.h basic.c basic.h file_lock.c file_lock.h io.c io.h dynamic_list.c dynamic_list.h parser.c  parser.h timer.c timer.h get_client.c get_client.h communication.h communication.c list_client.c list_client.h put_client.h put_client.c
	$(CC) $(CFLAGS) -pthread  -o  $@ $^ -lrt
