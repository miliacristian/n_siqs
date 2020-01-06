#ifndef BASIC_H
#define BASIC_H

#include <pthread.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

void handle_error_with_exit(char*error_string);
void destroy_mtx(pthread_mutex_t *mtx);
void initialize_mtx(pthread_mutex_t *mtx);
void lock_mtx(pthread_mutex_t *mtx);
void unlock_mtx(pthread_mutex_t *mtx);
FILE*open_file_log();
long get_file_size(char*path);
FILE*open_file(char*path);
int pin_thread_to_core(int core,cpu_set_t*oldset);
#endif
