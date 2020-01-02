#ifndef TIMING_H
#define TIMING_H
#include <stdio.h>
#include <time.h>
#include "parameters.h"
#include "basic.h"
void gettime(struct timespec*timer);

struct timespec diff_timespec(struct timespec time_current,struct timespec timer);
struct timespec print_time_elapsed(char*string);
void print_time_elapsed_local(char*string,struct timespec*timer_thread,unsigned int tid);
void print_total_time_elapsed(char*string,struct timespec timer_start);
#endif
