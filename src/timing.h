#include <stdio.h>
#include <time.h>
#include "parameters.h"
#include "basic.h"
void gettime(struct timespec*timer);

struct timespec diff_timespec(struct timespec time_current,struct timespec timer);
struct timespec print_time_elapsed(char*string);
