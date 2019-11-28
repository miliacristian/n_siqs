#include <pthread.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
void handle_error_with_exit(char*error_string);
FILE*open_file(char*path);