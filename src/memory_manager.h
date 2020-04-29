#ifndef MEM_MANAGER_H
#define MEM_MANAGER_H
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "atomic.h"
#include <sys/mman.h>
extern void* __real_malloc(size_t size);//l'implementazione vera di malloc
extern void __real_free(void*data_to_free);//l'implementazione vera di malloc
void* __wrap_malloc(size_t size);
void __wrap_free(void*data_to_free);
#endif
