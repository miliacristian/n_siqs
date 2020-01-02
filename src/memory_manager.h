#ifndef MEM_MANAGER_H
#define MEM_MANAGER_H
#include <stdio.h>
#include <stdlib.h>

void* __wrap_malloc(size_t size);
void __wrap_free(void*data_to_free);
#endif