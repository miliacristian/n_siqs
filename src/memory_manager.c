#include "memory_manager.h"

extern void* __real_malloc(size_t size);//l'implementazione vera di malloc
extern void __real_free(void*data_to_free);//l'implementazione vera di malloc
extern unsigned long num_times_malloc_called;
extern unsigned long num_times_free_called;
#define CACHE_LINE_SIZE_MIN_1 (CACHE_LINE_SIZE-1)
#if CUSTOM_MALLOC==1
#if LOCK_MEMORY==1
void __wrap_free(void*data_to_free){//malloc custom
	#if DEBUG==1
	atomic_inc_x86((atomic_t *)&num_times_free_called);
	#endif
	char *shift=data_to_free;
	shift--;
	shift=shift-(*shift)+1;
	__real_free(shift);
	return;
}

void* __wrap_malloc(size_t size){//malloc custom
	#if DEBUG==1
	atomic_inc_x86((atomic_t *)&num_times_malloc_called);
	if(CACHE_LINE_SIZE>=256){
		printf("char variables may overflow\n");
		return NULL;
	}
	#endif
	char* p= __real_malloc(size+CACHE_LINE_SIZE);
	int res=mlock(p,size+CACHE_LINE_SIZE);
	if(res!=0){
		handle_error_with_exit("impossible lock memory\n");
	}
	unsigned long address=(unsigned long)p;
	char remainder= address & CACHE_LINE_SIZE_MIN_1;
	if(remainder==0){
		address+=CACHE_LINE_SIZE_MIN_1;
		remainder=CACHE_LINE_SIZE;
	}
	else{
		address+=remainder-1;
	}
	p=(char*)address;
	*p=remainder;
	p++;
	return p;
	//do posix_memalign every malloc is very slow
}
#else
void __wrap_free(void*data_to_free){//malloc custom
	#if DEBUG==1
	atomic_inc_x86((atomic_t *)&num_times_free_called);
	#endif
	char *shift=data_to_free;
	shift--;
	shift=shift-(*shift)+1;
	__real_free(shift);
	return;
}

void* __wrap_malloc(size_t size){//malloc custom
	#if DEBUG==1
	atomic_inc_x86((atomic_t *)&num_times_malloc_called);
	if(CACHE_LINE_SIZE>=256){
		printf("char variables may overflow\n");
		return NULL;
	}
	#endif
	char* p= __real_malloc(size+CACHE_LINE_SIZE);
	unsigned long address=(unsigned long)p;
	char remainder= address & CACHE_LINE_SIZE_MIN_1;
	if(remainder==0){
		address+=CACHE_LINE_SIZE_MIN_1;
		remainder=CACHE_LINE_SIZE;
	}
	else{
		address+=remainder-1;
	}
	p=(char*)address;
	*p=remainder;
	p++;
	return p;
	//do posix_memalign every malloc is very slow
}
#endif
#endif
