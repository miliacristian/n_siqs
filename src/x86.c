
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "atomic.h"

inline bool CAS_x86(volatile unsigned long long *ptr, unsigned long long oldVal, unsigned long long newVal) {
	unsigned long res = 0;

	__asm__ __volatile__(
		LOCK "cmpxchgq %1, %2;"//ZF = 1 if succeeded
		"lahf;"  // to get the correct result even if oldVal == 0
		"bt $14, %%ax;" // is ZF set? (ZF is the 6'th bit in %ah, so it's the 14'th in ax)
		"adc %0, %0" // get the result
		: "=r"(res)
		: "r"(newVal), "m"(*ptr), "a"(oldVal), "0"(res)
		: "memory"
	);

	return (bool)res;
}

inline bool iCAS_x86(volatile unsigned int *ptr, unsigned int oldVal, unsigned int newVal) {
	unsigned long res = 0;

	__asm__ __volatile__(
		LOCK "cmpxchgl %1, %2;" //ZF = 1 if succeeded
		"lahf;"  // to get the correct result even if oldVal == 0
		"bt $14, %%ax;" // is ZF set? (ZF is the 6'th bit in %ah, so it's the 14'th in ax)
		"adc %0, %0" // get the result
		: "=r"(res)
		: "r"(newVal), "m"(*ptr), "a"(oldVal), "0"(res)
		: "memory"
	);

	return (bool)res;
}

inline int atomic_test_and_set_x86(int *b) {
    int result = 0;

	__asm__  __volatile__ (
		LOCK "bts $0, %1;\n\t"
		"adc %0, %0"
		: "=r" (result)
		: "m" (*b), "0" (result)
		: "memory"
	);

	return !result;
}

inline int atomic_test_and_reset_x86(int *b) {
	int result = 0;

	__asm__  __volatile__ (
		LOCK "btr $0, %1;\n\t"
		"adc %0, %0"
		: "=r" (result)
		: "m" (*b), "0" (result)
		: "memory"
	);

	return result;
}

inline void atomic_add_x86(atomic_t *v, int i) {
	__asm__ __volatile__(
		LOCK "addl %1,%0"
		: "=m" (v->count)
		: "ir" (i), "m" (v->count)
	);
}

inline void atomic_sub_x86(atomic_t *v, int i) {
	__asm__ __volatile__(
		LOCK "subl %1,%0"
		: "=m" (v->count)
		: "ir" (i), "m" (v->count)
	);
}

inline void atomic_inc_x86(atomic_t *v) {
	__asm__ __volatile__(
		LOCK "incl %0"
		: "=m" (v->count)
		: "m" (v->count)
	);
}

inline void atomic_dec_x86(atomic_t *v) {
	__asm__ __volatile__(
		LOCK "decl %0"
		: "=m" (v->count)
		: "m" (v->count)
	);
}

inline int atomic_inc_and_test_x86(atomic_t *v) {
	unsigned char c = 0;

	__asm__ __volatile__(
		LOCK "incl %0\n\t"
		"sete %1"
		: "=m" (v->count), "=qm" (c)
		: "m" (v->count)
		: "memory"
	);
	return c != 0;
}


#ifdef SPINLOCK_GIVES_COUNT

inline unsigned int spin_lock(spinlock_t *s) {

        int count = -1;

        __asm__ __volatile__(
                "1:\n\t"
                "movl $1,%%eax\n\t"
                LOCK "xchgl %%eax, %2\n\t"
                "addl $1, %0\n\t"
                "testl %%eax, %%eax\n\t"
                "jnz 1b"
                : "=c" (count)
                : "c" (count), "m" (s->lock)
                : "eax", "memory"
        );

        return (unsigned int)count;
}


#else

inline void spin_lock_x86(spinlock_t *s) {

	__asm__ __volatile__(
		"1:\n\t"
		"movl $1,%%eax\n\t"
		LOCK "xchgl %%eax, %0\n\t"
		"testl %%eax, %%eax\n\t"
		"jnz 1b"
		: /* no output */
		: "m" (s->lock)
		: "eax", "memory"
	);
}
#endif

inline bool spin_trylock_x86(spinlock_t *s) {
/*	unsigned int out = 0;
	unsigned int in = 1;

	__asm__ __volatile__(
		"movl $1,%%eax\n\t"
		"xchgl %%eax, %0\n\t"
		"testl %%eax, %%eax\n\t"
		: "=a" ((unsigned int)out)
		: "m" (s->lock)
		: "eax", "memory"
	);

	__asm__ __volatile__(
		LOCK "xchgl %0, %1"
		:"=r" ((unsigned int)out)
		:"m" (s->lock), "0" (in)
		:"memory");

	return (bool)out;
*/
	return atomic_test_and_set_x86((int *)&s->lock);
}

inline void spin_unlock_x86(spinlock_t *s) {

	__asm__ __volatile__(
		"mov $0, %%eax\n\t"
		LOCK "xchgl %%eax, %0"
		: /* no output */
		: "m" (s->lock)
		: "eax", "memory"
	);
}
