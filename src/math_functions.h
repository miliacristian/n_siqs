#ifndef MATH_H
#define MATH_H
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "basic.h"

long calculate_mod_n(long a,long n);
void reduce_int_mod_n(int *a,int n);
int reduce_mod_2(int a);

int max(int i,int j);
int min(int i,int j);
void reduce_mod_n(long *a,long n);
void xgcd(long result [3], long x, long y);
long inverse_mod_n(long x,long n);
int reduce_int_mod_n_v2(int a,int n);
long power_mod_n(long base,long exponent,long n);
long gcd(long a,long b);
int rand_int(int a,int b);
int rand_long(long a,long b);
long floor_sqrt(long x);
#endif
