#pragma once
#define WIN32_LEAN_AND_MEAN             //Exclude rarely used components from Windows headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

FILE* out;
int setupfiles(int argc, char const* argv[], FILE** in);
void put(const char* buffer);
void print(char* format, ...);
int valn(void* n);
int valacc(void* e);
void scan(const char* prompt, const char* format, void* x, int(*val)(void*));
int getarr(FILE* in, double* arr, int n);
int stdexit(const char* buffer, double* x, double* y, double* c, FILE* in);