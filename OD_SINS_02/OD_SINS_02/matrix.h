#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <math.h>

typedef struct  
{
	int size1;		//行数
	int size2;		//列数
	double *mp;		//矩阵首地址
} matrix;

matrix *matrix_calloc(int size1, int size2);
void matrix_free(matrix *mat);
void matrix_set(matrix *mat, int i, int j, double input);
double matrix_get(matrix *mat, int i, int j);
void matrix_transpose_memcpy(matrix *mat_t, matrix *mat);
void matrix_memcpy(matrix *s, matrix *d);
void matrix_add(matrix *a, matrix *b);
void matrix_sub(matrix *a, matrix *b);
void matrix_mul(matrix *a,matrix *b,matrix *c);
void matrix_mul(matrix *a,double b,matrix *c);
int matrix_inv(matrix *a);

#endif