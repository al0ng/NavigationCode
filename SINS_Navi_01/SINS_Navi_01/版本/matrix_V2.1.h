#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <math.h>

//�������ͽṹ��
typedef struct  
{
	int size1;		//����
	int size2;		//����
	double *mp;		//�����׵�ַ
} matrix;	

typedef struct 
{
	matrix **stk;
	int sp;
	int len;
} matstk;

matstk* matmas_init(int stklen);
matrix* matmas_calloc(int size1, int size2, matstk* mstk);
int matmas_push(matrix* mp, matstk* mstk);
int matmas_pop(matstk *mstk);
void matmas_clear(matstk* mstk);

matrix *matrix_calloc(int size1, int size2);
void matrix_free(matrix *mat);
void matrix_set(matrix *mat, int i, int j, double input);
double matrix_get(matrix *mat, int i, int j);
void matrix_getv(matrix *mat,int j, matrix *v);
void matrix_sumrow(matrix *mat, matrix *v);
void matrix_transpose_memcpy(matrix *mat_t, matrix *mat);
void matrix_memcpy(matrix *d, matrix *s);
void matrix_add(matrix *a, matrix *b);
void matrix_add(matrix *a, matrix *b, matrix *c);
void matrix_sub(matrix *a, matrix *b);
void matrix_sub(matrix *a, matrix *b, matrix *c);
void matrix_crossmul(matrix *a,matrix *b,matrix *c);
void matrix_mul(matrix *a,matrix *b,matrix *c);
void matrix_mul(matrix *a,double b,matrix *c);
double matrix_norm(matrix *a);
int matrix_inv(matrix *a);

#endif