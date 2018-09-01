/************************************************************************/
/* 矩阵基本运算                                                          */
/*     用于纯惯导解算和卡尔曼滤波的矩阵运算函数                            */
/* V2.0  2018-05-25  完成运算基本功能，加入简单的异常判断与抛出功能        */
/************************************************************************/
#include "matrix.h"

//异常号，用于判断矩阵计算的异常情况
int matrix_ErrorNum;

//==============================================
//	矩阵分配内存空间并初始化成0矩阵
//==============================================
matrix *matrix_calloc(int size1, int size2)
{
	double *mp = (double *)malloc(sizeof(double)*size1*size2);
	matrix *mat = (matrix *)malloc(sizeof(matrix));
	for (int k=0; k<size1*size2; k++)
	{
		mp[k] = 0.0;
	}
	mat->size1 = size1;
	mat->size2 = size2;
	mat->mp    = mp;

	return mat;
}

//==============================================
//	释放矩阵空间
//==============================================
void matrix_free(matrix *mat)
{
	free(mat->mp);
	free(mat);
}

//==============================================
//	矩阵元素赋值
//==============================================
void matrix_set(matrix *mat, int i, int j, double input)
{
	//判断给定元素是否在矩阵范围内
	if (i>mat->size1 || j>mat->size2)
	{
		matrix_ErrorNum = 2;
	}
	*(mat->mp+i*mat->size2 + j) = input;
}

//==============================================
//	矩阵元素取值
//==============================================
double matrix_get(matrix *mat, int i, int j)
{
	//判断元素是否在矩阵中
	return *(mat->mp+i*mat->size2 + j);
}

//==============================================
//	矩阵取列向量
//==============================================
void matrix_getv(matrix *mat,int j, matrix *v)
{
	//判断是否可以提取
	for (int i=0; i<mat->size1; i++)
	{
		matrix_set(v,i,0, matrix_get(mat,i,j));
	}
}

//==============================================
//	矩阵按行求和
//==============================================
void matrix_sumrow(matrix *mat, matrix *v)
{
	//判断是否可以提取
	double s;
	for (int i=0; i<mat->size1; i++)
	{
		s=0;
		for (int j=0; j<mat->size2; j++)
		{
			s += matrix_get(mat,i,j);
		}
		matrix_set(v,i,0,s);
	}
}

//=========================================================
//	矩阵转置并赋值
//=========================================================
void matrix_transpose_memcpy(matrix *mat_t, matrix *mat)
{
	//判断矩阵形状是否满足转置条件
	if (mat->size1!=mat_t->size2 || mat->size2!=mat_t->size1)
	{
		matrix_ErrorNum = 1;
	}
	for (int i=0; i<(mat->size1); i++)
	{
		for (int j=0; j<mat->size2; j++)
		{
			*(mat_t->mp+j*mat_t->size2 + i) = *(mat->mp+i*mat->size2 + j);
		}
	}
}

//==============================================
//	d=s
//==============================================
void matrix_memcpy(matrix *d, matrix *s)
{
	//判断矩阵形状是否相同
	if (d->size1!=s->size1 || d->size2!=s->size2)
	{
		matrix_ErrorNum = 3;
	}
	for (int i=0; i<s->size1; i++)
	{
		for (int j=0; j<s->size2; j++)
		{
			*(d->mp+i*d->size2 + j) = *(s->mp+i*s->size2 + j);
		}
	}
}

//==============================================
//	A=A+B
//==============================================
void matrix_add(matrix *a, matrix *b)
{
	//判断矩阵形状是否相同
	if (a->size1!=b->size1 || a->size2!=b->size2)
	{
		matrix_ErrorNum = 3;
	}
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<a->size2; j++)
		{
			*(a->mp+i*a->size2 + j) += *(b->mp+i*b->size2 + j);
		}
	}
}

//==============================================
//	C=A+B
//==============================================
void matrix_add(matrix *a, matrix *b, matrix *c)
{
	//判断矩阵形状是否相同
	if (a->size1!=b->size1 || a->size2!=b->size2)
	{
		matrix_ErrorNum = 3;
	}
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<a->size2; j++)
		{
			*(c->mp+i*c->size2 + j) = *(a->mp+i*a->size2 + j) + *(b->mp+i*b->size2 + j);
		}
	}
}

//==============================================
//	A=A-B
//==============================================
void matrix_sub(matrix *a, matrix *b)
{
	//判断矩阵形状是否相同
	if (a->size1!=b->size1 || a->size2!=b->size2)
	{
		matrix_ErrorNum = 3;
	}
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<a->size2; j++)
		{
			*(a->mp+i*a->size2 + j) -= *(b->mp+i*b->size2 + j);
		}
	}
}

//==============================================
//	C=A-B
//==============================================
void matrix_sub(matrix *a, matrix *b, matrix *c)
{
	//判断矩阵形状是否相同
	if (a->size1!=b->size1 || a->size2!=b->size2)
	{
		matrix_ErrorNum = 3;
	}
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<a->size2; j++)
		{
			*(c->mp+i*c->size2 + j) = *(a->mp+i*a->size2 + j) - *(b->mp+i*b->size2 + j);
		}
	}
}

//==============================================
//	C=A*B
//==============================================
void matrix_mul(matrix *a,matrix *b,matrix *c)
{
	//判断矩阵是否满足乘法
	if (a->size2!=b->size1)
	{
		matrix_ErrorNum = 4;
	}
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<b->size2; j++)
		{
			double sum=0.0;
			for (int k=0; k<b->size1; k++)
			{
				sum+=matrix_get(a,i,k)*matrix_get(b,k,j);
			}
			matrix_set(c,i,j,sum);
		}
	}
}

//==============================================
//	C=A×B (三维向量叉乘)
//==============================================
void matrix_crossmul(matrix *a,matrix *b,matrix *c)
{
	matrix *ax = matrix_calloc(3,3);
	matrix_set(ax, 0,1, -matrix_get(a, 2,0));
	matrix_set(ax, 0,2,  matrix_get(a, 1,0));
	matrix_set(ax, 1,0,  matrix_get(a, 2,0));
	matrix_set(ax, 1,2, -matrix_get(a, 0,0));
	matrix_set(ax, 2,0, -matrix_get(a, 1,0));
	matrix_set(ax, 2,1,  matrix_get(a, 0,0));

	matrix_mul(ax,b,c);

	matrix_free(ax);
}

//==============================================
//	C=A*b
//==============================================
void matrix_mul(matrix *a,double b,matrix *c)
{
	for (int i=0; i<a->size1; i++)
	{
		for (int j=0; j<a->size2; j++)
		{
			double res=0.0;
			res = matrix_get(a,i,j)*b;
			matrix_set(c,i,j,res);
		}
	}
}

//==============================================
//	norm = ||a||  (a为列向量)
//==============================================
double matrix_norm(matrix *a)
{
	int dim = a->size1;
	double norm = 0;
	for (int i=0; i<dim; i++)
	{
		norm += matrix_get(a,i,0)*matrix_get(a,i,0);
	}
	norm = sqrt(norm);
	return norm;
}

//==============================================
//	矩阵求逆
//==============================================
int matrix_inv(matrix *mat)
{
	//判断矩阵是否为方阵
	if (mat->size1!=mat->size2)
	{
		matrix_ErrorNum = 5;
	}
	double *a = mat->mp;
	int n = mat->size1;
	int *is = (int *)malloc(sizeof(double)*n);  
	int *js = (int *)malloc(sizeof(double)*n);  
	int i,j,k;  
	double d,p;  
	for ( k = 0; k < n; k++)  
	{  
		d = 0.0;  
		for (i=k; i<=n-1; i++)  
			for (j=k; j<=n-1; j++)  
			{  
				p=fabs(a[i*n+j]);  
				if (p>d) { d=p; is[k]=i; js[k]=j;}  
			}  
			if ( 0.0 == d )  
			{  
				free(is); free(js);  
				return -1;				//异常
			}  
			if (is[k]!=k)  
				for (j=0; j<=n-1; j++)  
				{  
					p=a[k*n+j];  
					a[k*n+j]=a[is[k]*n+j];  
					a[is[k]*n+j]=p;  
				}  
				if (js[k]!=k)  
					for (i=0; i<=n-1; i++)  
					{  
						p=a[i*n+k];  
						a[i*n+k]=a[i*n+js[k]];  
						a[i*n+js[k]]=p;  
					}  
					a[k*n+k] = 1.0/a[k*n+k];  
					for (j=0; j<=n-1; j++)  
						if (j!=k)  
						{  
							a[k*n+j] *= a[k*n+k];  
						}  
						for (i=0; i<=n-1; i++)  
							if (i!=k)  
								for (j=0; j<=n-1; j++)  
									if (j!=k)  
									{  
										a[i*n+j] -= a[i*n+k]*a[k*n+j];  
									}  
									for (i=0; i<=n-1; i++)  
										if (i!=k)  
										{  
											a[i*n+k] = -a[i*n+k]*a[k*n+k];  
										}  
	}  
	for ( k = n-1; k >= 0; k--)  
	{  
		if (js[k]!=k)  
			for (j=0; j<=n-1; j++)  
			{  
				p = a[k*n+j];  
				a[k*n+j] = a[js[k]*n+j];  
				a[js[k]*n+j]=p;  
			}  
			if (is[k]!=k)  
				for (i=0; i<=n-1; i++)  
				{   
					p = a[i*n+k];  
					a[i*n+k]=a[i*n+is[k]];  
					a[i*n+is[k]] = p;  
				}  
	}  
	free(is); free(js);   
	return 0;					//正常
}

//==============================================
//  矩阵计算异常处理（暂时不需要）
//==============================================
void matrix_errorprocess(int errornum)
{
	switch (errornum)
	{
	case 1:
		;
		break;
	case 2:
		;
		break;
	case 3:
		;
		break;
	default:
		;
		break;
	}
}

