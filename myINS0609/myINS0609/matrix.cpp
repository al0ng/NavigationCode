/************************************************************************/
/* �����������                                                          */
/*     ���ڿ������˲��ľ������㺯��                                       */
/* V1.0  2018-04-20  �������������ܣ�û���쳣�ж����׳�����              */
/************************************************************************/
#include "matrix.h"

int matrix_ErrorNum;

//==============================================
//	��������ڴ�ռ䲢��ʼ����0����
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
//	�ͷž���ռ�
//==============================================
void matrix_free(matrix *mat)
{
	free(mat->mp);
	free(mat);
}

//==============================================
//	����Ԫ�ظ�ֵ
//==============================================
void matrix_set(matrix *mat, int i, int j, double input)
{
	//�жϸ���Ԫ���Ƿ��ھ���Χ��
	if (i>mat->size1 || j>mat->size2)
	{
		matrix_ErrorNum = 2;
	}
	*(mat->mp+i*mat->size2 + j) = input;
}

//==============================================
//	����Ԫ����ֵ
//==============================================
double matrix_get(matrix *mat, int i, int j)
{
	//�ж�Ԫ���Ƿ��ھ�����
	return *(mat->mp+i*mat->size2 + j);
}

//=========================================================
//	����ת�ò���ֵ
//=========================================================
void matrix_transpose_memcpy(matrix *mat_t, matrix *mat)
{
	//�жϾ�����״�Ƿ�����ת������
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
	//�жϾ�����״�Ƿ���ͬ
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
	//�жϾ�����״�Ƿ���ͬ
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
//	A=A-B
//==============================================
void matrix_sub(matrix *a, matrix *b)
{
	//�жϾ�����״�Ƿ���ͬ
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
//	C=A*B
//==============================================
void matrix_mul(matrix *a,matrix *b,matrix *c)
{
	//�жϾ����Ƿ�����˷�
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
//	��������
//==============================================
int matrix_inv(matrix *mat)
{
	//�жϾ����Ƿ�Ϊ����
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
				return -1;				//�쳣
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
	return 0;					//����
}

//��������쳣������ʱ����Ҫ��
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

