#include <stdio.h>
#include <string.h>
#include "myNavigation.h"
#include "matrix.h"

int main()
{
	FILE *fpa, *fpa1, *fpq, *fpm;
	fpa = fopen("att.txt","r");
	fpa1 = fopen("att1.txt","w");
	fpq = fopen("qnb.txt","w");
	fpm = fopen("cnb.txt","w");
	matrix* Cnb = matrix_calloc(3,3);
	matrix* Qnb = matrix_calloc(4,1);
	double  att[3] = {0.0};
	float  qnb_n[4] = {0.0};
	float  Cnb_n[3][3] = {0.0};
	for (int i=0; i<8; i++)
	{
		fscanf(fpa, "%lf %lf %lf\n", att,att+1,att+2);
// 		fscanf(fpq, "%f %f %f %f\n", qnb_n,qnb_n+1,qnb_n+2,qnb_n+3);
// 		fscanf(fpm, "%f %f %f\n%f %f %f\n%f %f %f\n", Cnb_n,Cnb_n+1,Cnb_n+2,Cnb_n+3,Cnb_n+4,Cnb_n+5,Cnb_n+6,Cnb_n+7,Cnb_n+8);
		a2qua(att, Qnb);
		q2mat(Qnb,Cnb);
		q2att(Qnb,att);
		m2qua(Cnb,Qnb);
		fprintf(fpa1, "%lf %lf %lf\n", att[0],att[1],att[2]);
		fprintf(fpq, "%.6f %.6f %.6f %.6f\n", matrix_get(Qnb,0,0), matrix_get(Qnb,1,0),
			    matrix_get(Qnb,2,0),matrix_get(Qnb,3,0));
		fprintf(fpm, "%.6f %.6f %.6f\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n",
			    matrix_get(Cnb,0,0),matrix_get(Cnb,0,1),matrix_get(Cnb,0,2),
				matrix_get(Cnb,1,0),matrix_get(Cnb,1,1),matrix_get(Cnb,1,2),
				matrix_get(Cnb,2,0),matrix_get(Cnb,2,1),matrix_get(Cnb,2,2));
	}

	fclose(fpa);
	fclose(fpa1);
	fclose(fpq);
	fclose(fpm);

	return 0;
}