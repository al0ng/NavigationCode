#ifndef MY_NAVIGATION_H
#define MY_NAVIGATION_H

#include "matrix.h"
#include <math.h>

#define PI			3.141592653
#define Rad_Deg		(180/PI)
#define glv_Re		6378137.0
#define glv_f		(1/298.257)
#define glv_e		0.0818192214555232
#define glv_e2		0.00669438499958795
#define glv_wie		7.2921151467e-5
#define glv_Rp		6356752.29821597
#define glv_g0		9.7803267714
#define glv_mg		0.00978032677140000
#define glv_ug		9.78032677140000e-06
#define glv_ppm		1.00000000000000e-06
#define glv_deg		0.0174532925199433
#define glv_min		0.000290888208665722
#define glv_sec		4.84813681109536e-06
#define glv_dps		0.0174532925199433
#define glv_dph		4.84813681109536e-06
#define glv_dpsh	0.000290888208665722

typedef struct{
	float ts;
	int nn;
	float nnts;
	matrix *dvm[4], *dwm[4];
	matrix *vmm, *wmm;
	matrix *qnb;
	matrix *vn;
	matrix *vn1;
	matrix *pos;
} insdata;

typedef struct{
	double att[3];
	double vn[3];
	double pos[3];
} avpdata;

typedef struct{
	double sl,cl,tl;
	double sl2, sl4;
	double g;
	matrix *gn, *gcc;
	double RMh, RNh, clRNh;
	matrix *wnie;
	matrix *wnen;
	matrix *wnin;
	matrix *wnien;
}EARTH;

int m2qua(matrix* m, matrix* q);
int q2att(matrix* qnb, double att[]);
int a2qua(double *att, matrix *qnb);
int rv2q(matrix* rv, matrix* q);
int q2mat(matrix* qnb, matrix* Cnb);
int qmul(matrix* q1, matrix* q2, matrix* q);
int qmulv(matrix *qnb, matrix *vi, matrix *vo);
int qnormlz(matrix *qnb);

void earthupdate();
int insinitial(double avp0[9], int nn);

int attupdate();
int vnupdate();
int posupdate();
void insupdate(matrix *wm, matrix *vm);
void ins2avp();

#endif
