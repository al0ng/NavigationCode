#pragma once
#ifndef MYKALMANFILTER_H
#define MYKALMANFILTER_H

#include <time.h>
#include <stdio.h>
#include "Navigation.h"
#include "matrix.h"

#define W 15
#define V 3

typedef struct{
	matrix* Xk;
	matrix* Xkk_1;
	matrix* Yk;
	matrix* Pk;
	matrix* Pkk_1;
	matrix* Ft;
	matrix* Phi;
	matrix* Kk;
	matrix* Hk;
	matrix* G;
	matrix* Qk;
	matrix* Rk;
	int MeasureEnable;
	int MeasureUdtFinished;
	double ts;
} KF_Data;

void KalmanFilterStd();
void InitKalmanFilter();
void KalmanFilterSetPhi(double *imu_avp);
void KalmanFilterOut(FILE *fp_KFRes);

#endif