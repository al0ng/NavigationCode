#include "myKalmanFilter.h"
#include "Navigation.h"

extern KF_Data kfData;
extern int matrix_ErrorNum;

void InitKalmanFilter()
{
	matrix_ErrorNum = 0;
	kfData.MeasureEnable = 0;
	kfData.MeasureUdtFinished = 0;
	kfData.ts = 0.01;
	kfData.Xk    = matrix_calloc(W, 1);
	kfData.Xkk_1 = matrix_calloc(W, 1);
	kfData.Yk    = matrix_calloc(V, 1);
	kfData.Pk    = matrix_calloc(W, W);
	kfData.Pkk_1 = matrix_calloc(W, W);
	kfData.Pkk_1 = matrix_calloc(W, W);
	kfData.Ft    = matrix_calloc(W, W);
	kfData.Phi   = matrix_calloc(W, W);
	kfData.Kk    = matrix_calloc(W, V);
	kfData.Hk    = matrix_calloc(V, W);
	kfData.G     = matrix_calloc(W, 6);
	kfData.Qk    = matrix_calloc(6, 6);
	kfData.Rk    = matrix_calloc(V, V);

	// 为矩阵赋初值
	matrix_set(kfData.Pk, 0,0, pow(5*Deg_Rad, 2));
	matrix_set(kfData.Pk, 1,1, pow(5*Deg_Rad, 2));
	matrix_set(kfData.Pk, 2,2, pow(5*Deg_Rad, 2));
	matrix_set(kfData.Pk, 3,3, pow(1.0, 2));
	matrix_set(kfData.Pk, 4,4, pow(1.0, 2));
	matrix_set(kfData.Pk, 5,5, pow(1.0, 2));
	matrix_set(kfData.Pk, 6,6, pow(10/Re, 2));
	matrix_set(kfData.Pk, 7,7, pow(10/Re, 2));
	matrix_set(kfData.Pk, 8,8, pow(10.0, 2));
	matrix_set(kfData.Pk, 9,9,   pow(40*Dph_Rps, 2));
	matrix_set(kfData.Pk, 10,10, pow(40*Dph_Rps, 2));
	matrix_set(kfData.Pk, 11,11, pow(40*Dph_Rps, 2));
	matrix_set(kfData.Pk, 12,12, pow(10*mg, 2));
	matrix_set(kfData.Pk, 13,13, pow(10*mg, 2));
	matrix_set(kfData.Pk, 14,14, pow(10*mg, 2));

	matrix_set(kfData.Hk, 0,3, 1);
	matrix_set(kfData.Hk, 1,4, 1);
	matrix_set(kfData.Hk, 2,5, 1);
	matrix_set(kfData.Hk, 3,6, 0);
	matrix_set(kfData.Hk, 4,7, 0);
	matrix_set(kfData.Hk, 5,8, 0);

	matrix_set(kfData.G, 0,0, 1);
	matrix_set(kfData.G, 1,1, 1);
	matrix_set(kfData.G, 2,2, 1);
	matrix_set(kfData.G, 3,3, 1);
	matrix_set(kfData.G, 4,4, 1);
	matrix_set(kfData.G, 5,5, 1);

	matrix_set(kfData.Qk, 0,0, pow(0.1*dp05h, 2));
	matrix_set(kfData.Qk, 1,1, pow(0.1*dp05h, 2));
	matrix_set(kfData.Qk, 2,2, pow(0.1*dp05h, 2));
	matrix_set(kfData.Qk, 3,3, pow(0.01*mg, 2));
	matrix_set(kfData.Qk, 4,4, pow(0.01*mg, 2));
	matrix_set(kfData.Qk, 5,5, pow(0.01*mg, 2));

	matrix_set(kfData.Rk, 0,0, pow(0.1, 2));
	matrix_set(kfData.Rk, 1,1, pow(0.1, 2));
	matrix_set(kfData.Rk, 2,2, pow(0.1, 2));
	matrix_set(kfData.Rk, 3,3, pow(10.0/Re, 2));
	matrix_set(kfData.Rk, 4,4, pow(10.0/Re, 2));
	matrix_set(kfData.Rk, 5,5, pow(10.0, 2));
}

void KalmanFilterSetPhi(double *imu_avp)
{
	double imu[6]={0.0}, avp[9]={0.0};
	double tanL, sinL, cosL, secL, Wz, Wn, Rs, Rc;
	double VnR, VnRR, VeR, VeRR, VeLR, VeLRR, Rmh, Rnh;
	double Fg[3], Vg[3];
	double Cbn[3][3] = {0.0};
	EARTH_DATA earth;
	int i, j;
	for (i=0; i<6; i++)
		imu[i] = imu_avp[i];
	for(i=0; i<9; i++)
		avp[i] = imu_avp[i+6];

	for(i=0; i<3; i++)
	{
		Fg[i] = 0;
		for(j=0; j<3; j++)
		{
			Fg[i] += Cbn[i][j] * imu[j+3] * TM_COUNT;
		}		
	}

	Convert_Att_Cbn(avp, Cbn[0]);
	Calculate_Earth(avp+6, avp+3, Cbn[0], &earth);

	Wn  = earth.Wn;
	Wz  = earth.Wz;
	Rs  = earth.Rs;
	Rc  = earth.Rc;
	Rmh = earth.Rmh;
	Rnh = earth.Rnh;
	Vg[0] = avp[3];
	Vg[1] = avp[4];
	Vg[2] = avp[5];
	tanL = tan(avp[6]);   //Lat：纬度
	sinL = sin(avp[6]);
	cosL = cos(avp[6]);
	secL = 1/cos(avp[6]);
	VnR  = avp[3]/Rmh;
	VnRR = avp[3]/(Rmh*Rmh);
	VeR  = avp[4]/Rnh;
	VeRR = avp[4]/(Rnh*Rnh);
	VeLR = avp[4]*tanL/Rnh;
	VeLRR= avp[4]*tanL/(Rnh*Rnh);

	matrix_set(kfData.Ft, 0,1, Wz + VeLR);   //wie*sinL + Ve*tanl/RNh
	matrix_set(kfData.Ft, 0,2, -(Wn + VeR)); //-(wie*cosL + Ve/RNh)
	matrix_set(kfData.Ft, 1,0, -(Wz + VeLR));  //-(wie*sinL + Ve*tanl/RNh)
	matrix_set(kfData.Ft, 1,2, -VnR);   //-Vn/RMh
	matrix_set(kfData.Ft,2,0, Wn + VeR);   //(wie*cosL + Ve/RNh)
	matrix_set(kfData.Ft,2,1, VnR);   //Vn/RMh

	matrix_set(kfData.Ft,0,4, -1./Rmh);   //-1/RMh
	matrix_set(kfData.Ft,1,3, 1./Rnh);  //1/RNh		 	
	matrix_set(kfData.Ft,2,3, tanL/Rnh); //tanL/RNh

	matrix_set(kfData.Ft, 0,8, VnRR);
	matrix_set(kfData.Ft, 1,6, Wie*sinL);
	matrix_set(kfData.Ft, 1,8, -VnRR);
	matrix_set(kfData.Ft, 2,6, Wie*cosL+VeR*secL*secL);
	matrix_set(kfData.Ft, 2,8, -VeLRR);
							 	
	matrix_set(kfData.Ft,0,9,  - Cbn[0][0]);
	matrix_set(kfData.Ft,0,10, - Cbn[0][1]);
	matrix_set(kfData.Ft,0,11, - Cbn[0][2]);
	matrix_set(kfData.Ft,1,9,  - Cbn[1][0]);
	matrix_set(kfData.Ft,1,10, - Cbn[1][1]);
	matrix_set(kfData.Ft,1,11, - Cbn[1][2]);
	matrix_set(kfData.Ft,2,9,  - Cbn[2][0]);
	matrix_set(kfData.Ft,2,10, - Cbn[2][1]);
	matrix_set(kfData.Ft,2,11, - Cbn[2][2]);
							 
	matrix_set(kfData.Ft,3,1, -Fg[2]);   //-fz_n
	matrix_set(kfData.Ft,3,2,  Fg[1]);   //fy_n
	matrix_set(kfData.Ft,4,0, Fg[2]);   //fz_n
	matrix_set(kfData.Ft,4,2, -Fg[0]);  //-fx_n		
	matrix_set(kfData.Ft,5,0, -Fg[1]);   //-fy_n
	matrix_set(kfData.Ft,5,1,  Fg[0]);   //fx_n
							
	matrix_set(kfData.Ft,3,3, (avp[4] * tanL - avp[5])/Rnh);  //(Vn*tanL - VU)/RNh
	matrix_set(kfData.Ft,3,4, 2 * Wz + VeLR);   //2wie*sinl + Ve*tanL/RNh
	matrix_set(kfData.Ft,3,5, -(2*Wn + VeR));   //-(2*wie*cosL+ Ve/Rnh)			 
	matrix_set(kfData.Ft,4,3, -2 * (Wz + VeLR));  //-2*(wie*sinL + Ve*tanL/Rnh)
	matrix_set(kfData.Ft,4,4, -avp[5]/Rmh);   //-VU/RMh
	matrix_set(kfData.Ft,4,5, -VnR);    //-Vn/RMh	 
	matrix_set(kfData.Ft,5,3,  2*(Wn + VeR));   //2*(wie*cosL + Ve/Rnh)
	matrix_set(kfData.Ft,5,4,  2*VnR);   //2*Vn/RMh

	matrix_set(kfData.Ft,3,8, Vg[0]*(Vg[2] - Vg[1]*tanL)/Rnh/Rnh);
	matrix_set(kfData.Ft,3,6, 2*Wie*(Vg[2]*Rs+Vg[1]*Rc) + VeR*Vg[1]/Rc/Rc);
	matrix_set(kfData.Ft,4,6, -Vg[0]*(2*Wn  + VeR/Rc/Rc));
	matrix_set(kfData.Ft,4,8, Vg[2]*VnR/Rmh + Vg[0]*VeLR/Rnh);
	matrix_set(kfData.Ft,5,6, -2*Vg[0]*Wz);
	matrix_set(kfData.Ft,5,8, -(VnR*VnR + VeR*VeR));
						
	matrix_set(kfData.Ft,3,12, Cbn[0][0]);
	matrix_set(kfData.Ft,3,13, Cbn[0][1]);
	matrix_set(kfData.Ft,3,14, Cbn[0][2]);
	matrix_set(kfData.Ft,4,12,  Cbn[1][0]);
	matrix_set(kfData.Ft,4,13,  Cbn[1][1]);
	matrix_set(kfData.Ft,4,14,  Cbn[1][2]);
	matrix_set(kfData.Ft,5,12, Cbn[2][0]);
	matrix_set(kfData.Ft,5,13, Cbn[2][1]);
	matrix_set(kfData.Ft,5,14, Cbn[2][2]);

	matrix_set(kfData.Ft,6,4, 1./Rmh);  //1/RMh
	matrix_set(kfData.Ft,7,3, 1./Rc/Rnh);  //secL/Rnh
	matrix_set(kfData.Ft,8,5, 1);
	matrix_set(kfData.Ft,6,8, -VnR/Rmh);   //-Vn/(RMh)^2
	matrix_set(kfData.Ft,7,6, VeLR/Rc);
	matrix_set(kfData.Ft,7,8, -VeR/Rc/Rnh);

	//Phi = I15 + Ft*ts;
	matrix *eye_15 = matrix_calloc(W,W);
	matrix *tmp = matrix_calloc(W,W);
	for(int ii=0; ii<W; ii++)
		matrix_set(eye_15, ii,ii, 1);
	matrix_mul(kfData.Ft, kfData.ts*2, tmp);
	matrix_add(tmp, eye_15);
	matrix_memcpy(kfData.Phi, tmp);
	matrix_free(tmp);
	matrix_free(eye_15);
}

//===========================================================
//	标准Kalman滤波更新
//===========================================================
void KalmanFilterStd()
{
	matrix *tmp1, *tmp11, *tmp12, *tmp2, *tmp21, *tmp22;
	matrix *G_t, *Phi_t, *H_t, *tmp_t, *Kk_t;
	matrix *eye_2 = matrix_calloc(W,W);
	for(int ii=0; ii<W; ii++)
		matrix_set(eye_2, ii,ii, 1);
	G_t = matrix_calloc(kfData.G->size2, kfData.G->size1);
	Phi_t = matrix_calloc(kfData.Phi->size2, kfData.Phi->size1);
	H_t = matrix_calloc(kfData.Hk->size2, kfData.Hk->size1);
	Kk_t = matrix_calloc(kfData.Kk->size2, kfData.Kk->size1);

	//l = length(Xk_1); 
	int matSize = kfData.Xk->size1;
	tmp_t = matrix_calloc(matSize, matSize);

	//Pkk_1 = Phi*Pk_1*Phi' + G*Qk*G';
	tmp1 = matrix_calloc(matSize, matSize);
	tmp11 = matrix_calloc(matSize, matSize);
	tmp2 = matrix_calloc(matSize, matSize);
	tmp21 = matrix_calloc(matSize, 6);
	matrix_mul(kfData.Phi, kfData.Pk, tmp11);
	matrix_transpose_memcpy(Phi_t, kfData.Phi);
	matrix_mul(tmp11, Phi_t, tmp1);
	matrix_mul(kfData.G, kfData.Qk, tmp21);
	matrix_transpose_memcpy(G_t, kfData.G);
	matrix_mul(tmp21, G_t, tmp2);
	matrix_add(tmp1, tmp2);
	matrix_memcpy(kfData.Pkk_1, tmp1);
	matrix_free(tmp1);
	matrix_free(tmp11);
	matrix_free(tmp2);
	matrix_free(tmp21);
	matrix_free(Phi_t);
	matrix_free(G_t);

	//Xkk_1 = Phi*Xk_1;
	tmp1 = matrix_calloc(matSize, 1);
	matrix_mul(kfData.Phi, kfData.Xk, tmp1);
	matrix_memcpy(kfData.Xkk_1, tmp1);
	matrix_free(tmp1);

	// 判断是否进行量测更新
	if (kfData.MeasureEnable == 0)
	{
		matrix_memcpy(kfData.Pk, kfData.Pkk_1);
		matrix_memcpy(kfData.Xk, kfData.Xkk_1);
		matrix_free(eye_2);
		matrix_free(tmp_t);
		matrix_free(H_t);
		matrix_free(Kk_t);
		return;
	}

	kfData.MeasureEnable = 0;		// 将量测标志位置零
	//Kk = Pkk_1*H'*((H*Pkk_1*H' + Rk)^-1);
	tmp1 = matrix_calloc(matSize, 6);
	tmp11 = matrix_calloc(matSize, matSize);
	tmp2 = matrix_calloc(6, 6);
	tmp21 = matrix_calloc(6, matSize);
	matrix_transpose_memcpy(H_t, kfData.Hk);
	matrix_mul(kfData.Pkk_1, H_t, tmp1);
	matrix_mul(kfData.Hk, kfData.Pkk_1, tmp21);
	matrix_mul(tmp21, H_t, tmp2);
	matrix_add(tmp2, kfData.Rk);
	matrix_inv(tmp2);
	matrix_mul(tmp1, tmp2, kfData.Kk);
	matrix_free(tmp1);
	matrix_free(tmp11);
	matrix_free(tmp2);
	matrix_free(tmp21);
	matrix_free(H_t);

	//Xk = Phi*Xk_1 + Kk*(Yk - H*Phi*Xk_1);
	tmp1 = matrix_calloc(matSize, 1);
	tmp11 = matrix_calloc(matSize, 1);
	tmp2 = matrix_calloc(6, 1);
	tmp21 = matrix_calloc(6, matSize);
	tmp22 = matrix_calloc(6, 1);
	matrix_mul(kfData.Phi, kfData.Xkk_1, tmp1);
	matrix_mul(kfData.Hk, kfData.Phi, tmp21);
	matrix_mul(tmp21, kfData.Xkk_1, tmp22);
	matrix_memcpy(tmp2, kfData.Yk);
	matrix_sub(tmp2, tmp22);
	matrix_mul(kfData.Kk, tmp2, tmp11);
	matrix_add(tmp1, tmp11);
	matrix_memcpy(kfData.Xk, tmp1);
	matrix_free(tmp1);
	matrix_free(tmp11);
	matrix_free(tmp2);
	matrix_free(tmp21);
	matrix_free(tmp22);

	//Pk = (eye(l) - Kk*H)*Pkk_1*(eye(l) - Kk*H)' + Kk*Rk*Kk';
	tmp1 = matrix_calloc(matSize, matSize);
	tmp11 = matrix_calloc(matSize, matSize);
	tmp12 = matrix_calloc(matSize, matSize);
	tmp2 = matrix_calloc(matSize, matSize);
	tmp21 = matrix_calloc(matSize, 6);
	tmp22 = matrix_calloc(1, 1);
	matrix_mul(kfData.Kk, kfData.Hk, tmp11);
	matrix_sub(eye_2, tmp11);
	matrix_transpose_memcpy(tmp_t, eye_2);
	matrix_mul(eye_2, kfData.Pkk_1, tmp12);
	matrix_mul(tmp12, tmp_t, tmp1);
	matrix_transpose_memcpy(Kk_t, kfData.Kk);
	matrix_mul(kfData.Kk, kfData.Rk, tmp21);
	matrix_mul(tmp21, Kk_t, tmp2);
	matrix_add(tmp1, tmp2);
	matrix_memcpy(kfData.Pk, tmp1);
	matrix_free(tmp1);
	matrix_free(tmp11);
	matrix_free(tmp12);
	matrix_free(tmp2);
	matrix_free(tmp21);
	matrix_free(tmp22);
	matrix_free(tmp_t);
	matrix_free(Kk_t);
	matrix_free(eye_2);

	kfData.MeasureUdtFinished = 1;
}

void KalmanFilterOut(FILE *fp_KFRes)
{
	static int i = 1;
	for (int j=0; j<15; j++)
		fprintf(fp_KFRes, "%.10lf ", matrix_get(kfData.Xk, j,0));
	for (int j=0; j<15; j++)
		fprintf(fp_KFRes, "%.10lf ", matrix_get(kfData.Pk, j,j));
	fprintf(fp_KFRes, "%4d\n", i);
	i++;
}