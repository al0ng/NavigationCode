#include "myKalmanFilter.h"
#include "myNavigation.h"

KF_Data kfData;
extern insdata ins;
extern int matrix_ErrorNum;
extern double imu_eb, imu_db, imu_web, imu_wdb;
extern double mea_dvn[3], mea_dpos[3], phi0[3], dpos0[3], dvn0;

void InitKalmanFilter()
{
	matrix_ErrorNum = 0;
	kfData.MeasureEnable = 0;
	kfData.MeasureUdtFinished = 0;
	kfData.ts = ins.nnts;				//Ӧ����ins�ṹ���е�nntsȷʵ
	kfData.Xk    = matrix_calloc(W, 1);
	kfData.Xkk_1 = matrix_calloc(W, 1);
	kfData.Yk    = matrix_calloc(V, 1);
	kfData.Pk    = matrix_calloc(W, W);
	kfData.Pkk_1 = matrix_calloc(W, W);
	kfData.Ft    = matrix_calloc(W, W);
	kfData.Phi   = matrix_calloc(W, W);
	kfData.Kk    = matrix_calloc(W, V);
	kfData.Hk    = matrix_calloc(V, W);
	kfData.G     = matrix_calloc(W, 6);
	kfData.Qk    = matrix_calloc(6, 6);
	kfData.Rk    = matrix_calloc(V, V);
	kfData.phi   = matrix_calloc(3,1);
	kfData.dvn   = matrix_calloc(3,1);
	kfData.dpos  = matrix_calloc(3,1);

	// Ϊ���󸳳�ֵ
	matrix_set(kfData.Pk, 0,0, pow(phi0[0]*3, 2));
	matrix_set(kfData.Pk, 1,1, pow(phi0[1]*3, 2));
	matrix_set(kfData.Pk, 2,2, pow(phi0[2]*3, 2));
	matrix_set(kfData.Pk, 3,3, pow(dvn0*3, 2));
	matrix_set(kfData.Pk, 4,4, pow(dvn0*3, 2));
	matrix_set(kfData.Pk, 5,5, pow(dvn0*3, 2));
	matrix_set(kfData.Pk, 6,6, pow(dpos0[0]*3/glv_Re, 2));
	matrix_set(kfData.Pk, 7,7, pow(dpos0[1]*3/glv_Re, 2));
	matrix_set(kfData.Pk, 8,8, pow(dpos0[2]*3, 2));
	matrix_set(kfData.Pk, 9,9,   pow(imu_eb*3, 2));
	matrix_set(kfData.Pk, 10,10, pow(imu_eb*3, 2));
	matrix_set(kfData.Pk, 11,11, pow(imu_eb*3, 2));
	matrix_set(kfData.Pk, 12,12, pow(imu_db*3, 2));
	matrix_set(kfData.Pk, 13,13, pow(imu_db*3, 2));
	matrix_set(kfData.Pk, 14,14, pow(imu_db*3, 2));
// 	matrix_set(kfData.Pk, 0,0, pow(5*glv_deg, 2));
// 	matrix_set(kfData.Pk, 1,1, pow(5*glv_deg, 2));
// 	matrix_set(kfData.Pk, 2,2, pow(5*glv_deg, 2));
// 	matrix_set(kfData.Pk, 3,3, pow(1.0, 2));
// 	matrix_set(kfData.Pk, 4,4, pow(1.0, 2));
// 	matrix_set(kfData.Pk, 5,5, pow(1.0, 2));
// 	matrix_set(kfData.Pk, 6,6, pow(10/glv_Re, 2));
// 	matrix_set(kfData.Pk, 7,7, pow(10/glv_Re, 2));
// 	matrix_set(kfData.Pk, 8,8, pow(10.0, 2));
// 	matrix_set(kfData.Pk, 9,9,   pow(150*glv_dph, 2));
// 	matrix_set(kfData.Pk, 10,10, pow(150*glv_dph, 2));
// 	matrix_set(kfData.Pk, 11,11, pow(150*glv_dph, 2));
// 	matrix_set(kfData.Pk, 12,12, pow(1*glv_mg, 2));
// 	matrix_set(kfData.Pk, 13,13, pow(1*glv_mg, 2));
// 	matrix_set(kfData.Pk, 14,14, pow(1*glv_mg, 2));

	matrix_set(kfData.Hk, 0,3, 1);
	matrix_set(kfData.Hk, 1,4, 1);
	matrix_set(kfData.Hk, 2,5, 1);
// 	matrix_set(kfData.Hk, 3,6, 1);
// 	matrix_set(kfData.Hk, 4,7, 1);
// 	matrix_set(kfData.Hk, 5,8, 1);

	matrix_set(kfData.G, 0,0, 1);
	matrix_set(kfData.G, 1,1, 1);
	matrix_set(kfData.G, 2,2, 1);
	matrix_set(kfData.G, 3,3, 1);
	matrix_set(kfData.G, 4,4, 1);
	matrix_set(kfData.G, 5,5, 1);

// 	matrix_set(kfData.Qk, 0,0, pow(0.001*dp05h, 2));
// 	matrix_set(kfData.Qk, 1,1, pow(0.001*dp05h, 2));
// 	matrix_set(kfData.Qk, 2,2, pow(0.001*dp05h, 2));
// 	matrix_set(kfData.Qk, 3,3, pow(0.001*mg, 2));
// 	matrix_set(kfData.Qk, 4,4, pow(0.001*mg, 2));
// 	matrix_set(kfData.Qk, 5,5, pow(0.001*mg, 2));
	matrix_set(kfData.Qk, 0,0, pow(imu_web, 2));
	matrix_set(kfData.Qk, 1,1, pow(imu_web, 2));
	matrix_set(kfData.Qk, 2,2, pow(imu_web, 2));
	matrix_set(kfData.Qk, 3,3, pow(imu_wdb, 2));
	matrix_set(kfData.Qk, 4,4, pow(imu_wdb, 2));
	matrix_set(kfData.Qk, 5,5, pow(imu_wdb, 2));

	matrix_set(kfData.Rk, 0,0, pow(mea_dvn[0], 2));
	matrix_set(kfData.Rk, 1,1, pow(mea_dvn[1], 2));
	matrix_set(kfData.Rk, 2,2, pow(mea_dvn[2], 2));
 	matrix_set(kfData.Rk, 3,3, pow(mea_dpos[0]/glv_Re, 2));
 	matrix_set(kfData.Rk, 4,4, pow(mea_dpos[1]/glv_Re, 2));
 	matrix_set(kfData.Rk, 5,5, pow(mea_dpos[2], 2));
}

void KalmanFilterSetPhi(double *imu_avp)
{
	double imu[6]={0.0}, avp[9]={0.0};
	double tanL, sinL, cosL, secL, Wz, Wn;
	double VnR, VnRR, VeR, VeRR, VeLR, VeLRR, Rmh, Rnh;
	double Vg[3];
	matrix *qnb = matrix_calloc(4,1);
	matrix *Cnb = matrix_calloc(3,3);
	matrix *fb  = matrix_calloc(3,1);
	matrix *fg  = matrix_calloc(3,1);
	extern EARTH earth;
	extern insdata ins;	//�ߵ��������ݽṹ��
	int i;
	for (i=0; i<6; i++)
		imu[i] = imu_avp[i];
	for(i=0; i<9; i++)
		avp[i] = imu_avp[i+6];
	for(i=0; i<3; i++)
		matrix_set(fb,i,0, imu[i+3]);
	a2qua(avp, qnb);
	q2mat(qnb, Cnb);
	qmulv(qnb, fb, fg);
	matrix_mul(fg, 1.0/ins.ts, fg);

	tanL = tan(avp[6]);   //Lat��γ��
	sinL = sin(avp[6]);
	cosL = cos(avp[6]);
	secL = 1/cos(avp[6]);
 	Wn  = glv_wie*cosL;
 	Wz  = glv_wie*sinL;
	Rmh = earth.RMh;
	Rnh = earth.RNh;
	Vg[0] = avp[3];
	Vg[1] = avp[4];
	Vg[2] = avp[5];
	VnR  = avp[4]/Rmh;
	VnRR = avp[4]/(Rmh*Rmh);
	VeR  = avp[3]/Rnh;
	VeRR = avp[3]/(Rnh*Rnh);
	VeLR = avp[3]*tanL/Rnh;
	VeLRR= avp[3]*tanL/(Rnh*Rnh);

	matrix_set(kfData.Ft, 0,1, Wz + VeLR);   //wie*sinL + Ve*tanl/RNh
	matrix_set(kfData.Ft, 0,2, -(Wn + VeR)); //-(wie*cosL + Ve/RNh)
	matrix_set(kfData.Ft, 1,0, -(Wz + VeLR));  //-(wie*sinL + Ve*tanl/RNh)
	matrix_set(kfData.Ft, 1,2, -VnR);   //-Vn/RMh
	matrix_set(kfData.Ft, 2,0, Wn + VeR);   //(wie*cosL + Ve/RNh)
	matrix_set(kfData.Ft, 2,1, VnR);   //Vn/RMh

	matrix_set(kfData.Ft, 0,4, -1./Rmh);   //-1/RMh
	matrix_set(kfData.Ft, 1,3, 1./Rnh);  //1/RNh		 	
	matrix_set(kfData.Ft, 2,3, tanL/Rnh); //tanL/RNh

	matrix_set(kfData.Ft, 0,8, VnRR);
	matrix_set(kfData.Ft, 1,6, glv_wie*sinL);
	matrix_set(kfData.Ft, 1,8, -VnRR);
	matrix_set(kfData.Ft, 2,6, glv_wie*cosL+VeR*secL*secL);
	matrix_set(kfData.Ft, 2,8, -VeLRR);
							 	
	matrix_set(kfData.Ft,0,9,  -matrix_get(Cnb,0,0));
	matrix_set(kfData.Ft,0,10, -matrix_get(Cnb,0,1));
	matrix_set(kfData.Ft,0,11, -matrix_get(Cnb,0,2));
	matrix_set(kfData.Ft,1,9,  -matrix_get(Cnb,1,0));
	matrix_set(kfData.Ft,1,10, -matrix_get(Cnb,1,1));
	matrix_set(kfData.Ft,1,11, -matrix_get(Cnb,1,2));
	matrix_set(kfData.Ft,2,9,  -matrix_get(Cnb,2,0));
	matrix_set(kfData.Ft,2,10, -matrix_get(Cnb,2,1));
	matrix_set(kfData.Ft,2,11, -matrix_get(Cnb,2,2));
							 
	matrix_set(kfData.Ft,3,1, -matrix_get(fg,2,0));   //-fz_n
	matrix_set(kfData.Ft,3,2,  matrix_get(fg,1,0));   //fy_n
	matrix_set(kfData.Ft,4,0,  matrix_get(fg,2,0));   //fz_n
	matrix_set(kfData.Ft,4,2, -matrix_get(fg,0,0));  //-fx_n		
	matrix_set(kfData.Ft,5,0, -matrix_get(fg,1,0));   //-fy_n
	matrix_set(kfData.Ft,5,1,  matrix_get(fg,0,0));   //fx_n
							
	matrix_set(kfData.Ft,3,3, (avp[4] * tanL - avp[5])/Rnh);  //(Vn*tanL - VU)/RNh
	matrix_set(kfData.Ft,3,4, 2 * Wz + VeLR);   //2wie*sinl + Ve*tanL/RNh
	matrix_set(kfData.Ft,3,5, -(2*Wn + VeR));   //-(2*wie*cosL+ Ve/Rnh)			 
	matrix_set(kfData.Ft,4,3, -2 * (Wz + VeLR));  //-2*(wie*sinL + Ve*tanL/Rnh)
	matrix_set(kfData.Ft,4,4, -avp[5]/Rmh);   //-VU/RMh
	matrix_set(kfData.Ft,4,5, -VnR);    //-Vn/RMh	 
	matrix_set(kfData.Ft,5,3,  2*(Wn + VeR));   //2*(wie*cosL + Ve/Rnh)
	matrix_set(kfData.Ft,5,4,  2*VnR);   //2*Vn/RMh

	matrix_set(kfData.Ft,3,8, Vg[0]*(Vg[2] - Vg[1]*tanL)/Rnh/Rnh);
	matrix_set(kfData.Ft,3,6, 2*glv_wie*(Vg[2]*sinL+Vg[1]*cosL) + VeR*Vg[1]/cosL/cosL);
	matrix_set(kfData.Ft,4,6, -Vg[0]*(2*Wn  + VeR/cosL/cosL));
	matrix_set(kfData.Ft,4,8, Vg[2]*VnR/Rmh + Vg[0]*VeLR/Rnh);
	matrix_set(kfData.Ft,5,6, -2*Vg[0]*Wz);
	matrix_set(kfData.Ft,5,8, -(VnR*VnR + VeR*VeR));
						
	matrix_set(kfData.Ft,3,12, matrix_get(Cnb,0,0));
	matrix_set(kfData.Ft,3,13, matrix_get(Cnb,0,1));
	matrix_set(kfData.Ft,3,14, matrix_get(Cnb,0,2));
	matrix_set(kfData.Ft,4,12, matrix_get(Cnb,1,0));
	matrix_set(kfData.Ft,4,13, matrix_get(Cnb,1,1));
	matrix_set(kfData.Ft,4,14, matrix_get(Cnb,1,2));
	matrix_set(kfData.Ft,5,12, matrix_get(Cnb,2,0));
	matrix_set(kfData.Ft,5,13, matrix_get(Cnb,2,1));
	matrix_set(kfData.Ft,5,14, matrix_get(Cnb,2,2));

	matrix_set(kfData.Ft,6,4, 1./Rmh);  //1/RMh
	matrix_set(kfData.Ft,7,3, 1./cosL/Rnh);  //secL/Rnh
	matrix_set(kfData.Ft,8,5, 1);

	matrix_set(kfData.Ft,6,8, -VnR/Rmh);   //-Vn/(RMh)^2
	matrix_set(kfData.Ft,7,6, VeLR/cosL);
	matrix_set(kfData.Ft,7,8, -VeR/cosL/Rnh);

	matrix_free(qnb);
	matrix_free(Cnb);
	matrix_free(fb);
	matrix_free(fg);

	//Phi = I15 + Ft*ts;
	matrix *eye_15 = matrix_calloc(W,W);
	matrix *tmp = matrix_calloc(W,W);
	for(int ii=0; ii<W; ii++)
		matrix_set(eye_15, ii,ii, 1);
	matrix_mul(kfData.Ft, kfData.ts, tmp);
	matrix_add(tmp, eye_15);
	matrix_memcpy(kfData.Phi, tmp);
	matrix_free(tmp);
	matrix_free(eye_15);
}

//===========================================================
//	��׼Kalman�˲�����
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

	// �ж��Ƿ�����������
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

	kfData.MeasureEnable = 0;		// �������־λ����
	//Kk = Pkk_1*H'*((H*Pkk_1*H' + Rk)^-1);
	tmp1 = matrix_calloc(matSize, V);
	tmp11 = matrix_calloc(matSize, matSize);
	tmp2 = matrix_calloc(V, V);
	tmp21 = matrix_calloc(V, matSize);
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

	//Xk = Xkk_1 + Kk*(Yk - H*Xkk_1);
	tmp1 = matrix_calloc(matSize, 1);
	tmp11 = matrix_calloc(matSize, 1);
	tmp2 = matrix_calloc(V, 1);
	tmp22 = matrix_calloc(V, 1);
	matrix_memcpy(kfData.Xkk_1, tmp1);
	matrix_mul(kfData.Hk, kfData.Xkk_1, tmp22);
	matrix_memcpy(tmp2, kfData.Yk);
	matrix_sub(tmp2, tmp22);
	matrix_mul(kfData.Kk, tmp2, tmp11);
	matrix_add(tmp1, tmp11);
	matrix_memcpy(kfData.Xk, tmp1);
	matrix_free(tmp1);
	matrix_free(tmp11);
	matrix_free(tmp2);
	matrix_free(tmp22);

	//Pk = (eye(l) - Kk*H)*Pkk_1*(eye(l) - Kk*H)' + Kk*Rk*Kk';
	tmp1 = matrix_calloc(matSize, matSize);
	tmp11 = matrix_calloc(matSize, matSize);
	tmp12 = matrix_calloc(matSize, matSize);
	tmp2 = matrix_calloc(matSize, matSize);
	tmp21 = matrix_calloc(matSize, V);
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

	//״̬�����������̬�ٶ�λ�����
	for (int i=0; i<3; i++)
	{
		matrix_set(kfData.phi,  i,0, matrix_get(kfData.Xk, i,0));
		matrix_set(kfData.dvn,  i,0, matrix_get(kfData.Xk, i+3,0));
		matrix_set(kfData.dpos, i,0, matrix_get(kfData.Xk, i+6,0));
	}
	kfData.MeasureUdtFinished = 1;
}

void KalmanFilterOut(FILE *fp_KFRes)
{
	static int i = 1;
	for (int j=0; j<15; j++)
		fprintf(fp_KFRes, "%.16lf ", matrix_get(kfData.Xk, j,0));
	for (int j=0; j<15; j++)
		fprintf(fp_KFRes, "%.16lf ", matrix_get(kfData.Pk, j,j));
	for (int j=0; j<6; j++)
		fprintf(fp_KFRes, "%.16lf ", matrix_get(kfData.Yk, j,0));
	fprintf(fp_KFRes, "%4d\n", i);
// 	for (int i=0; i<15; i++)
// 	{
// 		for (int j=0; j<15;j++)
// 		{
// 			fprintf(fp_KFRes, "%.16lf ", matrix_get(kfData.Phi, i,j));
// 		}
// 		fprintf(fp_KFRes, "\n");
// 	}
// 	fprintf(fp_KFRes, "\n");
	i++;
}