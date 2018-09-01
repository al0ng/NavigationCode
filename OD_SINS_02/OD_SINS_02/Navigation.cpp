//Navigation.c
#define NAVIGATION_EXTERN

#include "stdio.h"
#include "string.h"
#include "Navigation.h"

EARTH_DATA  Earth;
SYSTEM_DATA SystemData;

//NAVI_OUT   OutputData; 

////////////////////////////
void Navigation_Init(double MINS[9])
{
	int i;

    SystemData.OutData.Att[0] = MINS[0];
    SystemData.OutData.Att[1] = MINS[1];
    SystemData.OutData.Att[2] = MINS[2];

 
	Convert_Att_To_Q(SystemData.OutData.Att,SystemData.Q );
    Convert_Q_To_Cbn(SystemData.Q,SystemData.Cb_n[0]);   //Cb2n



	SystemData.OutData.Vg[0] = MINS[3];   //�ٶ�
    SystemData.OutData.Vg[1] = MINS[4];
    SystemData.OutData.Vg[2] = MINS[5];

	SystemData.OutData.Pos[0] = MINS[6];   //λ��
    SystemData.OutData.Pos[1] = MINS[7];
    SystemData.OutData.Pos[2] = MINS[8];


	SystemData.T100MS_Flag = 0;
	SystemData.Tm          = SAMPLE_TIME_S * 2;


	for(i=0; i<3; i++)
	{
		SystemData.ThetaA2[i] = 0;
		SystemData.DeltaV2[i] = 0;
		SystemData.DV[i] = 0;
	}

    SystemData.X2Over_Count= 0;


	
	//EARTH_DATA�ṹ������
	Earth.Rs  = 0;
	Earth.Rc  = 0;
	Earth.Rmh = 0;
	Earth.Rnh = 0;
	Earth.g   = 0;
	Earth.Wz  = 0;
	Earth.Wn  = 0;
	for(i=0; i<3; i++)
	{
		Earth.Wen_n[i]  = 0;
		Earth.Win_n[i]  = 0;
		Earth.Win_b[i]  = 0;
		Earth.A_Cori[i] = 0;
	}

}

void Delta_Cal(double DTheta[3],   double Delta_Vf[3], EARTH_DATA *pEarth,
			   double DTheta_b[3], double Delta_Vg[3], 
			   double Att[3], 	    double Cb_n[3][3],  double Q[4])
{
	int i, j;
	
	/********************** �ٶ��������� ************************/

	//Delta_Vg = Cb_n * Delta_Vf
	for(i=0; i<3; i++)
    {
    	Delta_Vg[i] = 0;
    	
    	for(j=0; j<3; j++)
    	{
    		Delta_Vg[i] += Cb_n[i][j] * Delta_Vf[j];
    	}
    }

	//�������к����ٶȲ���
	Delta_Vg[0] -=  pEarth->A_Cori[0] * SystemData.Tm;
	Delta_Vg[1] -=  pEarth->A_Cori[1] * SystemData.Tm;
	Delta_Vg[2] -= (pEarth->A_Cori[2] + pEarth->g) * SystemData.Tm;

	/**********************   ��̬����    ***********************/

	//��Ԫ������
	for(i=0; i<3; i++)
	{
		DTheta_b[i] = DTheta[i] - pEarth->Win_b[i] * SystemData.Tm;
	}
	
	Refresh_Q(DTheta_b, Q);

	//��̬�������
	Convert_Q_To_Cbn(Q, Cb_n[0]);
	
	//������̬��
	Convert_Cbn_To_Att(Cb_n[0], Att);
}

/********** �����ߵ����� **********/
void SINS_Update(double DTheta[3],  double Delta_Vf[3],
			  	 double Pos[3],     double Vg[3],  double Att[3],
			  	 double Cb_n[3][3], double Q[4],   double Dins[3], EARTH_DATA *pEarth)
{
	int i;
	double Vg0[3], DTheta_b[3], Delta_Vg[3];

	//���������ز���
	Calculate_Earth(Pos, Vg, Cb_n[0], pEarth);

	//��������
	Delta_Cal(DTheta, Delta_Vf, pEarth, DTheta_b, Delta_Vg, Att, Cb_n, Q);
	
	/************************  �ٶȸ��� ************************/
	
	//��¼ǰ�����ٶ�
	for(i=0; i<3; i++)
	{
		Vg0[i] = Vg[i];
	}

	//�ٶȸ���
	for(i=0; i<3; i++)
	{
		Vg[i] += Delta_Vg[i];
	}
	
	/************************  λ�ø��� ************************/
	for(i=0; i<3; i++)
	{
		Dins[i] = (Vg0[i] + Vg[i])/2 * SystemData.Tm;
	}

	Pos[0] += Dins[0] / (pEarth->Rnh * pEarth->Rc);
	Pos[1] += Dins[1] / pEarth->Rmh;

	Pos[0] = Limit_Angle(Pos[0], 1); 
	Pos[1] = Limit_Angle(Pos[1], 0); 

	Pos[2] += Dins[2];
}

///////////////////////////////////////////////////
void Navigation(NAVI_OUT *pNaviData, double IMU1[6], double IMU2[6])
{
	int i;
    double DTheta_b[3] = {0.0},Vg0[3]= {0.0},Delta_Vg[3]= {0.0},Delta_Vgend[3]= {0.0};
    double Dw1[3] = {0.0,0.0,0.0}, Dw2[3] = {0.0,0.0,0.0};
    double Dv1[3] = {0.0,0.0,0.0}, Dv2[3] = {0.0,0.0,0.0};
    double Dw[3] = {0.0},DV[3]= {0.0};
    double dVrotm[3] = {0.0},dVsculm[3] = {0.0},dvG_Corm[3] = {0.0};
    double Temp1[3] = {0.0},Temp2[3] = {0.0},Temp3[3] = {0.0},Temp4[3] = {0.0};
	double Dins[3] = {0.0};
	double Qtmp[4] = {0.0};
	double Cn2n[3][3] = {0.0};
	double VWin_n[3] = {0.0};
    double Cb2nk_1[3][3] = {0.0};


    for(i=0; i<3; i++)
	{
	   Dw1[i] = IMU1[i];
	   Dw2[i] = IMU2[i];
	   Dv1[i] = IMU1[i+3];
	   Dv2[i] = IMU2[i+3];
	
	}

	//���������ز���
	Calculate_Earth(SystemData.OutData.Pos, SystemData.OutData.Vg, SystemData.Cb_n[0], &Earth);
    Convert_Q_To_Cbn(SystemData.Q, Cb2nk_1[0]);  //ǰһʱ�̵���̬�����ٶȽ���
/*	printf("%.16f\n%",
			   Earth.g);
*/
////////////////////////////////////////////////////////////////////////////////
	Cross(Dw1,Dw2,Temp1);  //Temp1 = cross(Dw1,Dw2);

	for(i=0; i<3; i++)
	{
		SystemData.ThetaA2[i] = Dw1[i] + Dw2[i]  + 2.0/3*Temp1[i]; //phi = IMUk.dw1 + IMUk.dw2 + 2/3*cross(IMUk.dw1,IMUk.dw2);
	}


	/**********************   ��̬����    ***********************/
	//��Ԫ������
	for(i=0; i<3; i++)
	{
		DTheta_b[i] = SystemData.ThetaA2[i] - Earth.Win_b[i] * SystemData.Tm;   //Phi + VMulQuat(Phin2nk_1,Qn2bk_1)
	}
	Refresh_Q(DTheta_b, SystemData.Q);   //
	
	//��̬�������	
	Convert_Q_To_Cbn(SystemData.Q, SystemData.Cb_n[0]);  //  Cb2n;

	//������̬��
	Convert_Cbn_To_Att(SystemData.Cb_n[0], SystemData.OutData.Att);   //Cb_n = Cb2n
	
	/********************** �ٶȸ���******************************/
	for(i=0; i<3; i++)
	{
		Vg0[i] = SystemData.OutData.Vg[i];   //ǰ�����ٶ�ֵ
	}
	
	for(i=0; i<3; i++)
    {
		Dw[i] = Dw1[i] + Dw2[i];   //Theta
        DV[i] = Dv1[i] + Dv2[i];   //dV
    }

    Cross(Dw,DV,Temp1);   //��תЧӦ���������
    dVrotm[0] = 0.5*Temp1[0];
    dVrotm[1] = 0.5*Temp1[1];
    dVrotm[2] = 0.5*Temp1[2];   //1/2*cross(Theta,dV)
    
    Cross(Dv1,Dw2,Temp2);   //cross(DV{1},DW{2})
    Cross(Dw1,Dv2,Temp3);   //cross(DW{1},DV{2})
	for(i=0; i<3; i++)
	{   
		dVsculm[i] = 2./3*(Temp2[i] + Temp3[i]);
		SystemData.DeltaV2[i] = DV[i] + dVrotm[i] + dVsculm[i] ;  //dV+dVrotm + dVsculm;
		SystemData.DV[i] = DV[i];   //Ϊ����ϵ�������fb��
	}

    Matrix31_Mult(Cb2nk_1[0],SystemData.DeltaV2,Delta_Vg);   //Cb2n*(dV+dVrotm + dVsculm)
    
    VWin_n[0] = Earth.Win_n[0];
    VWin_n[1] = Earth.Win_n[1];
    VWin_n[2] = Earth.Win_n[2];

    VMulf(&VWin_n[0], -SystemData.Tm/2.0, &Temp4[0]);
	Rv2Quat(&Temp4[0], &Qtmp[0]);

	Convert_Q_To_Cbn(Qtmp,Cn2n[0]);
    Matrix31_Mult(Cn2n[0],Delta_Vg,Delta_Vgend);   //
    /////////////////////////////////////
	dvG_Corm[0] =  - Earth.A_Cori[0] * SystemData.Tm;
	dvG_Corm[1] =  - Earth.A_Cori[1] * SystemData.Tm;
	dvG_Corm[2] =   (- Earth.g - Earth.A_Cori[2]) * SystemData.Tm;

	for(i=0; i<3; i++)
	{
		SystemData.OutData.Vg[i] += Delta_Vgend[i] + dvG_Corm[i];
	}

/*	printf("%.16f %.16f %.16f \n%",
			   SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2]);
*/
	/************************  λ�ø��� ************************/
	for(i=0; i<3; i++)
	{
		Dins[i] = (Vg0[i] + SystemData.OutData.Vg[i])/2 * SystemData.Tm;
	}

	SystemData.OutData.Pos[0] += Dins[1] / Earth.Rmh;
	SystemData.OutData.Pos[1] += Dins[0] / (Earth.Rnh * Earth.Rc);

	SystemData.OutData.Pos[0] = Limit_Angle(SystemData.OutData.Pos[0], 0); 
	SystemData.OutData.Pos[1] = Limit_Angle(SystemData.OutData.Pos[1], 1); 

	SystemData.OutData.Pos[2] += Dins[2];

	//
	*pNaviData = SystemData.OutData;   //�������
}

/**********�����Һ���**********/
double Limit_asin(double data)
{
	if(data > 1.)
	{
		data = 1.;
	}
	else if(data < -1.)
	{
		data = -1.;
	}
	
	return asin( data );
}

/**********�ǶȺ���ֵת��**********/
double Limit_Angle(double Ang, int Flag)
{
	double data = Ang;
	
	switch( Flag )
	{
	case 0:								//ת��Ϊ-PI/2 ~ +PI/2
		if(data > PI/2)
		{
			data = PI - data;
		}
		else if(data < -PI/2)
		{
			data = -PI - data;
		}
		break;
		
	case 1:								//ת��Ϊ-PI ~ +PI
		if(data > PI)
		{
			data -= 2*PI;
		}
		else if(data < -PI)
		{
			data += 2*PI;
		}
		break;
		
	case 2:								//ת��Ϊ 0 ~ +2*PI
		if(data > 2*PI)
		{
			data -= 2*PI;
		}
		else if(data < 0)
		{
			data += 2*PI;
		}
		break;
		
	default:
		break;
	}
	
	return data;
}

/**********�������**********/
void Cross(double A[3], double B[3], double C[3])
{
	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];
	C[2] = A[0]*B[1] - A[1]*B[0];
}

/**********�����ķ��Գ���**********/
void LCross(double V[3], double *pM)
{
//	M[0][0] =  0;
//	M[0][1] = -V[2];
//	M[0][2] =  V[1];

	*(pM+0) =  0;
	*(pM+1) = -V[2];
	*(pM+2) =  V[1];
	
//	M[1][0] =  V[2];
//	M[1][1] =  0;
//	M[1][2] = -V[0];

	*(pM+3) =  V[2];
	*(pM+4) =  0;
	*(pM+5) = -V[0];
	
//	M[2][0] = -V[1];
//	M[2][1] =  V[0];
//	M[2][2] =  0; 

	*(pM+6) = -V[1];
	*(pM+7) =  V[0];
	*(pM+8) =  0; 
}

/**********������׼��**********/
void norm(double Vec[3])
{
	int i;
	double m;

	m = sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);

	for(i=0; i<3; i++)
	{
		Vec[i] /= m;
	}
}

/********** 3*3����ת�� **********/
void Matrix33_Tran(double *pM)
{
	int i,j;
	double Temp[3][3] = {0.0};

	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			Temp[i][j] = *(pM+j*3+i);
		}
	}

	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			*(pM+i*3+j) = Temp[i][j];
		}
	}

}


///////////////////////////////////////////////////
double VNorm(double *Vin)
{
	return sqrt(Vin[0] * Vin[0] + Vin[1] * Vin[1] + Vin[2] * Vin[2]);
}
///////////////////////////////////////////////////

double* Rv2Quat(double *Rv, double *Qbk2bk_1)
{
	double norm = VNorm(Rv), tmp;
	
	if (norm > 1.e-20)
	{
        tmp = sin(norm/2) / (norm/2);
	}
	else
	{
		tmp = 1;
	}
	Qbk2bk_1[0] = cos(norm/2.0);
	Qbk2bk_1[1] = -tmp * Rv[0] / 2.;
	Qbk2bk_1[2] = -tmp * Rv[1] / 2.;
	Qbk2bk_1[3] = -tmp * Rv[2] / 2.;

    tmp = sqrt(Qbk2bk_1[0]*Qbk2bk_1[0] + Qbk2bk_1[1]*Qbk2bk_1[1] + Qbk2bk_1[2]*Qbk2bk_1[2] + Qbk2bk_1[3]*Qbk2bk_1[3]);

	Qbk2bk_1[0] = Qbk2bk_1[0]/tmp;
	Qbk2bk_1[1] = Qbk2bk_1[1]/tmp;
	Qbk2bk_1[2] = Qbk2bk_1[2]/tmp;
	Qbk2bk_1[3] = Qbk2bk_1[3]/tmp;

	return Qbk2bk_1;
}
///////////////////////////////////////////////////
double* VMulf(double *Vin, double f,double *Vout)
{
	Vout[0] = Vin[0]*f;
	Vout[1] = Vin[1]*f;
	Vout[2] = Vin[2]*f;

	return Vout;
}

/********** 3*3����˷� **********/
void Matrix33_Mult(double *pa, double *pb, double *pc)
{
	int i, j, k;

	for(i=0; i<3; i++)    
	{
		for(j=0; j<3; j++)
		{
			*(pc+i*3+j) = 0;
			for(k=0; k<3; k++)
			{
			    *(pc+i*3+j) += *(pa+i*3+k) * (*(pb+k*3+j));
			}
		}
	}
}

/********** 3*1����˷� **********/
void Matrix31_Mult(double *pa, double b[3], double c[3])
{
	int i, j;

	for(i=0; i<3; i++)
	{
		c[i] = 0;
		for(j=0; j<3; j++)
		{
			c[i] += *(pa+i*3+j) * b[j];
		}
	}
}

/**********��������������ٶ���Ϣ**********/
double Calculae_G(double Pos[3])
{
	double Rs,  g;
	double Rs2, Rs4, Lat, Rm, Rmh, h;
	
	Lat = Pos[1];
	h   = Pos[2];
	
	//������������
	Rs = sin(Lat);
	
	Rs2 = Rs  * Rs;
	Rs4 = Rs2 * Rs2;

	Rm  = Re * (1 - 2*e + 3*e*Rs2);
	Rmh = Rm + h;

	//�������ٶȼ���
	g = 9.7803267714*(1 + 5.27094e-3*Rs2 + 2.32718e-5*Rs4) - 3.086e-6*h;
	//g = 9.780318 * (1 + 0.0053024 * Rs2) * Rm*Rm /(Rmh*Rmh);

	return g;
}

/**********  ���������ز���  **********/
void Calculate_Earth(double Pos[3], double Vg[3], double *pCb_n, EARTH_DATA *pEarth)
{
	//pCb_n   : ��ʾCn2b
	int    i, j;
	double Rs2, Rs4, Lat, Rm, h;
	double Wie_n[3], Wcori[3];
	


	//λ������
	Lat = Pos[0];   //γ��
	h   = Pos[2];   //�߶�
	
	//�����������
	pEarth->Rs = sin(Lat);
	pEarth->Rc = cos(Lat);
	
	Rs2 = pEarth->Rs  * pEarth->Rs;
	Rs4 = Rs2 * Rs2;
    
	Rm = Re * (1 - 2*e + 3*e*Rs2);
    pEarth->Rmh = Rm + h;
    pEarth->Rnh = Re * (1 + e*Rs2) + h;

    pEarth->Wn = Wie * pEarth->Rc;
    pEarth->Wz = Wie * pEarth->Rs;
    
 	//�������ٶȼ���
	pEarth->g = 9.7803267714*(1 + 5.27094e-3*Rs2 + 2.32718e-5*Rs4) - 3.086e-6*h;

	//Wie_n����
	Wie_n[0] = 0;
	Wie_n[1] = pEarth->Wn;
	Wie_n[2] = pEarth->Wz;
	
	//Wen_n����
	pEarth->Wen_n[0] = -Vg[1]/pEarth->Rmh;
    pEarth->Wen_n[1] =  Vg[0]/pEarth->Rnh;
    pEarth->Wen_n[2] =  Vg[0]/pEarth->Rnh * tan(Lat);

	//Win_n����
	for(i=0; i<3; i++)
	{
		pEarth->Win_n[i] = Wie_n[i] + pEarth->Wen_n[i];
    }

    //Win_b����
    for(i=0; i<3; i++)
	{
        pEarth->Win_b[i] = 0;
	    for(j=0; j<3; j++)
	    {
	       pEarth->Win_b[i] += *(pCb_n+j*3+i) * pEarth->Win_n[j];	     
		}
	}
	//�к����ٶȼ���
	for(i=0; i<3; i++)
	{
	    Wcori[i] = 2*Wie_n[i] + pEarth->Wen_n[i];
    }
    
    Cross(Wcori, Vg, pEarth->A_Cori);

}
/**********����̬������̬����Cnb**********/
void Convert_Att_Cbn(double Att[3], double *pCbn) 
{
	double sinP, cosP, sinR, cosR, sinH, cosH;

	sinP = sin(Att[0]);		cosP = cos(Att[0]);
	sinR = sin(Att[1]);		cosR = cos(Att[1]);
	sinH = sin(Att[2]);		cosH = cos(Att[2]);
	
	*pCbn =  cosR*cosH + sinP*sinR*sinH;	
	*(pCbn+3) = -cosR*sinH + sinP*sinR*cosH;	
	*(pCbn+6) = -cosP*sinR;
	
	*(pCbn+1) =  cosP*sinH; 	
	*(pCbn+4) =  cosP*cosH;
	*(pCbn+7) =  sinP;
	
	*(pCbn+2) =  sinR*cosH - sinP*cosR*sinH;
	*(pCbn+5) = -sinR*sinH - sinP*cosR*cosH;
	*(pCbn+8) =  cosP*cosR;
}

/**********��̬����Ԫ��Qn2b��ת��**********/
void Convert_Att_To_Q(double Att[3], double q[4])
{
	int	i, sign_q;
	double	norm;
    double	cosH2, cosP2, cosR2, sinH2, sinP2, sinR2;
    
	//������Ԫ��ֵ
    cosP2 = cos(Att[0] / 2);
    cosR2 = cos(Att[1] / 2);
	cosH2 = cos(Att[2] / 2);
    sinP2 = sin(Att[0] / 2);
    sinR2 = sin(Att[1] / 2);
	sinH2 = sin(Att[2] / 2);
	
	q[0] = cosH2 * cosP2 * cosR2 + sinH2 * sinP2 * sinR2;
	q[1] = cosH2 * sinP2 * cosR2 + sinH2 * cosP2 * sinR2;
	q[2] = cosH2 * cosP2 * sinR2 - sinH2 * sinP2 * cosR2;
	q[3] = cosH2 * sinP2 * sinR2 - sinH2 * cosP2 * cosR2;
	
	norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	
	if(norm < 1.e-12)
	{
		norm = 1;
	}

	if (q[0]<0)
	{
		sign_q = -1;
	}
	else
	{
		sign_q = 1;
	}
	
	for(i=0; i<4; i++)
	{
		q[i] = q[i]*sign_q / norm;
	}
}

/**********��Ԫ������̬�����ת��**********/
void Convert_Q_To_Cbn(double q[4], double *pCbn)
{
	*pCbn = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	*(pCbn+1) = 2*(q[1]*q[2] - q[0]*q[3]);
	*(pCbn+2) = 2*(q[1]*q[3] + q[0]*q[2]);
	
	*(pCbn+3) = 2*(q[1]*q[2] + q[0]*q[3]);
	*(pCbn+4) = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	*(pCbn+5) = 2*(q[2]*q[3] - q[0]*q[1]);
	
	*(pCbn+6) = 2*(q[1]*q[3] - q[0]*q[2]);
	*(pCbn+7) = 2*(q[2]*q[3] + q[0]*q[1]);
	*(pCbn+8) = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

/**********����̬����Cb2n����ȡ��̬**********/
void Convert_Cbn_To_Att(double *pCbn, double Att[3])
{
	Att[0] = Limit_asin(*(pCbn+7));
    Att[1] = atan2(-(*(pCbn+6)), *(pCbn+8));
    Att[2] = atan2( *(pCbn+1), *(pCbn+4));

    if(Att[2] < 0.)
    {
		Att[2] += 2*PI;   //��������0~360��Χ
	}
}

/********** ��Ԫ������ **********/
void Refresh_Q(double DTheta_b[3], double Q[4])
{
	int i;
	double f, norm;//double f, norm;
	double MQ[4], q1[4], q[4];//double MQ[4], q1[4], q[4];
	
	//��תʸ��ת��Ϊ��Ԫ��
	norm = sqrt(DTheta_b[0]*DTheta_b[0] 
		      + DTheta_b[1]*DTheta_b[1]
			  + DTheta_b[2]*DTheta_b[2]);
	
	if(norm > 1.e-20)
	{
		f = sin(norm/2) / (norm/2);
	}
	else
	{
		f = 1;
	}
	
	MQ[0] = cos(norm/2);
	MQ[1] = DTheta_b[0] * f / 2.;
	MQ[2] = DTheta_b[1] * f / 2.;
	MQ[3] = DTheta_b[2] * f / 2.;

	for(i=0; i<4; i++)
	{
		q[i] = Q[i];
	}
	
	//��Ԫ���˷�Qk+1=Qk��q(h)=M'(q(h))*Qk
	q1[0] = MQ[0]*q[0] - MQ[1]*q[1] - MQ[2]*q[2] - MQ[3]*q[3];
	q1[1] = MQ[1]*q[0] + MQ[0]*q[1] + MQ[3]*q[2] - MQ[2]*q[3];
	q1[2] = MQ[2]*q[0] - MQ[3]*q[1] + MQ[0]*q[2] + MQ[1]*q[3];
	q1[3] = MQ[3]*q[0] + MQ[2]*q[1] - MQ[1]*q[2] + MQ[0]*q[3];
	
	//��Ԫ���淶������
	norm = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2] + q1[3]*q1[3]);
	if(norm < 1.e-20)
	{
		norm = 1;
	}
	
	for(i=0; i<4; i++)
	{
		Q[i] = q1[i] / norm;
	}
}

/**********  ��̬��������  **********/
void Attitude_Amend(double Fai[3], double *pCb_n, double Att[3], double Q[4])
{
    int i, j, k;
    double Fai_e, Fai_n, Fai_z;
	double C_g_g[3][3], Cbg[3][3];

	//Fai_enu
	Fai_e = Fai[0];
	Fai_n = Fai[1];
	Fai_z = Fai[2];
	
	//I+Fai_x
	C_g_g[0][0] = C_g_g[1][1] = C_g_g[2][2] = 1.0;
	C_g_g[0][1] = -Fai_z;
	C_g_g[0][2] =  Fai_n;
	
	C_g_g[1][0] =  Fai_z;
	C_g_g[1][2] = -Fai_e;
	
	C_g_g[2][0] = -Fai_n;
	C_g_g[2][1] =  Fai_e;
	
	//Cbg = C_g_g * Cb_n
	for(i=0; i<3; i++)    
	{
		for(j=0; j<3; j++)
		{
			Cbg[i][j] = 0;
			for(k=0; k<3; k++)
			{
//			    Cbg[i][j] += C_g_g[i][k] * Cb_n[k][j];
				Cbg[i][j] += C_g_g[i][k] * (*(pCb_n+k*3+j));
			}
		}
	}

	//Cbn = Cbg
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
//			Cb_n[i][j] = Cbg[i][j];
			*(pCb_n+i*3+j) = Cbg[i][j];
		}
	}
	
	//���������¼�����̬��
    Convert_Cbn_To_Att(pCb_n, Att);
	
	//���¼�����Ԫ��
	Convert_Att_To_Q(Att, Q);
}


/**********��������(6�׾���)**********/
int Matrix_Inversesix(double M[6][6])
{
	int i, j, k, i0, j0, p0[6], q0[6];
	double w, z;

	//����
	for(k=1; k<=6; k++)
	{
		w = 0.0;
		
		for(i=k; i<=6; i++)
		{
			for(j=k; j<=6; j++)
			{
				if( fabs(M[i-1][j-1]) > fabs(w))
				{
					w  = M[i-1][j-1];
					i0 = i;
					j0 = j;
				}
			}
		}
		
		if(fabs(w) <= 1e-20)
		{
			return 0;
		}
		
		//exchange line elements
	 	if(i0 != k)
	 	{
			for(j=1; j<=6; j++)
			{
				z = M[i0-1][j-1];
				M[i0-1][j-1] = M[k-1][j-1];
				M[k-1][j-1]  = z;
			}
		}
	
		//exchange column elements
		if(j0 != k)
		{
			for(i=1; i<=6; i++)
			{
				z = M[i-1][j0-1];
				M[i-1][j0-1] = M[i-1][k-1];
				M[i-1][k-1]  = z;
			}
		}
	
		p0[k-1] = i0;
		q0[k-1] = j0;
		M[k-1][k-1] = 1 / M[k-1][k-1];
	
	    for(j=1; j<=6; j++)
	    {
			if(j != k)
			{
				M[k-1][j-1]=M[k-1][j-1] * M[k-1][k-1];
			}
	    }
    
	    for(i=1; i<=6; i++)
	    {
			if(i != k )
			{
				for(j=1; j<=6; j++)
				{
					if(j != k)
					{
						M[i-1][j-1]=M[i-1][j-1] - M[i-1][k-1] * M[k-1][j-1];
					}
				}
			}
	    }
    
		for(i=1; i<=6; i++)
		{
			if(i != k)
			{
				M[i-1][k-1] = -M[i-1][k-1] * M[k-1][k-1];
			}
		}
	}

	//Restore the elements according to the record in matrix p0&q0
	for(k=6; k>=1; k--)
	{
		i0 = p0[k-1];
		j0 = q0[k-1];
		
		if(j0 != k)
		{
			for(j=1; j<=6; j++)
			{
				z = M[k-1][j-1];
				M[k-1][j-1]  = M[j0-1][j-1];
				M[j0-1][j-1] = z;
			}
		}
		
		if(i0 != k)
		{
			for(i=1; i<=6; i++)
			{
				z = M[i-1][i0-1];
				M[i-1][i0-1] = M[i-1][k-1];
				M[i-1][k-1]  = z;
			}
		}
	}
	
	return 1;
}




/**********��������**********/
int Matrix_Inverse(double *pM)
{
	int i, j, k, i0, j0, p0[3], q0[3];
	double w, z;

	//����
	for(k=1; k<=3; k++)
	{
		w = 0.0;
		
		for(i=k; i<=3; i++)
		{
			for(j=k; j<=3; j++)
			{ 
				
				if( fabs(*(pM +(i-1)*3+(j-1))) > fabs(w))
				{
					w  = *(pM +(i-1)*3+(j-1)); 
					i0 = i;
					j0 = j;
				}
			}
		}
		
		if(fabs(w) <= 1e-20)
		{
			return 0;
		}
		
		//exchange line elements
	 	if(i0 != k)
	 	{
			for(j=1; j<=3; j++)
			{
				z = *(pM+(i0-1)*3+(j-1));
				*(pM+(i0-1)*3+(j-1)) = *(pM+(k-1)*3+(j-1));
				*(pM+(k-1)*3+(j-1))  = z;
			}
		}
	
		//exchange column elements
		if(j0 != k)
		{
			for(i=1; i<=3; i++)
			{
				z = *(pM+(i-1)*3+(j0-1)); //M[i-1][j0-1];
				*(pM+(i-1)*3+(j0-1)) = *(pM+(i-1)*3+(k-1)); //M[i-1][j0-1] = M[i-1][k-1];
				*(pM+(i-1)*3+(k-1))  = z; //M[i-1][k-1] = z;
			}
		}
	
		p0[k-1] = i0;
		q0[k-1] = j0;
		*(pM+(k-1)*3+(k-1)) = 1.0 / *(pM+(k-1)*3+(k-1)); //M[k-1][k-1] = 1 / M[k-1][k-1];
	
	    for(j=1; j<=3; j++)
	    {
			if(j != k)
			{
				*(pM+(k-1)*3+(j-1)) = *(pM+(k-1)*3+(j-1)) * (*(pM+(k-1)*3+(k-1)));
				//M[k-1][j-1]=M[k-1][j-1] * M[k-1][k-1];
			}
	    }
    
	    for(i=1; i<=3; i++)
	    {
			if(i != k )
			{
				for(j=1; j<=3; j++)
				{
					if(j != k)
					{
						*(pM+(i-1)*3+(j-1)) = *(pM+(i-1)*3+(j-1)) - (*(pM+(i-1)*3+(k-1))) * (*(pM+(k-1)*3+(j-1)));
//						M[i-1][j-1]=M[i-1][j-1] - M[i-1][k-1] * M[k-1][j-1];
					}
				}
			}
	    }
    
		for(i=1; i<=3; i++)
		{
			if(i != k)
			{
				*(pM+(i-1)*3+(k-1)) = -(*(pM+(i-1)*3+(k-1))) * (*(pM+(k-1)*3+(k-1)));
				//M[i-1][k-1] = -M[i-1][k-1] * M[k-1][k-1];
			}
		}
	}

	//Restore the elements according to the record in matrix p0&q0
	for(k=3; k>=1; k--)
	{
		i0 = p0[k-1];
		j0 = q0[k-1];
		
		if(j0 != k)
		{
			for(j=1; j<=3; j++)
			{
				z = *(pM + (k-1)*3 + (j-1));  //z = M[k-1][j-1];
				*(pM + (k-1)*3 + (j-1)) = *(pM + (j0-1)*3 + (j-1));  //M[k-1][j-1]  = M[j0-1][j-1];
				*(pM + (j0-1)*3 + (j-1)) = z;  //M[j0-1][j-1] = z;
			}
		}
		
		if(i0 != k)
		{
			for(i=1; i<=3; i++)
			{
				z = *(pM + (i-1)*3 + (i0-1));  //z = M[i-1][i0-1];
				*(pM + (i-1)*3 + (i0-1)) = *(pM + (i-1)*3 + (k-1)); //M[i-1][i0-1] = M[i-1][k-1];
				*(pM + (i-1)*3 + (k-1)) = z;  //M[i-1][k-1]  = z;
			}
		}
	}
	
	return 1;
}



/********** ˮƽλ�þ������ **********/
double HPos_Dis(double Pos0[3], double Pos1[3])
{
	double Reh, Rc;
	double Dpos[2];
	
	Reh = Re + Pos0[2];
	Rc  = cos(Pos0[1]);
	
	Dpos[0] = Limit_Angle(Pos0[0] - Pos1[0], 1) * Reh * Rc;
	Dpos[1] = (Pos0[1] - Pos1[1]) * Reh;
	
	return sqrt(Dpos[0]*Dpos[0] + Dpos[1]*Dpos[1]);
}


/********** ��Ԫ���˷�����Q3=Q2*Q1 **********/
void Quat_Mul(double Q1[4], double Q2[4], double Q3[4])
{
	double QR[4] = {0}, norm = 0;
	//��Ԫ���˷�Qk+1=Qk��q(h)=M'(q(h))*Qk
	QR[0] = Q1[0]*Q2[0] - Q1[1]*Q2[1] - Q1[2]*Q2[2] - Q1[3]*Q2[3];
	QR[1] = Q1[1]*Q2[0] + Q1[0]*Q2[1] + Q1[3]*Q2[2] - Q1[2]*Q2[3];
	QR[2] = Q1[2]*Q2[0] - Q1[3]*Q2[1] + Q1[0]*Q2[2] + Q1[1]*Q2[3];
	QR[3] = Q1[3]*Q2[0] + Q1[2]*Q2[1] - Q1[1]*Q2[2] + Q1[0]*Q2[3];
	
	//��Ԫ���淶������ / norm
/**/
	norm = sqrt(QR[0]*QR[0] + QR[1]*QR[1] + QR[2]*QR[2] + QR[3]*QR[3]);
	if(norm < 1.e-20)
	{
		norm = 1;
	}
	
	for(int i=0; i<4; i++)
	{
		Q3[i] = QR[i] / norm;
	}

}
/*******************************************************/
/*******************************************************/

// �������׷��������ʽ
double sdet_33(double *pM)
{
	double det = 0.0;

	det = *(pM) * ( *(pM+4) * (*(pM+8)) - *(pM+5) * (*(pM+7)) ) - 
		  *(pM+1) * ( *(pM+3) * (*(pM+8)) - *(pM+5) * (*(pM+6)) ) + 
		  *(pM+2) * ( *(pM+3) * (*(pM+7)) - *(pM+4) * (*(pM+6)) );

	return det;
}

/*******************************************************/
/*******************************************************/

/*	Fun describe: Phi = Quat2Rv(Qb2n*Qn2bc)
	Para in		: Qb2nc, Qb2n
	Para out	: Phi 	*/
void QQ2Phi(double Qb2nc[4], double Qb2n[4], double Phi[3])//Phi = Qb2n*Qb2nc'
{
	double Qtmp[4];
	double dQ[4];
    Qtmp[0] = Qb2nc[0];
	Qtmp[1] = -Qb2nc[1];
	Qtmp[2] = -Qb2nc[2];
	Qtmp[3] = -Qb2nc[3];
		
//	Quat_Mul(Qb2n, Qtmp, dQ);

	dQ[0] = Qb2n[0]*Qtmp[0] - Qb2n[1]*Qtmp[1] - Qb2n[2]*Qtmp[2] - Qb2n[3]*Qtmp[3];
	dQ[1] = Qb2n[1]*Qtmp[0] + Qb2n[0]*Qtmp[1] - Qb2n[3]*Qtmp[2] + Qb2n[2]*Qtmp[3];
	dQ[2] = Qb2n[2]*Qtmp[0] + Qb2n[3]*Qtmp[1] + Qb2n[0]*Qtmp[2] - Qb2n[1]*Qtmp[3];
	dQ[3] = Qb2n[3]*Qtmp[0] - Qb2n[2]*Qtmp[1] + Qb2n[1]*Qtmp[2] + Qb2n[0]*Qtmp[3];
	
	

/*  printf("%lf %lf %lf %lf \n%lf %lf %lf %lf \n%lf %lf %lf %lf \n",
      Qb2n[0],Qb2n[1],Qb2n[2],Qb2n[3], 
	  Qtmp[0],Qtmp[1],Qtmp[2],Qtmp[3], 
	  dQ[0],dQ[1],dQ[2],dQ[3]);*/
		
	Quat2Rv(dQ, Phi);	
}

/*	Fun describe: Rv = Quat2Rv(Qb2n)
	Para in		: Qb2n
	Para out	: Rv 	*/
//�ο�����Ԫ�����Ե�����p306, 9.3.4�������
void Quat2Rv(double Qb2n[4], double Rv[3])//Quaternion to rv
{
	double n2, k;

	if(Qb2n[0]<0)
	{
		Qb2n[0] = -Qb2n[0];
		Qb2n[1] = -Qb2n[1];
		Qb2n[2] = -Qb2n[2];
		Qb2n[3] = -Qb2n[3];
	}
	
    n2 = acos(Qb2n[0]);
    if (n2 > 1e-100)
        k = 2*n2/sin(n2);
    else
        k = 2;
	
    //VInit(Rv, k*Qb2n[1], k*Qb2n[2], k*Qb2n[3]);
    Rv[0] = k*Qb2n[1];
    Rv[1] = k*Qb2n[2];
	Rv[2] = k*Qb2n[3];
//	return Rv;
}
