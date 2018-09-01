//Navigation.c
#define NAVIGATION_EXTERN

#include "stdio.h"
#include "string.h"
#include "Navigation.h"

EARTH_DATA  Earth;
SYSTEM_DATA SystemData;

//NAVI_OUT   OutputData; 

////////////////////////////
void Navigation_Init(DBL MINS[9])
{
	int i;

    SystemData.OutData.Att[0] = MINS[0];
    SystemData.OutData.Att[1] = MINS[1];
    SystemData.OutData.Att[2] = MINS[2];

 
	Convert_Att_To_Q(SystemData.OutData.Att,SystemData.Q );
    Convert_Q_To_Cbn(SystemData.Q,SystemData.Cb_n[0]);   //Cb2n



	SystemData.OutData.Vg[0] = MINS[3];   //速度
    SystemData.OutData.Vg[1] = MINS[4];
    SystemData.OutData.Vg[2] = MINS[5];

	SystemData.OutData.Pos[0] = MINS[6];   //位置
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


	
	//EARTH_DATA结构体清零
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

/********** Kalman滤波器初始化 **********/
void Init_KF_Filter()
{
	int i, j;

	//变量清零
	SystemData.TransCount = 0;
	SystemData.TransFlag  = 0;   //滤波器时间更新标志

	for(i=0; i<STS15; i++)
	{
		SystemData.Xk[i] = 0;
		SystemData.Qt[i] = 0;
        SystemData.DFai[i] = 0;
		for(j=0; j<STS15; j++)
		{
			SystemData.Pk[i][j] = 0;
			SystemData.Ft[i][j] = 0;
		}
	}

	for(i=0; i<WAT; i++)   //WAT3:量测维数
	{
		for(j=0; j<WAT; j++)
		{
			SystemData.R1[i][j] = 0;   //INS/GPS
		}
		
		for(j=0; j<STS15; j++)
		{
			SystemData.H1[i][j] = 0;   //INS/GPS
		}
	}


	
	SystemData.NoFilterTime = MAX_TIME;   //连续未滤波时间0~1：对准、INS/GPS
	SystemData.FilterCount  = 0;//滤波次数
	SystemData.MesFlag[0]      = 0;//滤波计算标志
	
    //P阵初值
    SystemData.Pk[0][0]   = pow(5*Deg_Rad, 2);
	SystemData.Pk[1][1]   = pow(5*Deg_Rad, 2);
	SystemData.Pk[2][2]   = pow(10*Deg_Rad, 2);
	
	SystemData.Pk[3][3]   = pow(1.0, 2);
	SystemData.Pk[4][4]   = pow(1.0, 2);
	SystemData.Pk[5][5]   = pow(1.0, 2);
	
	SystemData.Pk[6][6]   = pow(10*Dph_Rps, 2);
	SystemData.Pk[7][7]   = pow(10*Dph_Rps, 2);
	SystemData.Pk[8][8]   = pow(10*Dph_Rps, 2);
	
	SystemData.Pk[9][9]   = pow(1*mg, 2);
	SystemData.Pk[10][10] = pow(1*mg, 2);
	SystemData.Pk[11][11] = pow(1*mg, 2);

	SystemData.Pk[12][12] = 5.0*Deg_Rad;
	SystemData.Pk[13][13] = 5.0*Deg_Rad;
	SystemData.Pk[14][14] = 5.0*Deg_Rad;

	//Q阵初值
	SystemData.Qt[0]   = pow(0.1*dp05h, 2);
	SystemData.Qt[1]   = pow(0.1*dp05h, 2);
    SystemData.Qt[2]   = pow(0.1*dp05h, 2);

	SystemData.Qt[3]   = pow(0.1*mg, 2);
	SystemData.Qt[4]   = pow(0.1*mg, 2);
	SystemData.Qt[5]   = pow(0.1*mg, 2);
    
 /*   SystemData.Qt[6]   = pow(0./Re, 2);
	SystemData.Qt[7]   = pow(0./Re, 2);
	SystemData.Qt[8]   = pow(0., 2);  */


	//INS/GPS组合
	
	SystemData.H1[0][0]  =  1;
	SystemData.H1[1][1]  =  1;
	SystemData.H1[2][2]  =  1;
	
	SystemData.H1[3][3]  =  1;
	SystemData.H1[4][4]  =  1;
	SystemData.H1[5][5]  =  1;


	SystemData.R1[0][0] = pow(2*Min_Rad, 2);
	SystemData.R1[1][1] = pow(2*Min_Rad, 2);
	SystemData.R1[2][2] = pow(4*Min_Rad, 2);
//
	SystemData.R1[3][3] = pow(0.1, 2);
	SystemData.R1[4][4] = pow(0.1, 2);
	SystemData.R1[5][5] = pow(0.1, 2);
}

/********** Kalman滤波器F阵计算 **********/
void Cal_KF_FNav(DBL Lat, DBL Vg[3], DBL DeltaVb[3],
				 DBL Cbn[3][3],      EARTH_DATA earth)   //Cbn = Cb2n
{    
	int i, j;

	DBL tanL, Wz, Wn, Rs, Rc;
	DBL VnR, VeR, VeLR, Rmh, Rnh;
	DBL Fg[3];


	for(i=0; i<3; i++)
	{
		Fg[i] = 0;
		
		for(j=0; j<3; j++)
		{
			Fg[i] += Cbn[i][j] * DeltaVb[j] * TM_COUNT;
		}		
	}

/*	printf("%.16f %.16f %.16f %.16f\n%",
			   Fg[0],Fg[1],Fg[2],DeltaVb[2]);*/



	memset((void*)SystemData.Ft, 0, STS15*STS15*sizeof(DBL));

	//控制每100ms进行
	if( SystemData.T100MS_Flag )
	{

		//变量引用
		Wn  = earth.Wn;   //
		Wz  = earth.Wz;
		Rs  = earth.Rs;
		Rc  = earth.Rc;
		Rmh = earth.Rmh;
		Rnh = earth.Rnh;

		tanL = tan( Lat );   //Lat：纬度

		VnR  = Vg[1]/Rmh;
		VeR  = Vg[0]/Rnh;
		VeLR = Vg[0]*tanL/Rnh;
		

        SystemData.Ft[0][1]   = Wz + VeLR;   //wie*sinL + Ve*tanl/RNh
        SystemData.Ft[0][2]   = -(Wn + VeR); //-(wie*cosL + Ve/RNh)


        SystemData.Ft[1][0]   =  -(Wz + VeLR);  //-(wie*sinL + Ve*tanl/RNh)
        SystemData.Ft[1][2]   =  -VnR;   //-Vn/RMh

        SystemData.Ft[2][0]   =  Wn + VeR;   //(wie*cosL + Ve/RNh)
        SystemData.Ft[2][1]   =  VnR;   //Vn/RMh

		SystemData.Ft[0][4]   = -1./Rmh;   //-1/RMh
        SystemData.Ft[1][3]   = 1./Rnh;  //1/RNh

        SystemData.Ft[2][3]   = tanL/Rnh; //tanL/RNh

      //  SystemData.Ft[0][8]   =  VnR/Rmh;
       // SystemData.Ft[1][6]   = -Wz;
      //  SystemData.Ft[1][8]   = -VeR/Rnh;

		SystemData.Ft[0][6]  = - Cbn[0][0];
        SystemData.Ft[0][7]  = - Cbn[0][1];
        SystemData.Ft[0][8]  = - Cbn[0][2];

		SystemData.Ft[1][6]  = - Cbn[1][0];
		SystemData.Ft[1][7]  = - Cbn[1][1];
		SystemData.Ft[1][8]  = - Cbn[1][2];

		SystemData.Ft[2][6]  = - Cbn[2][0];
		SystemData.Ft[2][7]  = - Cbn[2][1];
		SystemData.Ft[2][8]  = - Cbn[2][2];


	//	SystemData.Ft[2][6]   =  Wn + VeR/Rc/Rc;
	//	SystemData.Ft[2][8]   = -VeLR/Rnh;
		
        SystemData.Ft[3][1]   = -Fg[2];   //-fz_n
        SystemData.Ft[3][2]   =  Fg[1];   //fy_n

        SystemData.Ft[4][0]   =  Fg[2];   //fz_n
        SystemData.Ft[4][2]   = -Fg[0];  //-fx_n

	    SystemData.Ft[5][0]   = -Fg[1];   //-fy_n
		SystemData.Ft[5][1]   =  Fg[0];   //fx_n


        SystemData.Ft[3][3]   = (Vg[1] * tanL - Vg[2])/Rnh;  //(Vn*tanL - VU)/RNh
        SystemData.Ft[3][4]   = 2 * Wz + VeLR;   //2wie*sinl + Ve*tanL/RNh
        SystemData.Ft[3][5]   = -(2*Wn + VeR);   //-(2*wie*cosL+ Ve/Rnh)

        SystemData.Ft[4][3]   = -2 * (Wz + VeLR);  //-2*(wie*sinL + Ve*tanL/Rnh)
		SystemData.Ft[4][4]   = -Vg[2]/Rmh;   //-VU/RMh
		SystemData.Ft[4][5]   = -VnR;    //-Vn/RMh


	    SystemData.Ft[5][3]   =  2*(Wn + VeR);   //2*(wie*cosL + Ve/Rnh)
		SystemData.Ft[5][4]   =  2*VnR;   //2*Vn/RMh

	    SystemData.Ft[3][9]  = Cbn[0][0];
		SystemData.Ft[3][10]  = Cbn[0][1];
	    SystemData.Ft[3][11]  = Cbn[0][2];

	    SystemData.Ft[4][9]  =  Cbn[1][0];
		SystemData.Ft[4][10]  =  Cbn[1][1];
		SystemData.Ft[4][11]  =  Cbn[1][2];

		SystemData.Ft[5][9]   = Cbn[2][0];
	    SystemData.Ft[5][10]  = Cbn[2][1];
		SystemData.Ft[5][11]  = Cbn[2][2];

	//	SystemData.Ft[6][4]   = 1./Rmh;   //1/RMh
     //   SystemData.Ft[7][3]   = 1./Rc/Rnh;  //secL/Rnh
     //   SystemData.Ft[8][5]   = 1;


      //  SystemData.Ft[6][8]   = -VnR/Rmh;   //-Vn/(RMh)^2
	//	SystemData.Ft[7][6]   = VeLR/Rc;
	//	SystemData.Ft[7][8]   = -VeR/Rc/Rnh;

    //    SystemData.Ft[3][8]   = Vg[0]*(Vg[2] - Vg[1]*tanL)/Rnh/Rnh;
    //    SystemData.Ft[3][6]   = 2*Wie*(Vg[2]*Rs+Vg[1]*Rc) + VeR*Vg[1]/Rc/Rc;

	 //   SystemData.Ft[4][6]   = -Vg[0]*(2*Wn  + VeR/Rc/Rc);
	//	SystemData.Ft[4][8]   = Vg[2]*VnR/Rmh + Vg[0]*VeLR/Rnh;

	//	SystemData.Ft[5][6]   = -2*Vg[0]*Wz;
	//	SystemData.Ft[5][8]   = -(VnR*VnR + VeR*VeR);


		/*
		SystemData.Ft[9][9]   = -1/GYRO_TAO;
		SystemData.Ft[10][10] = -1/GYRO_TAO;
		SystemData.Ft[11][11] = -1/GYRO_TAO;
		*/

/*printf("%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f %.16f %.16f %.16f\n  %.16f %.16f %.16f %.16f %.16f %.16f\n "
		"%.16f %.16f %.16f\n"
        "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n",
		

	    Cbn[0][0],Cbn[0][1],Cbn[0][2],
        Cbn[1][0],Cbn[1][1],Cbn[1][2],
		Cbn[2][0],Cbn[2][1],Cbn[2][2],
	    SystemData.Ft[0][0],SystemData.Ft[0][1],SystemData.Ft[0][2],SystemData.Ft[0][3],
		SystemData.Ft[0][4],SystemData.Ft[0][5],SystemData.Ft[0][6],
		SystemData.Ft[0][7],SystemData.Ft[0][8],SystemData.Ft[1][0],SystemData.Ft[1][1],SystemData.Ft[1][2],SystemData.Ft[1][3],
		SystemData.Ft[1][4],SystemData.Ft[1][5],SystemData.Ft[1][6],SystemData.Ft[1][7],SystemData.Ft[1][8],SystemData.Ft[2][0],
		SystemData.Ft[2][1],SystemData.Ft[2][2],SystemData.Ft[2][3],SystemData.Ft[2][4],SystemData.Ft[2][5],SystemData.Ft[2][6],
		SystemData.Ft[2][7],SystemData.Ft[2][8],SystemData.Ft[3][0],SystemData.Ft[3][1],SystemData.Ft[3][2],
		SystemData.Ft[3][3],SystemData.Ft[3][4],SystemData.Ft[3][5],SystemData.Ft[3][9],SystemData.Ft[3][10],SystemData.Ft[3][11],
		SystemData.Ft[6][0],SystemData.Ft[6][1],SystemData.Ft[6][2],SystemData.Ft[6][3],SystemData.Ft[6][4],SystemData.Ft[6][5],
		SystemData.Ft[7][0],SystemData.Ft[7][1],SystemData.Ft[7][2],SystemData.Ft[7][3],SystemData.Ft[7][4],SystemData.Ft[7][5],
		SystemData.Ft[8][5]);   */
				 
		SystemData.TransFlag  = 1;   //状态转移标志每0.1s置1一次为了时间更新准备
	}
}

/********** Kalman滤波器时间更新 **********/
void KF_TimeUpdate()
{
	int i, j, k;
	
	DBL I0, Xkk_1[STS15];
	DBL Fk[STS15][STS15];
	DBL pkmid1[STS15][STS15];

 /*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n"

		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n\n\n"
		
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n",
	    SystemData.Pk[0][0],SystemData.Pk[0][1],SystemData.Pk[0][2],SystemData.Pk[0][3],
		SystemData.Pk[0][4],SystemData.Pk[0][5],SystemData.Pk[1][0],
		SystemData.Pk[1][1],SystemData.Pk[1][2],SystemData.Pk[1][3],SystemData.Pk[1][4],SystemData.Pk[1][5],SystemData.Pk[2][0],
		SystemData.Pk[2][1],SystemData.Pk[2][2],SystemData.Pk[2][3],SystemData.Pk[2][4],SystemData.Pk[2][5],SystemData.Pk[3][0],
		SystemData.Pk[3][1],SystemData.Pk[3][2],SystemData.Pk[3][3],SystemData.Pk[3][4],SystemData.Pk[3][5],SystemData.Pk[4][0],
		SystemData.Pk[4][1],SystemData.Pk[4][2],SystemData.Pk[4][3],SystemData.Pk[4][4],SystemData.Pk[4][5],
		SystemData.Pk[5][0],SystemData.Pk[5][1],SystemData.Pk[5][2],SystemData.Pk[5][3],SystemData.Pk[5][4],SystemData.Pk[5][5],
		SystemData.Pk[6][0],SystemData.Pk[6][1],SystemData.Pk[6][2],SystemData.Pk[6][3],SystemData.Pk[6][4],SystemData.Pk[6][5],
		SystemData.Pk[7][0],SystemData.Pk[7][1],SystemData.Pk[7][2],SystemData.Pk[7][3],SystemData.Pk[7][4],SystemData.Pk[7][5],
		SystemData.Pk[8][5],
		SystemData.Pk[9][0],SystemData.Pk[9][1],SystemData.Pk[9][2],SystemData.Pk[9][3],SystemData.Pk[9][4],SystemData.Pk[9][5],
		SystemData.Pk[10][0],SystemData.Pk[10][1],SystemData.Pk[10][2],SystemData.Pk[10][3],SystemData.Pk[10][4],SystemData.Pk[10][5],
		SystemData.Pk[11][0],SystemData.Pk[11][1],SystemData.Pk[11][2],SystemData.Pk[11][3],SystemData.Pk[11][4],SystemData.Pk[11][5],
		SystemData.Pk[12][0],SystemData.Pk[12][1],SystemData.Pk[12][2],SystemData.Pk[12][3],SystemData.Pk[12][4],SystemData.Pk[12][5],
		SystemData.Pk[13][0],SystemData.Pk[13][1],SystemData.Pk[13][2],SystemData.Pk[13][3],SystemData.Pk[13][4],SystemData.Pk[13][5],
		SystemData.Pk[14][0],SystemData.Pk[14][1],SystemData.Pk[14][2],SystemData.Pk[14][3],SystemData.Pk[14][4],SystemData.Pk[14][5],
		
		SystemData.Pk[0][6],SystemData.Pk[0][7],SystemData.Pk[0][8],SystemData.Pk[0][9],SystemData.Pk[0][10],SystemData.Pk[0][11],
		SystemData.Pk[1][6],SystemData.Pk[1][7],SystemData.Pk[1][8],SystemData.Pk[1][9],SystemData.Pk[1][10],SystemData.Pk[1][11],
		SystemData.Pk[2][6],SystemData.Pk[2][7],SystemData.Pk[2][8],SystemData.Pk[2][9],SystemData.Pk[2][10],SystemData.Pk[2][11],
		SystemData.Pk[3][6],SystemData.Pk[3][7],SystemData.Pk[3][8],SystemData.Pk[3][9],SystemData.Pk[3][10],SystemData.Pk[3][11],
		SystemData.Pk[4][6],SystemData.Pk[4][7],SystemData.Pk[4][8],SystemData.Pk[4][9],SystemData.Pk[4][10],SystemData.Pk[4][11],
		SystemData.Pk[5][6],SystemData.Pk[5][7],SystemData.Pk[5][8],SystemData.Pk[5][9],SystemData.Pk[5][10],SystemData.Pk[5][11],
		SystemData.Pk[6][6],SystemData.Pk[6][7],SystemData.Pk[6][8],SystemData.Pk[6][9],SystemData.Pk[6][10],SystemData.Pk[6][11],
		SystemData.Pk[7][6],SystemData.Pk[7][7],SystemData.Pk[7][8],SystemData.Pk[7][9],SystemData.Pk[7][10],SystemData.Pk[7][11],
		SystemData.Pk[8][6],SystemData.Pk[8][7],SystemData.Pk[8][8],SystemData.Pk[8][9],SystemData.Pk[8][10],SystemData.Pk[8][11],
		SystemData.Pk[9][6],SystemData.Pk[9][7],SystemData.Pk[9][8],SystemData.Pk[9][9],SystemData.Pk[9][10],SystemData.Pk[9][11],
		SystemData.Pk[10][6],SystemData.Pk[10][7],SystemData.Pk[10][8],SystemData.Pk[10][9],SystemData.Pk[10][10],SystemData.Pk[10][11],
		SystemData.Pk[11][6],SystemData.Pk[11][7],SystemData.Pk[11][8],SystemData.Pk[11][9],SystemData.Pk[11][10],SystemData.Pk[11][11],
		SystemData.Pk[12][6],SystemData.Pk[12][7],SystemData.Pk[12][8],SystemData.Pk[12][9],SystemData.Pk[12][10],SystemData.Pk[12][11],
		SystemData.Pk[13][6],SystemData.Pk[13][7],SystemData.Pk[13][8],SystemData.Pk[13][9],SystemData.Pk[13][10],SystemData.Pk[13][11],
		SystemData.Pk[14][6],SystemData.Pk[14][7],SystemData.Pk[14][8],SystemData.Pk[14][9],SystemData.Pk[14][10],SystemData.Pk[14][11],
		

	    SystemData.Pk[0][12],SystemData.Pk[0][13],SystemData.Pk[0][14],
		SystemData.Pk[1][12],SystemData.Pk[1][13],SystemData.Pk[1][14],
		SystemData.Pk[2][12],SystemData.Pk[2][13],SystemData.Pk[2][14],
		SystemData.Pk[3][12],SystemData.Pk[3][13],SystemData.Pk[3][14],
		SystemData.Pk[4][12],SystemData.Pk[4][13],SystemData.Pk[4][14],
		SystemData.Pk[5][12],SystemData.Pk[5][13],SystemData.Pk[5][14],
		SystemData.Pk[6][12],SystemData.Pk[6][13],SystemData.Pk[6][14],
		SystemData.Pk[7][12],SystemData.Pk[7][13],SystemData.Pk[7][14],
		SystemData.Pk[8][12],SystemData.Pk[8][13],SystemData.Pk[8][14],
		SystemData.Pk[9][12],SystemData.Pk[9][13],SystemData.Pk[9][14],
		SystemData.Pk[10][12],SystemData.Pk[10][13],SystemData.Pk[10][14],
		SystemData.Pk[11][12],SystemData.Pk[11][13],SystemData.Pk[11][14],
		SystemData.Pk[12][12],SystemData.Pk[12][13],SystemData.Pk[12][14],
		SystemData.Pk[13][12],SystemData.Pk[13][13],SystemData.Pk[13][14],
		SystemData.Pk[14][12],SystemData.Pk[14][13],SystemData.Pk[14][14]);

*/





	if( SystemData.TransFlag )//SystemData.TransFlag :状态转移标志(时间跟新标志)每0.1s时间更新周期到了被置1
	{
		//状态转移矩阵离散化(Fk = I + F*T)
		for(i=0; i<STS15; i++)
		{
			for(j=0; j<STS15; j++)
			{
				I0 = (i == j);
				
				Fk[i][j] = I0 + SystemData.Ft[i][j] * TF;
			}
		}
//	printf("%.16f\n", TF);

/*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n"

		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n\n\n"
		
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n",
	    Fk[0][0],Fk[0][1],Fk[0][2],Fk[0][3],
		Fk[0][4],Fk[0][5],Fk[1][0],
		Fk[1][1],Fk[1][2],Fk[1][3],Fk[1][4],Fk[1][5],Fk[2][0],
		Fk[2][1],Fk[2][2],Fk[2][3],Fk[2][4],Fk[2][5],Fk[3][0],
		Fk[3][1],Fk[3][2],Fk[3][3],Fk[3][4],Fk[3][5],Fk[4][0],
		Fk[4][1],Fk[4][2],Fk[4][3],Fk[4][4],Fk[4][5],
		Fk[5][0],Fk[5][1],Fk[5][2],Fk[5][3],Fk[5][4],Fk[5][5],
		Fk[6][0],Fk[6][1],Fk[6][2],Fk[6][3],Fk[6][4],Fk[6][5],
		Fk[7][0],Fk[7][1],Fk[7][2],Fk[7][3],Fk[7][4],Fk[7][5],
		Fk[8][5],
		Fk[9][0],Fk[9][1],Fk[9][2],Fk[9][3],Fk[9][4],Fk[9][5],
		Fk[10][0],Fk[10][1],Fk[10][2],Fk[10][3],Fk[10][4],Fk[10][5],
		Fk[11][0],Fk[11][1],Fk[11][2],Fk[11][3],Fk[11][4],Fk[11][5],
		Fk[12][0],Fk[12][1],Fk[12][2],Fk[12][3],Fk[12][4],Fk[12][5],
		Fk[13][0],Fk[13][1],Fk[13][2],Fk[13][3],Fk[13][4],Fk[13][5],
		Fk[14][0],Fk[14][1],Fk[14][2],Fk[14][3],Fk[14][4],Fk[14][5],
		
		Fk[0][6],Fk[0][7],Fk[0][8],Fk[0][9],Fk[0][10],Fk[0][11],
		Fk[1][6],Fk[1][7],Fk[1][8],Fk[1][9],Fk[1][10],Fk[1][11],
		Fk[2][6],Fk[2][7],Fk[2][8],Fk[2][9],Fk[2][10],Fk[2][11],
		Fk[3][6],Fk[3][7],Fk[3][8],Fk[3][9],Fk[3][10],Fk[3][11],
		Fk[4][6],Fk[4][7],Fk[4][8],Fk[4][9],Fk[4][10],Fk[4][11],
		Fk[5][6],Fk[5][7],Fk[5][8],Fk[5][9],Fk[5][10],Fk[5][11],
		Fk[6][6],Fk[6][7],Fk[6][8],Fk[6][9],Fk[6][10],Fk[6][11],
		Fk[7][6],Fk[7][7],Fk[7][8],Fk[7][9],Fk[7][10],Fk[7][11],
		Fk[8][6],Fk[8][7],Fk[8][8],Fk[8][9],Fk[8][10],Fk[8][11],
		Fk[9][6],Fk[9][7],Fk[9][8],Fk[9][9],Fk[9][10],Fk[9][11],
		Fk[10][6],Fk[10][7],Fk[10][8],Fk[10][9],Fk[10][10],Fk[10][11],
		Fk[11][6],Fk[11][7],Fk[11][8],Fk[11][9],Fk[11][10],Fk[11][11],
		Fk[12][6],Fk[12][7],Fk[12][8],Fk[12][9],Fk[12][10],Fk[12][11],
		Fk[13][6],Fk[13][7],Fk[13][8],Fk[13][9],Fk[13][10],Fk[13][11],
		Fk[14][6],Fk[14][7],Fk[14][8],Fk[14][9],Fk[14][10],Fk[14][11],
		

	    Fk[0][12],Fk[0][13],Fk[0][14],
		Fk[1][12],Fk[1][13],Fk[1][14],
		Fk[2][12],Fk[2][13],Fk[2][14],
		Fk[3][12],Fk[3][13],Fk[3][14],
		Fk[4][12],Fk[4][13],Fk[4][14],
		Fk[5][12],Fk[5][13],Fk[5][14],
		Fk[6][12],Fk[6][13],Fk[6][14],
		Fk[7][12],Fk[7][13],Fk[7][14],
		Fk[8][12],Fk[8][13],Fk[8][14],
		Fk[9][12],Fk[9][13],Fk[9][14],
		Fk[10][12],Fk[10][13],Fk[10][14],
		Fk[11][12],Fk[11][13],Fk[11][14],
		Fk[12][12],Fk[12][13],Fk[12][14],
		Fk[13][12],Fk[13][13],Fk[13][14],
		Fk[14][12],Fk[14][13],Fk[14][14]);

*/



		//均方差一步预测 Pk_k_1 = F*Pk_1*Fk+Qk
		for(i=0; i<STS15; i++)    
		{
			for(j=0; j<STS15; j++)
			{
				pkmid1[i][j] = 0;
				for(k=0; k<STS15; k++)
				{
					pkmid1[i][j] += Fk[i][k] * SystemData.Pk[k][j];
				}
			}
		}

/* printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n"

		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n\n\n"
		
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n",
	    pkmid1[0][0],pkmid1[0][1],pkmid1[0][2],pkmid1[0][3],
		pkmid1[0][4],pkmid1[0][5],pkmid1[1][0],
		pkmid1[1][1],pkmid1[1][2],pkmid1[1][3],pkmid1[1][4],pkmid1[1][5],pkmid1[2][0],
		pkmid1[2][1],pkmid1[2][2],pkmid1[2][3],pkmid1[2][4],pkmid1[2][5],pkmid1[3][0],
		pkmid1[3][1],pkmid1[3][2],pkmid1[3][3],pkmid1[3][4],pkmid1[3][5],pkmid1[4][0],
		pkmid1[4][1],pkmid1[4][2],pkmid1[4][3],pkmid1[4][4],pkmid1[4][5],
		pkmid1[5][0],pkmid1[5][1],pkmid1[5][2],pkmid1[5][3],pkmid1[5][4],pkmid1[5][5],
		pkmid1[6][0],pkmid1[6][1],pkmid1[6][2],pkmid1[6][3],pkmid1[6][4],pkmid1[6][5],
		pkmid1[7][0],pkmid1[7][1],pkmid1[7][2],pkmid1[7][3],pkmid1[7][4],pkmid1[7][5],
		pkmid1[8][5],
		pkmid1[9][0],pkmid1[9][1],pkmid1[9][2],pkmid1[9][3],pkmid1[9][4],pkmid1[9][5],
		pkmid1[10][0],pkmid1[10][1],pkmid1[10][2],pkmid1[10][3],pkmid1[10][4],pkmid1[10][5],
		pkmid1[11][0],pkmid1[11][1],pkmid1[11][2],pkmid1[11][3],pkmid1[11][4],pkmid1[11][5],
		pkmid1[12][0],pkmid1[12][1],pkmid1[12][2],pkmid1[12][3],pkmid1[12][4],pkmid1[12][5],
		pkmid1[13][0],pkmid1[13][1],pkmid1[13][2],pkmid1[13][3],pkmid1[13][4],pkmid1[13][5],
		pkmid1[14][0],pkmid1[14][1],pkmid1[14][2],pkmid1[14][3],pkmid1[14][4],pkmid1[14][5],
		
		pkmid1[0][6],pkmid1[0][7],pkmid1[0][8],pkmid1[0][9],pkmid1[0][10],pkmid1[0][11],
		pkmid1[1][6],pkmid1[1][7],pkmid1[1][8],pkmid1[1][9],pkmid1[1][10],pkmid1[1][11],
		pkmid1[2][6],pkmid1[2][7],pkmid1[2][8],pkmid1[2][9],pkmid1[2][10],pkmid1[2][11],
		pkmid1[3][6],pkmid1[3][7],pkmid1[3][8],pkmid1[3][9],pkmid1[3][10],pkmid1[3][11],
		pkmid1[4][6],pkmid1[4][7],pkmid1[4][8],pkmid1[4][9],pkmid1[4][10],pkmid1[4][11],
		pkmid1[5][6],pkmid1[5][7],pkmid1[5][8],pkmid1[5][9],pkmid1[5][10],pkmid1[5][11],
		pkmid1[6][6],pkmid1[6][7],pkmid1[6][8],pkmid1[6][9],pkmid1[6][10],pkmid1[6][11],
		pkmid1[7][6],pkmid1[7][7],pkmid1[7][8],pkmid1[7][9],pkmid1[7][10],pkmid1[7][11],
		pkmid1[8][6],pkmid1[8][7],pkmid1[8][8],pkmid1[8][9],pkmid1[8][10],pkmid1[8][11],
		pkmid1[9][6],pkmid1[9][7],pkmid1[9][8],pkmid1[9][9],pkmid1[9][10],pkmid1[9][11],
		pkmid1[10][6],pkmid1[10][7],pkmid1[10][8],pkmid1[10][9],pkmid1[10][10],pkmid1[10][11],
		pkmid1[11][6],pkmid1[11][7],pkmid1[11][8],pkmid1[11][9],pkmid1[11][10],pkmid1[11][11],
		pkmid1[12][6],pkmid1[12][7],pkmid1[12][8],pkmid1[12][9],pkmid1[12][10],pkmid1[12][11],
		pkmid1[13][6],pkmid1[13][7],pkmid1[13][8],pkmid1[13][9],pkmid1[13][10],pkmid1[13][11],
		pkmid1[14][6],pkmid1[14][7],pkmid1[14][8],pkmid1[14][9],pkmid1[14][10],pkmid1[14][11],
		

	    pkmid1[0][12],pkmid1[0][13],pkmid1[0][14],
		pkmid1[1][12],pkmid1[1][13],pkmid1[1][14],
		pkmid1[2][12],pkmid1[2][13],pkmid1[2][14],
		pkmid1[3][12],pkmid1[3][13],pkmid1[3][14],
		pkmid1[4][12],pkmid1[4][13],pkmid1[4][14],
		pkmid1[5][12],pkmid1[5][13],pkmid1[5][14],
		pkmid1[6][12],pkmid1[6][13],pkmid1[6][14],
		pkmid1[7][12],pkmid1[7][13],pkmid1[7][14],
		pkmid1[8][12],pkmid1[8][13],pkmid1[8][14],
		pkmid1[9][12],pkmid1[9][13],pkmid1[9][14],
		pkmid1[10][12],pkmid1[10][13],pkmid1[10][14],
		pkmid1[11][12],pkmid1[11][13],pkmid1[11][14],
		pkmid1[12][12],pkmid1[12][13],pkmid1[12][14],
		pkmid1[13][12],pkmid1[13][13],pkmid1[13][14],
		pkmid1[14][12],pkmid1[14][13],pkmid1[14][14]);
*/
		for(i=0; i<STS15; i++)    
		{
			for(j=0; j<STS15; j++)
			{
				SystemData.Pk[i][j] = 0.0;
				for(k=0; k<STS15; k++)
				{
					SystemData.Pk[i][j] += pkmid1[i][k] * Fk[j][k];
				}
			}
		}
		
		for(i=0; i<STS15; i++)    
		{
			SystemData.Pk[i][i] += SystemData.Qt[i] * TF;
		}


 /*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n"

		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n\n\n"
		
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n",
	    SystemData.Pk[0][0],SystemData.Pk[0][1],SystemData.Pk[0][2],SystemData.Pk[0][3],
		SystemData.Pk[0][4],SystemData.Pk[0][5],SystemData.Pk[1][0],
		SystemData.Pk[1][1],SystemData.Pk[1][2],SystemData.Pk[1][3],SystemData.Pk[1][4],SystemData.Pk[1][5],SystemData.Pk[2][0],
		SystemData.Pk[2][1],SystemData.Pk[2][2],SystemData.Pk[2][3],SystemData.Pk[2][4],SystemData.Pk[2][5],SystemData.Pk[3][0],
		SystemData.Pk[3][1],SystemData.Pk[3][2],SystemData.Pk[3][3],SystemData.Pk[3][4],SystemData.Pk[3][5],SystemData.Pk[4][0],
		SystemData.Pk[4][1],SystemData.Pk[4][2],SystemData.Pk[4][3],SystemData.Pk[4][4],SystemData.Pk[4][5],
		SystemData.Pk[5][0],SystemData.Pk[5][1],SystemData.Pk[5][2],SystemData.Pk[5][3],SystemData.Pk[5][4],SystemData.Pk[5][5],
		SystemData.Pk[6][0],SystemData.Pk[6][1],SystemData.Pk[6][2],SystemData.Pk[6][3],SystemData.Pk[6][4],SystemData.Pk[6][5],
		SystemData.Pk[7][0],SystemData.Pk[7][1],SystemData.Pk[7][2],SystemData.Pk[7][3],SystemData.Pk[7][4],SystemData.Pk[7][5],
		SystemData.Pk[8][5],
		SystemData.Pk[9][0],SystemData.Pk[9][1],SystemData.Pk[9][2],SystemData.Pk[9][3],SystemData.Pk[9][4],SystemData.Pk[9][5],
		SystemData.Pk[10][0],SystemData.Pk[10][1],SystemData.Pk[10][2],SystemData.Pk[10][3],SystemData.Pk[10][4],SystemData.Pk[10][5],
		SystemData.Pk[11][0],SystemData.Pk[11][1],SystemData.Pk[11][2],SystemData.Pk[11][3],SystemData.Pk[11][4],SystemData.Pk[11][5],
		SystemData.Pk[12][0],SystemData.Pk[12][1],SystemData.Pk[12][2],SystemData.Pk[12][3],SystemData.Pk[12][4],SystemData.Pk[12][5],
		SystemData.Pk[13][0],SystemData.Pk[13][1],SystemData.Pk[13][2],SystemData.Pk[13][3],SystemData.Pk[13][4],SystemData.Pk[13][5],
		SystemData.Pk[14][0],SystemData.Pk[14][1],SystemData.Pk[14][2],SystemData.Pk[14][3],SystemData.Pk[14][4],SystemData.Pk[14][5],
		
		SystemData.Pk[0][6],SystemData.Pk[0][7],SystemData.Pk[0][8],SystemData.Pk[0][9],SystemData.Pk[0][10],SystemData.Pk[0][11],
		SystemData.Pk[1][6],SystemData.Pk[1][7],SystemData.Pk[1][8],SystemData.Pk[1][9],SystemData.Pk[1][10],SystemData.Pk[1][11],
		SystemData.Pk[2][6],SystemData.Pk[2][7],SystemData.Pk[2][8],SystemData.Pk[2][9],SystemData.Pk[2][10],SystemData.Pk[2][11],
		SystemData.Pk[3][6],SystemData.Pk[3][7],SystemData.Pk[3][8],SystemData.Pk[3][9],SystemData.Pk[3][10],SystemData.Pk[3][11],
		SystemData.Pk[4][6],SystemData.Pk[4][7],SystemData.Pk[4][8],SystemData.Pk[4][9],SystemData.Pk[4][10],SystemData.Pk[4][11],
		SystemData.Pk[5][6],SystemData.Pk[5][7],SystemData.Pk[5][8],SystemData.Pk[5][9],SystemData.Pk[5][10],SystemData.Pk[5][11],
		SystemData.Pk[6][6],SystemData.Pk[6][7],SystemData.Pk[6][8],SystemData.Pk[6][9],SystemData.Pk[6][10],SystemData.Pk[6][11],
		SystemData.Pk[7][6],SystemData.Pk[7][7],SystemData.Pk[7][8],SystemData.Pk[7][9],SystemData.Pk[7][10],SystemData.Pk[7][11],
		SystemData.Pk[8][6],SystemData.Pk[8][7],SystemData.Pk[8][8],SystemData.Pk[8][9],SystemData.Pk[8][10],SystemData.Pk[8][11],
		SystemData.Pk[9][6],SystemData.Pk[9][7],SystemData.Pk[9][8],SystemData.Pk[9][9],SystemData.Pk[9][10],SystemData.Pk[9][11],
		SystemData.Pk[10][6],SystemData.Pk[10][7],SystemData.Pk[10][8],SystemData.Pk[10][9],SystemData.Pk[10][10],SystemData.Pk[10][11],
		SystemData.Pk[11][6],SystemData.Pk[11][7],SystemData.Pk[11][8],SystemData.Pk[11][9],SystemData.Pk[11][10],SystemData.Pk[11][11],
		SystemData.Pk[12][6],SystemData.Pk[12][7],SystemData.Pk[12][8],SystemData.Pk[12][9],SystemData.Pk[12][10],SystemData.Pk[12][11],
		SystemData.Pk[13][6],SystemData.Pk[13][7],SystemData.Pk[13][8],SystemData.Pk[13][9],SystemData.Pk[13][10],SystemData.Pk[13][11],
		SystemData.Pk[14][6],SystemData.Pk[14][7],SystemData.Pk[14][8],SystemData.Pk[14][9],SystemData.Pk[14][10],SystemData.Pk[14][11],
		

	    SystemData.Pk[0][12],SystemData.Pk[0][13],SystemData.Pk[0][14],
		SystemData.Pk[1][12],SystemData.Pk[1][13],SystemData.Pk[1][14],
		SystemData.Pk[2][12],SystemData.Pk[2][13],SystemData.Pk[2][14],
		SystemData.Pk[3][12],SystemData.Pk[3][13],SystemData.Pk[3][14],
		SystemData.Pk[4][12],SystemData.Pk[4][13],SystemData.Pk[4][14],
		SystemData.Pk[5][12],SystemData.Pk[5][13],SystemData.Pk[5][14],
		SystemData.Pk[6][12],SystemData.Pk[6][13],SystemData.Pk[6][14],
		SystemData.Pk[7][12],SystemData.Pk[7][13],SystemData.Pk[7][14],
		SystemData.Pk[8][12],SystemData.Pk[8][13],SystemData.Pk[8][14],
		SystemData.Pk[9][12],SystemData.Pk[9][13],SystemData.Pk[9][14],
		SystemData.Pk[10][12],SystemData.Pk[10][13],SystemData.Pk[10][14],
		SystemData.Pk[11][12],SystemData.Pk[11][13],SystemData.Pk[11][14],
		SystemData.Pk[12][12],SystemData.Pk[12][13],SystemData.Pk[12][14],
		SystemData.Pk[13][12],SystemData.Pk[13][13],SystemData.Pk[13][14],
		SystemData.Pk[14][12],SystemData.Pk[14][13],SystemData.Pk[14][14]);

*/


		
		//状态矩阵一步预测 Xk_k_1 = F*Xk_1
		for(i=0; i<STS15; i++)
		{
	        Xkk_1[i] = 0;
			for(j=0; j<STS15; j++)
			{
				Xkk_1[i] += Fk[i][j] * SystemData.Xk[j];	       
			}
		}
		
		for(i=0; i<STS15; i++)
		{
			SystemData.Xk[i] = Xkk_1[i];
		}
		
		//状态转移有效标志
		SystemData.TransCount ++;

		//标志清零
		SystemData.TransFlag = 0;   //时间更新完毕后清零，以准备下一次时间更新周期到来
	}
}

/********** 卡尔曼滤波计算 **********/
void KF_Filter(DBL Xk[STS15],        DBL Pk[STS15][STS15],
		   	   DBL Hk[6][STS15],  DBL Zk[6],    DBL R[6][6], 
		   	   UI  FilterTime, 	 	 UI *pMesFlag,    DBL multiple)
{
    int i, j, k;
    DBL Kk[STS15][WAT];
    DBL filter_m1[WAT][STS15];
    DBL filter_m2[WAT][WAT];
    DBL filter_m4[WAT][WAT];
    DBL filter_m5[STS15][WAT];
    DBL pkmid1[STS15][STS15];
	DBL pkmid2[STS15][STS15];
	DBL I0;
	DBL Z_HX[WAT];

	//无量测不滤波
	if(*pMesFlag != 1)
	{
		return;
	}
	*pMesFlag = 0;

	//滤波方差阵(H*Pk_k_1*Ht+R)
	for(i=0; i<WAT; i++)
	{
		for(j=0; j<STS15; j++)
		{
			filter_m1[i][j] = 0.0;
			for(k=0; k<STS15; k++)
			{
			    filter_m1[i][j] += Hk[i][k] * Pk[k][j];
			}
		}
	} 


/*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n"

		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n\n\n"
		
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n %.16f %.16f %.16f\n %.16f %.16f %.16f\n",
	    Pk[0][0],Pk[0][1],Pk[0][2],Pk[0][3],
		Pk[0][4],Pk[0][5],Pk[1][0],
		Pk[1][1],Pk[1][2],Pk[1][3],Pk[1][4],Pk[1][5],Pk[2][0],
		Pk[2][1],Pk[2][2],Pk[2][3],Pk[2][4],Pk[2][5],Pk[3][0],
		Pk[3][1],Pk[3][2],Pk[3][3],Pk[3][4],Pk[3][5],Pk[4][0],
		Pk[4][1],Pk[4][2],Pk[4][3],Pk[4][4],Pk[4][5],
		Pk[5][0],Pk[5][1],Pk[5][2],Pk[5][3],Pk[5][4],Pk[5][5],
		Pk[6][0],Pk[6][1],Pk[6][2],Pk[6][3],Pk[6][4],Pk[6][5],
		Pk[7][0],Pk[7][1],Pk[7][2],Pk[7][3],Pk[7][4],Pk[7][5],
		Pk[8][5],
		Pk[9][0],Pk[9][1],Pk[9][2],Pk[9][3],Pk[9][4],Pk[9][5],
		Pk[10][0],Pk[10][1],Pk[10][2],Pk[10][3],Pk[10][4],Pk[10][5],
		Pk[11][0],Pk[11][1],Pk[11][2],Pk[11][3],Pk[11][4],Pk[11][5],
		Pk[12][0],Pk[12][1],Pk[12][2],Pk[12][3],Pk[12][4],Pk[12][5],
		Pk[13][0],Pk[13][1],Pk[13][2],Pk[13][3],Pk[13][4],Pk[13][5],
		Pk[14][0],Pk[14][1],Pk[14][2],Pk[14][3],Pk[14][4],Pk[14][5],
		
		Pk[0][6],Pk[0][7],Pk[0][8],Pk[0][9],Pk[0][10],Pk[0][11],
		Pk[1][6],Pk[1][7],Pk[1][8],Pk[1][9],Pk[1][10],Pk[1][11],
		Pk[2][6],Pk[2][7],Pk[2][8],Pk[2][9],Pk[2][10],Pk[2][11],
		Pk[3][6],Pk[3][7],Pk[3][8],Pk[3][9],Pk[3][10],Pk[3][11],
		Pk[4][6],Pk[4][7],Pk[4][8],Pk[4][9],Pk[4][10],Pk[4][11],
		Pk[5][6],Pk[5][7],Pk[5][8],Pk[5][9],Pk[5][10],Pk[5][11],
		Pk[6][6],Pk[6][7],Pk[6][8],Pk[6][9],Pk[6][10],Pk[6][11],
		Pk[7][6],Pk[7][7],Pk[7][8],Pk[7][9],Pk[7][10],Pk[7][11],
		Pk[8][6],Pk[8][7],Pk[8][8],Pk[8][9],Pk[8][10],Pk[8][11],
		Pk[9][6],Pk[9][7],Pk[9][8],Pk[9][9],Pk[9][10],Pk[9][11],
		Pk[10][6],Pk[10][7],Pk[10][8],Pk[10][9],Pk[10][10],Pk[10][11],
		Pk[11][6],Pk[11][7],Pk[11][8],Pk[11][9],Pk[11][10],Pk[11][11],
		Pk[12][6],Pk[12][7],Pk[12][8],Pk[12][9],Pk[12][10],Pk[12][11],
		Pk[13][6],Pk[13][7],Pk[13][8],Pk[13][9],Pk[13][10],Pk[13][11],
		Pk[14][6],Pk[14][7],Pk[14][8],Pk[14][9],Pk[14][10],Pk[14][11],
		

	    Pk[0][12],Pk[0][13],Pk[0][14],
		Pk[1][12],Pk[1][13],Pk[1][14],
		Pk[2][12],Pk[2][13],Pk[2][14],
		Pk[3][12],Pk[3][13],Pk[3][14],
		Pk[4][12],Pk[4][13],Pk[4][14],
		Pk[5][12],Pk[5][13],Pk[5][14],
		Pk[6][12],Pk[6][13],Pk[6][14],
		Pk[7][12],Pk[7][13],Pk[7][14],
		Pk[8][12],Pk[8][13],Pk[8][14],
		Pk[9][12],Pk[9][13],Pk[9][14],
		Pk[10][12],Pk[10][13],Pk[10][14],
		Pk[11][12],Pk[11][13],Pk[11][14],
		Pk[12][12],Pk[12][13],Pk[12][14],
		Pk[13][12],Pk[13][13],Pk[13][14],
		Pk[14][12],Pk[14][13],Pk[14][14]);*/



/*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n\n\n",
	    filter_m1[0][0],filter_m1[0][1],filter_m1[0][2],filter_m1[0][3],
		filter_m1[0][4],filter_m1[0][5],filter_m1[1][0],
		filter_m1[1][1],filter_m1[1][2],filter_m1[1][3],filter_m1[1][4],filter_m1[1][5],filter_m1[2][0],
		filter_m1[2][1],filter_m1[2][2],filter_m1[2][3],filter_m1[2][4],filter_m1[2][5],filter_m1[3][0],
		filter_m1[3][1],filter_m1[3][2],filter_m1[3][3],filter_m1[3][4],filter_m1[3][5],filter_m1[4][0],
		filter_m1[4][1],filter_m1[4][2],filter_m1[4][3],filter_m1[4][4],filter_m1[4][5],
		filter_m1[5][0],filter_m1[5][1],filter_m1[5][2],filter_m1[5][3],filter_m1[5][4],filter_m1[5][5],
		
		filter_m1[0][6],filter_m1[0][7],filter_m1[0][8],filter_m1[0][9],filter_m1[0][10],filter_m1[0][11],
		filter_m1[1][6],filter_m1[1][7],filter_m1[1][8],filter_m1[1][9],filter_m1[1][10],filter_m1[1][11],
		filter_m1[2][6],filter_m1[2][7],filter_m1[2][8],filter_m1[2][9],filter_m1[2][10],filter_m1[2][11],
		filter_m1[3][6],filter_m1[3][7],filter_m1[3][8],filter_m1[3][9],filter_m1[3][10],filter_m1[3][11],
		filter_m1[4][6],filter_m1[4][7],filter_m1[4][8],filter_m1[4][9],filter_m1[4][10],filter_m1[4][11],
		filter_m1[5][6],filter_m1[5][7],filter_m1[5][8],filter_m1[5][9],filter_m1[5][10],filter_m1[5][11],


	    filter_m1[0][12],filter_m1[0][13],filter_m1[0][14],
		filter_m1[1][12],filter_m1[1][13],filter_m1[1][14],
		filter_m1[2][12],filter_m1[2][13],filter_m1[2][14],
		filter_m1[3][12],filter_m1[3][13],filter_m1[3][14],
		filter_m1[4][12],filter_m1[4][13],filter_m1[4][14],
		filter_m1[5][12],filter_m1[5][13],filter_m1[5][14]);
*/




/*printf("%.16f %.16f %.16f\n" 
	    "%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n"
		"%.16f %.16f %.16f\n",
        Hk[0][3],Hk[1][4],Hk[2][5],Hk[3][6],Hk[4][7],Hk[5][8],
	    Hk[0][0],Hk[0][1],Hk[0][2],Hk[0][3],
		Hk[0][4],Hk[0][5],Hk[1][0],
		Hk[1][1],Hk[1][2],Hk[1][3],Hk[1][4],Hk[1][5],Hk[2][0],
		Hk[2][1],Hk[2][2],Hk[2][3],Hk[2][4],Hk[2][5],Hk[3][0],
		Hk[3][1],Hk[3][2],Hk[3][3],Hk[3][4],Hk[3][5],Hk[4][0],
		Hk[4][1],Hk[4][2],Hk[4][3],Hk[4][4],Hk[4][5],
		Hk[5][0],Hk[5][1],Hk[5][2],Hk[5][3],Hk[5][4],Hk[5][5],
		Hk[0][6],Hk[0][7],Hk[0][8],Hk[0][9],Hk[0][10],Hk[0][11],
		Hk[0][12],Hk[0][12],Hk[0][13],Hk[0][14],Hk[1][6],Hk[1][7],
		Hk[1][8]);
*/





	for(i=0; i<WAT; i++)    
	{
		for(j=0; j<WAT; j++)
		{
			filter_m2[i][j] = 0.0;
			for(k=0; k<STS15; k++)
			{
			    filter_m2[i][j] += filter_m1[i][k] * Hk[j][k];
			}
		}
	}         
	for(i=0; i<WAT; i++)    
	{
		for(j=0; j<WAT; j++)
		{
			filter_m4[i][j]  = filter_m2[i][j] + R[i][j];
		}
	}


/*	printf("%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n",
		  
	filter_m4[0][0],filter_m4[0][1],filter_m4[0][2],
	filter_m4[1][0],filter_m4[1][1],filter_m4[1][2],
	filter_m4[2][0],filter_m4[2][1],filter_m4[2][2]);
*/


	//Z-H*X
	for(i=0; i<WAT; i++)
	{
		Z_HX[i] = 0;
		for(j=0; j<STS15; j++)
		{
			Z_HX[i] += Hk[i][j] * Xk[j];
		}

		Z_HX[i] = Zk[i] - Z_HX[i];
	}

/*	printf("%.16f\n %.16f\n %.16f\n",
	Z_HX[0],Z_HX[1],Z_HX[2]);
*/
	
	//滤波校准前判断
	if(FilterTime >= 2)
	{
		for(i=0; i<WAT; i++)
		{
		    if( fabs(Z_HX[i]) > multiple * sqrt(filter_m4[i][i]) )   //检验是否溢出
		    {				
		        return;
		    }
	    }
    }

     
	//求逆阵得m4
	if( Matrix_Inversesix( filter_m4 ) == 0 )    //6维求逆
    {
		return;
    }

    //滤波增益阵 K = Pk_k_1*Ht*inv(H*Pk_k_1*Ht+R)
    for(i=0; i<STS15; i++)    
	{
		for(j=0; j<WAT; j++)
		{
			filter_m5[i][j] = 0.0;
			for(k=0;k<STS15; k++)
			{
			    filter_m5[i][j] += Pk[i][k] * Hk[j][k];
			}
		}
	}


/*	printf("\n%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n",
		  
	filter_m4[0][0],filter_m4[0][1],filter_m4[0][2],
	filter_m4[1][0],filter_m4[1][1],filter_m4[1][2],
	filter_m4[2][0],filter_m4[2][1],filter_m4[2][2]);
*/


/*	printf("%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
           "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n\n\n"
           "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n\n\n\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
           "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
           "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
           "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n",

	Hk[5][0],Hk[5][1],Hk[5][2],Hk[5][3],Hk[5][4],Hk[5][5],Hk[5][6],Hk[5][7],Hk[5][8],Hk[5][9],Hk[5][10],Hk[5][11],
    Hk[5][12],Hk[5][13],Hk[5][14],

	Pk[14][0], Pk[14][1], Pk[14][2], Pk[14][3], Pk[14][4], Pk[14][5], Pk[14][6], Pk[14][7], Pk[14][8], 
	Pk[14][9], Pk[14][10], Pk[14][11], Pk[14][12], Pk[14][13], Pk[14][14], 
	filter_m5[0][0],filter_m5[0][1],filter_m5[0][2],filter_m5[0][3],
	filter_m5[0][4],filter_m5[0][5],
	filter_m5[1][0],
	filter_m5[1][1],filter_m5[1][2],filter_m5[1][3],filter_m5[1][4],filter_m5[1][5],filter_m5[2][0],
	filter_m5[2][1],filter_m5[2][2],filter_m5[2][3],filter_m5[2][4],filter_m5[2][5],filter_m5[3][0],
	filter_m5[3][1],filter_m5[3][2],filter_m5[3][3],filter_m5[3][4],filter_m5[3][5],filter_m5[4][0],
	filter_m5[4][1],filter_m5[4][2],filter_m5[4][3],filter_m5[4][4],filter_m5[4][5],
	filter_m5[5][0],filter_m5[5][1],filter_m5[5][2],filter_m5[5][3],filter_m5[5][4],filter_m5[5][5],
	
	filter_m5[6][0],filter_m5[6][1],filter_m5[6][2],filter_m5[6][3],filter_m5[6][4],filter_m5[6][5],
	filter_m5[7][0],filter_m5[7][1],filter_m5[7][2],filter_m5[7][3],filter_m5[7][4],filter_m5[7][5],
	filter_m5[8][0],filter_m5[8][1],filter_m5[8][2],filter_m5[8][3],filter_m5[8][4],filter_m5[8][5],
	filter_m5[9][0],filter_m5[9][1],filter_m5[9][2],filter_m5[9][3],filter_m5[9][4],filter_m5[9][5],
	filter_m5[10][0],filter_m5[10][1],filter_m5[10][2],filter_m5[10][3],filter_m5[10][4],filter_m5[10][5],
	filter_m5[11][0],filter_m5[11][1],filter_m5[11][2],filter_m5[11][3],filter_m5[11][4],filter_m5[11][5],
	filter_m5[12][0],filter_m5[12][1],filter_m5[12][2],filter_m5[12][3],filter_m5[12][4],filter_m5[12][5],
	filter_m5[13][0],filter_m5[13][1],filter_m5[13][2],filter_m5[13][3],filter_m5[13][4],filter_m5[13][5],
	filter_m5[14][0],filter_m5[14][1],filter_m5[14][2],filter_m5[14][3],filter_m5[14][4],filter_m5[14][5]);

*/

	for(i=0; i<STS15; i++)    
	{
		for(j=0; j<WAT; j++)
		{
			Kk[i][j] = 0;
			for(k=0;k<WAT; k++)
			{
			    Kk[i][j] += filter_m5[i][k] * filter_m4[k][j];
			}
		}
	}


/*	printf("%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n"

		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n "
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n"
		   "%.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n"
		   "%.16f\n %.16f\n %.16f\n %.16f\n %.16f\n",
	Kk[0][0],Kk[0][1],Kk[0][2],Kk[0][3],
	Kk[0][4],Kk[0][5],
	Kk[1][0],
	Kk[1][1],Kk[1][2],Kk[1][3],Kk[1][4],Kk[1][5],Kk[2][0],
	Kk[2][1],Kk[2][2],Kk[2][3],Kk[2][4],Kk[2][5],Kk[3][0],
	Kk[3][1],Kk[3][2],Kk[3][3],Kk[3][4],Kk[3][5],Kk[4][0],
	Kk[4][1],Kk[4][2],Kk[4][3],Kk[4][4],Kk[4][5],
	Kk[5][0],Kk[5][1],Kk[5][2],Kk[5][3],Kk[5][4],Kk[5][5],
	
	Kk[6][0],Kk[6][1],Kk[6][2],Kk[6][3],
	Kk[6][4],Kk[6][5],
	Kk[7][0],
	Kk[7][1],Kk[7][2],Kk[7][3],Kk[7][4],Kk[7][5],Kk[8][0],
	Kk[8][1],Kk[8][2],Kk[8][3],Kk[8][4],Kk[8][5],Kk[9][0],
	Kk[9][1],Kk[9][2],Kk[9][3],Kk[9][4],Kk[9][5],Kk[10][0],
	Kk[10][1],Kk[10][2],Kk[10][3],Kk[10][4],Kk[10][5],
	Kk[11][0],Kk[11][1],Kk[11][2],Kk[11][3],Kk[11][4],Kk[11][5],

	Kk[12][0],Kk[12][1],Kk[12][2],Kk[12][3],
	Kk[12][4],Kk[12][5],
	Kk[13][0],
	Kk[13][1],Kk[13][2],Kk[13][3],Kk[13][4],Kk[13][5],Kk[14][0],
	Kk[14][1],Kk[14][2],Kk[14][3],Kk[14][4],Kk[14][5]);
*/


	//状态估计方程 Xk = Xk_k_1+K*(Z-H*Xk_k_1)
	for(i=0; i<STS15; i++)    
	{
		for(j=0; j<WAT; j++)
		{
			Xk[i] += Kk[i][j] * Z_HX[j];
		}
	}         
    
	//估计均方差 Pk=(I-K*H)*Pk_k_1
	for(i=0; i<STS15; i++)    
	{
		for(j=0; j<STS15; j++)
		{
			pkmid1[i][j] = 0.0;
			for(k=0; k<WAT; k++)
			{
			    pkmid1[i][j] += Kk[i][k] * Hk[k][j];
			}
		}
	}

	for(i=0; i<STS15; i++)    
	{
		for(j=0; j<STS15; j++)
		{
			if(i == j)
			{
				I0 = 1;
			}
			else
			{
				I0 = 0;
			}
			
			pkmid1[i][j] = I0 - pkmid1[i][j];    
		}
	}
	
	for(i=0; i<STS15; i++)    
	{
		for(j=0; j<STS15; j++)
		{
			pkmid2[i][j] = 0.0;
			for(k=0;k<STS15; k++)
			{
			    pkmid2[i][j] += pkmid1[i][k] * Pk[k][j];
			}
		}
	}

	//对称化处理
	for(i=1; i<STS15; i++)
	{
		for(j=0; j<i; j++)
		{
			pkmid2[i][j] = pkmid2[j][i];
		}
	}
	
	for(i=0; i<STS15; i++)
	{
	    for(j=0; j<STS15; j++)
	    {
	        Pk[i][j] = pkmid2[i][j];
	    }
	}

	*pMesFlag = 2;   //滤波成功标志
}


/********** 卡尔曼滤波校正 **********/
void KF_Correct()
{
	int i;
	DBL Fai[3];

	//滤波未结束
	if(SystemData.MesFlag[0] == 1)   //为1表示量测量有效需要进行量测更新；2：量测更新完毕(滤波成功)
		{
			return;
		}

	//滤波解算成功判断
	if(SystemData.MesFlag[0] == 2)
		{
			SystemData.MesFlag[0]      = 0;   //修正完毕后清零
		//	SystemData.NoFilterTime = 0;
			SystemData.FilterCount ++;  //滤波成功次数
		
/*	printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %d \n%",
			   SystemData.Xk[0],SystemData.Xk[1],SystemData.Xk[2],SystemData.Xk[3],SystemData.Xk[4],SystemData.Xk[5],
			   SystemData.Xk[6],SystemData.Xk[7],SystemData.Xk[8],
			   SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],SystemData.MesFlag[0]);
*/
	//滤波解算成功

		//速度修正
		for(i=0; i<3; i++)
		{
			SystemData.OutData.Vg[i] -= SystemData.Xk[i+3];
			SystemData.Xk[i+3] = 0;
		}
		
		//姿态修正
		Fai[0] = SystemData.Xk[0];
		Fai[1] = SystemData.Xk[1];
		Fai[2] = SystemData.Xk[2];
	//	Fai[2] = 0;



		Attitude_Amend(Fai, SystemData.Cb_n[0], SystemData.OutData.Att, SystemData.Q);	  //修正姿态	
		SystemData.Xk[0] = 0;
		SystemData.Xk[1] = 0;
		SystemData.Xk[2] = 0;
		//SystemData.Xk[5] = 0;

		//累加姿态修正量
		for(i=0; i<3; i++)
		{
			SystemData.DFai[i] += Fai[i];
		}

/*		//位置修正
		for(i=0; i<3; i++)
		{
			SystemData.OutData.Pos[i] -= SystemData.Xk[i+6];	
			SystemData.Xk[i+6] = 0;
		} 
*/		
		//SystemData.OutData.Pos[0] = Limit_Angle(SystemData.OutData.Pos[0], 1); 
		//SystemData.OutData.Pos[1] = Limit_Angle(SystemData.OutData.Pos[1], 0);

/*			printf("\n\n%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n%",
			   SystemData.Xk[0],SystemData.Xk[1],SystemData.Xk[2],SystemData.Xk[3],SystemData.Xk[4],SystemData.Xk[5],
			   SystemData.Xk[6],SystemData.Xk[7],SystemData.Xk[8],
			   SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2]);
*/

	}
}

/********** INS/GPS滤波器量测更新 **********/
/*void GPS_KF_MesUpdate(DBL GPS[6])
{
	int i;
		for(i=0; i<3; i++)
		{
			SystemData.Z1[i] = SystemData.OutData.Vg[i] - GPS[i];
            SystemData.Z1[i+3] = SystemData.OutData.Pos[i] - GPS[i+3];

		}	
}
*/



/********** 由机体角增量更新四元数和速度增量 **********/
void Delta_Cal(DBL DTheta[3],   DBL Delta_Vf[3], EARTH_DATA *pEarth,
			   DBL DTheta_b[3], DBL Delta_Vg[3], 
			   DBL Att[3], 	    DBL Cb_n[3][3],  DBL Q[4])
{
	int i, j;
	
	/********************** 速度增量计算 ************************/

	//Delta_Vg = Cb_n * Delta_Vf
	for(i=0; i<3; i++)
    {
    	Delta_Vg[i] = 0;
    	
    	for(j=0; j<3; j++)
    	{
    		Delta_Vg[i] += Cb_n[i][j] * Delta_Vf[j];
    	}
    }

	//重力和有害加速度补偿
	Delta_Vg[0] -=  pEarth->A_Cori[0] * SystemData.Tm;
	Delta_Vg[1] -=  pEarth->A_Cori[1] * SystemData.Tm;
	Delta_Vg[2] -= (pEarth->A_Cori[2] + pEarth->g) * SystemData.Tm;

	/**********************   姿态更新    ***********************/

	//四元数更新
	for(i=0; i<3; i++)
	{
		DTheta_b[i] = DTheta[i] - pEarth->Win_b[i] * SystemData.Tm;
	}
	
	Refresh_Q(DTheta_b, Q);

	//姿态矩阵计算
	Convert_Q_To_Cbn(Q, Cb_n[0]);
	
	//计算姿态角
	Convert_Cbn_To_Att(Cb_n[0], Att);
}

/********** 捷联惯导解算 **********/
void SINS_Update(DBL DTheta[3],  DBL Delta_Vf[3],
			  	 DBL Pos[3],     DBL Vg[3],  DBL Att[3],
			  	 DBL Cb_n[3][3], DBL Q[4],   DBL Dins[3], EARTH_DATA *pEarth)
{
	int i;
	DBL Vg0[3], DTheta_b[3], Delta_Vg[3];

	//计算地球相关参数
	Calculate_Earth(Pos, Vg, Cb_n[0], pEarth);

	//增量计算
	Delta_Cal(DTheta, Delta_Vf, pEarth, DTheta_b, Delta_Vg, Att, Cb_n, Q);
	
	/************************  速度更新 ************************/
	
	//记录前周期速度
	for(i=0; i<3; i++)
	{
		Vg0[i] = Vg[i];
	}

	//速度更新
	for(i=0; i<3; i++)
	{
		Vg[i] += Delta_Vg[i];
	}
	
	/************************  位置更新 ************************/
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
void Navigation(NAVI_OUT *pNaviData, DBL IMU1[6], DBL IMU2[6])
{
	int i;
    DBL DTheta_b[3] = {0.0},Vg0[3]= {0.0},Delta_Vg[3]= {0.0},Delta_Vgend[3]= {0.0};
    DBL Dw1[3] = {0.0,0.0,0.0}, Dw2[3] = {0.0,0.0,0.0};
    DBL Dv1[3] = {0.0,0.0,0.0}, Dv2[3] = {0.0,0.0,0.0};
    DBL Dw[3] = {0.0},DV[3]= {0.0};
    DBL dVrotm[3] = {0.0},dVsculm[3] = {0.0},dvG_Corm[3] = {0.0};
    DBL Temp1[3] = {0.0},Temp2[3] = {0.0},Temp3[3] = {0.0},Temp4[3] = {0.0};
	DBL Dins[3] = {0.0};
	DBL Qtmp[4] = {0.0};
	DBL Cn2n[3][3] = {0.0};
	DBL VWin_n[3] = {0.0};
    DBL Cb2nk_1[3][3] = {0.0};


    for(i=0; i<3; i++)
	{
	   Dw1[i] = IMU1[i];
	   Dw2[i] = IMU2[i];
	   Dv1[i] = IMU1[i+3];
	   Dv2[i] = IMU2[i+3];
	
	}

	//计算地球相关参数
	Calculate_Earth(SystemData.OutData.Pos, SystemData.OutData.Vg, SystemData.Cb_n[0], &Earth);
    Convert_Q_To_Cbn(SystemData.Q, Cb2nk_1[0]);  //前一时刻的姿态用于速度解算
/*	printf("%.16f\n%",
			   Earth.g);
*/
////////////////////////////////////////////////////////////////////////////////
	Cross(Dw1,Dw2,Temp1);  //Temp1 = cross(Dw1,Dw2);

	for(i=0; i<3; i++)
	{
		SystemData.ThetaA2[i] = Dw1[i] + Dw2[i]  + 2.0/3*Temp1[i]; //phi = IMUk.dw1 + IMUk.dw2 + 2/3*cross(IMUk.dw1,IMUk.dw2);
	}


	/**********************   姿态更新    ***********************/
	//四元数更新
	for(i=0; i<3; i++)
	{
		DTheta_b[i] = SystemData.ThetaA2[i] - Earth.Win_b[i] * SystemData.Tm;   //Phi + VMulQuat(Phin2nk_1,Qn2bk_1)
	}
	Refresh_Q(DTheta_b, SystemData.Q);   //
	
	//姿态矩阵计算	
	Convert_Q_To_Cbn(SystemData.Q, SystemData.Cb_n[0]);  //  Cb2n;

	//计算姿态角
	Convert_Cbn_To_Att(SystemData.Cb_n[0], SystemData.OutData.Att);   //Cb_n = Cb2n
	
	/********************** 速度更新******************************/
	for(i=0; i<3; i++)
	{
		Vg0[i] = SystemData.OutData.Vg[i];   //前周期速度值
	}
	
	for(i=0; i<3; i++)
    {
		Dw[i] = Dw1[i] + Dw2[i];   //Theta
        DV[i] = Dv1[i] + Dv2[i];   //dV
    }

    


    Cross(Dw,DV,Temp1);   //旋转效应补偿项计算
    dVrotm[0] = 0.5*Temp1[0];
    dVrotm[1] = 0.5*Temp1[1];
    dVrotm[2] = 0.5*Temp1[2];   //1/2*cross(Theta,dV)
    
    Cross(Dv1,Dw2,Temp2);   //cross(DV{1},DW{2})
    Cross(Dw1,Dv2,Temp3);   //cross(DW{1},DV{2})
	for(i=0; i<3; i++)
	{   
		dVsculm[i] = 2./3*(Temp2[i] + Temp3[i]);
		SystemData.DeltaV2[i] = DV[i] + dVrotm[i] + dVsculm[i] ;  //dV+dVrotm + dVsculm;
		SystemData.DV[i] = DV[i];   //为了组合导航计算fb用
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
	/************************  位置更新 ************************/
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
	*pNaviData = SystemData.OutData;   //返回输出
}



/********** 滤波解算函数 **********/
void Kalman_Filter()
{
	UI Mes;
	//计算F阵
    Cal_KF_FNav(SystemData.OutData.Pos[0], SystemData.OutData.Vg, SystemData.DV, SystemData.Cb_n, Earth);
	
	if( SystemData.TransFlag )//SystemData.TransFlag :状态转移标志(时间跟新标志)每0.1s时间更新周期到了被置1
	{
    	//时间更新
	    KF_TimeUpdate();
        Mes = SystemData.MesFlag[0];     //1:量测有效

	    //INS/GPS滤波解算
	    KF_Filter(SystemData.Xk,  SystemData.Pk, 
			      SystemData.H1,  SystemData.Z1, SystemData.R1, 
			      SystemData.FilterCount, &SystemData.MesFlag[0], 6);

   //    printf("%d\n",
//			   SystemData.MesFlag[0]);

        //滤波结果修正
	    KF_Correct();	
		//溢出计数

    	if( (Mes == 1) && (SystemData.MesFlag == 0) )   //说明溢出了
		{
			SystemData.X2Over_Count ++; //SystemData.X2Over_Count
		}

	}
		
}



/**********反正弦函数**********/
DBL Limit_asin(DBL data)
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

/**********角度合理值转换**********/
DBL Limit_Angle(DBL Ang, int Flag)
{
	DBL data = Ang;
	
	switch( Flag )
	{
	case 0:								//转换为-PI/2 ~ +PI/2
		if(data > PI/2)
		{
			data = PI - data;
		}
		else if(data < -PI/2)
		{
			data = -PI - data;
		}
		break;
		
	case 1:								//转换为-PI ~ +PI
		if(data > PI)
		{
			data -= 2*PI;
		}
		else if(data < -PI)
		{
			data += 2*PI;
		}
		break;
		
	case 2:								//转换为 0 ~ +2*PI
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

/**********向量叉乘**********/
void Cross(DBL A[3], DBL B[3], DBL C[3])
{
	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];
	C[2] = A[0]*B[1] - A[1]*B[0];
}

/**********向量的反对称阵**********/
void LCross(DBL V[3], DBL *pM)
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

/**********向量标准化**********/
void norm(DBL Vec[3])
{
	int i;
	DBL m;

	m = sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);

	for(i=0; i<3; i++)
	{
		Vec[i] /= m;
	}
}

/********** 3*3矩阵转置 **********/
void Matrix33_Tran(DBL *pM)
{
	int i,j;
	DBL Temp[3][3] = {0.0};

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
double VNorm(DBL *Vin)
{
	return sqrt(Vin[0] * Vin[0] + Vin[1] * Vin[1] + Vin[2] * Vin[2]);
}
///////////////////////////////////////////////////

DBL* Rv2Quat(DBL *Rv, DBL *Qbk2bk_1)
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
DBL* VMulf(DBL *Vin, DBL f,DBL *Vout)
{
	Vout[0] = Vin[0]*f;
	Vout[1] = Vin[1]*f;
	Vout[2] = Vin[2]*f;

	return Vout;
}

/********** 3*3矩阵乘法 **********/
void Matrix33_Mult(DBL *pa, DBL *pb, DBL *pc)
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

/********** 3*1矩阵乘法 **********/
void Matrix31_Mult(DBL *pa, DBL b[3], DBL c[3])
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

/**********计算地球重力加速度信息**********/
DBL Calculae_G(DDL Pos[3])
{
	DBL Rs,  g;
	DBL Rs2, Rs4, Lat, Rm, Rmh, h;
	
	Lat = Pos[1];
	h   = Pos[2];
	
	//基本参数计算
	Rs = sin(Lat);
	
	Rs2 = Rs  * Rs;
	Rs4 = Rs2 * Rs2;

	Rm  = Re * (1 - 2*e + 3*e*Rs2);
	Rmh = Rm + h;

	//重力加速度计算
	g = 9.7803267714*(1 + 5.27094e-3*Rs2 + 2.32718e-5*Rs4) - 3.086e-6*h;
	//g = 9.780318 * (1 + 0.0053024 * Rs2) * Rm*Rm /(Rmh*Rmh);

	return g;
}

/**********  计算地球相关参数  **********/
void Calculate_Earth(DDL Pos[3], DBL Vg[3], DBL *pCb_n, EARTH_DATA *pEarth)
{
	//pCb_n   : 表示Cn2b
	int    i, j;
	DBL Rs2, Rs4, Lat, Rm, h;
	DBL Wie_n[3], Wcori[3];
	


	//位置数据
	Lat = Pos[0];   //纬度
	h   = Pos[2];   //高度
	
	//地球参数计算
	pEarth->Rs = sin(Lat);
	pEarth->Rc = cos(Lat);
	
	Rs2 = pEarth->Rs  * pEarth->Rs;
	Rs4 = Rs2 * Rs2;
    
	Rm = Re * (1 - 2*e + 3*e*Rs2);
    pEarth->Rmh = Rm + h;
    pEarth->Rnh = Re * (1 + e*Rs2) + h;

    pEarth->Wn = Wie * pEarth->Rc;
    pEarth->Wz = Wie * pEarth->Rs;
    
 	//重力加速度计算
	pEarth->g = 9.7803267714*(1 + 5.27094e-3*Rs2 + 2.32718e-5*Rs4) - 3.086e-6*h;

	//Wie_n计算
	Wie_n[0] = 0;
	Wie_n[1] = pEarth->Wn;
	Wie_n[2] = pEarth->Wz;
	
	//Wen_n计算
	pEarth->Wen_n[0] = -Vg[1]/pEarth->Rmh;
    pEarth->Wen_n[1] =  Vg[0]/pEarth->Rnh;
    pEarth->Wen_n[2] =  Vg[0]/pEarth->Rnh * tan(Lat);

	//Win_n计算
	for(i=0; i<3; i++)
	{
		pEarth->Win_n[i] = Wie_n[i] + pEarth->Wen_n[i];
    }

    //Win_b计算
    for(i=0; i<3; i++)
	{
        pEarth->Win_b[i] = 0;
	    for(j=0; j<3; j++)
	    {
	       pEarth->Win_b[i] += *(pCb_n+j*3+i) * pEarth->Win_n[j];	     
		}
	}
	//有害加速度计算
	for(i=0; i<3; i++)
	{
	    Wcori[i] = 2*Wie_n[i] + pEarth->Wen_n[i];
    }
    
    Cross(Wcori, Vg, pEarth->A_Cori);

}
/**********由姿态计算姿态矩阵Cnb**********/
void Convert_Att_Cbn(DBL Att[3], DBL *pCbn) 
{
	DBL sinP, cosP, sinR, cosR, sinH, cosH;

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

/**********姿态到四元数Qn2b的转换**********/
void Convert_Att_To_Q(DBL Att[3], DDL q[4])
{
	int	i, sign_q;
	DBL	norm;
    DBL	cosH2, cosP2, cosR2, sinH2, sinP2, sinR2;
    
	//建立四元数值
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

/**********四元数到姿态矩阵的转换**********/
void Convert_Q_To_Cbn(DDL q[4], DBL *pCbn)
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

/**********从姿态矩阵Cb2n中提取姿态**********/
void Convert_Cbn_To_Att(DBL *pCbn, DBL Att[3])
{
	Att[0] = Limit_asin(*(pCbn+7));
    Att[1] = atan2(-(*(pCbn+6)), *(pCbn+8));
    Att[2] = atan2( *(pCbn+1), *(pCbn+4));

    if(Att[2] < 0.)
    {
		Att[2] += 2*PI;   //将航向定在0~360范围
	}
}

/********** 四元数更新 **********/
void Refresh_Q(DBL DTheta_b[3], DDL Q[4])
{
	int i;
	DDL f, norm;//DBL f, norm;
	DDL MQ[4], q1[4], q[4];//DBL MQ[4], q1[4], q[4];
	
	//旋转矢量转换为四元数
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
	
	//四元数乘法Qk+1=Qk×q(h)=M'(q(h))*Qk
	q1[0] = MQ[0]*q[0] - MQ[1]*q[1] - MQ[2]*q[2] - MQ[3]*q[3];
	q1[1] = MQ[1]*q[0] + MQ[0]*q[1] + MQ[3]*q[2] - MQ[2]*q[3];
	q1[2] = MQ[2]*q[0] - MQ[3]*q[1] + MQ[0]*q[2] + MQ[1]*q[3];
	q1[3] = MQ[3]*q[0] + MQ[2]*q[1] - MQ[1]*q[2] + MQ[0]*q[3];
	
	//四元数规范化处理
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

/**********  姿态矩阵修正  **********/
void Attitude_Amend(DBL Fai[3], DBL *pCb_n, DBL Att[3], DDL Q[4])
{
    int i, j, k;
    DBL Fai_e, Fai_n, Fai_z;
	DBL C_g_g[3][3], Cbg[3][3];

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
	
	//修正后重新计算姿态角
    Convert_Cbn_To_Att(pCb_n, Att);
	
	//重新计算四元数
	Convert_Att_To_Q(Att, Q);
}


/**********矩阵求逆(6阶矩阵)**********/
int Matrix_Inversesix(DBL M[6][6])
{
	int i, j, k, i0, j0, p0[6], q0[6];
	DBL w, z;

	//求逆
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




/**********矩阵求逆**********/
int Matrix_Inverse(DBL *pM)
{
	int i, j, k, i0, j0, p0[3], q0[3];
	DBL w, z;

	//求逆
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



/********** 水平位置距离计算 **********/
DBL HPos_Dis(DBL Pos0[3], DBL Pos1[3])
{
	DBL Reh, Rc;
	DBL Dpos[2];
	
	Reh = Re + Pos0[2];
	Rc  = cos(Pos0[1]);
	
	Dpos[0] = Limit_Angle(Pos0[0] - Pos1[0], 1) * Reh * Rc;
	Dpos[1] = (Pos0[1] - Pos1[1]) * Reh;
	
	return sqrt(Dpos[0]*Dpos[0] + Dpos[1]*Dpos[1]);
}


/********** 四元数乘法运算Q3=Q2*Q1 **********/
void Quat_Mul(DBL Q1[4], DBL Q2[4], DBL Q3[4])
{
	DBL QR[4] = {0}, norm = 0;
	//四元数乘法Qk+1=Qk×q(h)=M'(q(h))*Qk
	QR[0] = Q1[0]*Q2[0] - Q1[1]*Q2[1] - Q1[2]*Q2[2] - Q1[3]*Q2[3];
	QR[1] = Q1[1]*Q2[0] + Q1[0]*Q2[1] + Q1[3]*Q2[2] - Q1[2]*Q2[3];
	QR[2] = Q1[2]*Q2[0] - Q1[3]*Q2[1] + Q1[0]*Q2[2] + Q1[1]*Q2[3];
	QR[3] = Q1[3]*Q2[0] + Q1[2]*Q2[1] - Q1[1]*Q2[2] + Q1[0]*Q2[3];
	
	//四元数规范化处理 / norm
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

// 计算三阶方阵的行列式
DBL sdet_33(DBL *pM)
{
	DBL det = 0.0;

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
void QQ2Phi(DBL Qb2nc[4], DBL Qb2n[4], DBL Phi[3])//Phi = Qb2n*Qb2nc'
{
	DBL Qtmp[4];
	DBL dQ[4];
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
//参考秦永元《惯性导航》p306, 9.3.4的逆过程
void Quat2Rv(DBL Qb2n[4], DBL Rv[3])//Quaternion to rv
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
