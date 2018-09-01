#include <stdio.h>
#include <string.h>
#include "Navigation.h"

ALIGN_DATA Align;
extern  SYSTEM_DATA SystemData;
NAVI_OUT   OutputData;

DBL MINSAtt[3];
DBL MINSQ[4];
DBL phi[3];
DBL Cbm2n[3][3];


/********* 标定子函数声明 **********************************************/
void Calibration(DBL IMU1[6], DBL IMU2[6]);
/********* 标定子函数声明 **********************************************/

/********* 对准结构体初始化函数声明 ************************************/
void Initial(ALIGN_DATA *Align, DBL Pos_0[6]);
/********* 对准结构体初始化函数声明 ************************************/

/********* 对准结构体更新函数声明 **************************************/
void Align_IMU_Update(ALIGN_DATA *Align, DBL IMU1[6], DBL IMU2[6]);
/********* 对准结构体更新函数声明 **************************************/

/********* 对准结构体更新函数声明 **************************************/
void Align_GPS_Update(ALIGN_DATA *Align, DBL GPS_1[6], DBL GPS_2[6]);
/********* 对准结构体更新函数声明 **************************************/

void main()
{
	int Aligntime = 30;
    unsigned int AlignInitFlag = 1;
    unsigned int INSGPSFlag = 0;
    unsigned int Alignflag = 1;
	FILE *fileIMU, *fileGPS, *fileResult,*fileINSGPSResult;
	int IMU_IDX = 0, IMU_cnt = 0,GPS_IDX_1 = 0, GPS_cnt = 0,  GPS_IDX_2 = 0, IMU_Count = 0, resFileRead = 0;
	int i = 0, j = 0, k = 0, GPS_Lost = 0, GPS_num;
	DBL IMU_1[6] = {0.0}, IMU_2[6] = {0.0}, GPS_1[6] = {0.0}, GPS_2[6] = {0.0}, GPS_Save[6] = {0.0};
	DBL IMU_Read[25] = {0.0}, Head_c = 0.0, Head_t = 0.0, GPS_second = 0,GPS_GDOP = 0.0;
//	DBL Euler[3] = {-0.002407087820183,0.003164411518702,-2.529863184465090};
    DBL MINS[9] = {0.0};
	DBL Cbn[3][3]={0.0}, dL[3]={0.0, 0.0, 0.0}, res_cross[3]={0.0}, wib_b[3]={0.0}, dV[3] = {0.0}; 
	DBL Csm[3][3] = {1.0, 0.0, 0.0,   0.0, 1.0, 0.0,   0.0, 0.0, 1.0};
	DBL Cmn[3][3] = {0.0}, tmp[3][3] = {0.0};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	fileIMU = fopen("Tr_RLG.txt","r");
	fileResult  = fopen("Align_R.dat","w");
	fileINSGPSResult = fopen("INSResult.dat","w");
	//TA = 0.01;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i=0;i<1;i++)
	{
		fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			   &IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
			   &IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
			   &IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
			   &IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21],&IMU_Read[22],&IMU_Read[23],&IMU_Read[24]);
	}
    
//	Initial(&Align, GPS_1);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("对准解算进度：%5.0f s",0);
	while ((resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			   &IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
			   &IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
			   &IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
			   &IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21],&IMU_Read[22],&IMU_Read[23],&IMU_Read[24])) != EOF)
	{

	// 循环判断中读取第一帧数据，有效则赋值；然后接着读取第二帧数据，有效则赋值。更新IMU数据读取个数计数器
		for (i=0;i<6;i++)
		{
			IMU_1[i] = IMU_Read[i+16];
		}

		for (i=0;i<9;i++)
		{
			MINS[i] = IMU_Read[i];  //主惯导姿态、速度、位置
		}

	/**********姿态以及杆臂效应修正***************/
		Convert_Att_Cbn(MINS, Cmn[0]);
		Matrix33_Mult(Cmn[0], Csm[0], tmp[0]);
		Convert_Cbn_To_Att(tmp[0], MINS);

		for (i=0;i<3;i++)
		{
			wib_b[i] = IMU_2[i]/400;
		}
		Convert_Att_Cbn(MINS, Cbn[0]);
		Cross(wib_b, dL, res_cross);
		Matrix31_Mult(Cbn[0], res_cross, dV);
		for (i=0;i<3;i++)
		{
			MINS[i+3] += dV[i];
		}
	/*********************************************/

       //初始化
		if (AlignInitFlag)
		{
			Navigation_Init(MINS);   //导航初始化
			Init_KF_Filter();   //滤波器初始化
			AlignInitFlag = 0;   //初始化只进行一次 
		}
		

		resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			   &IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
			   &IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
			   &IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
			   &IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21],&IMU_Read[22],&IMU_Read[23],&IMU_Read[24]);

		if (EOF == resFileRead)
		{
			break;   //数据读取完毕，程序结束
		}
		else
		{
			for (i=0;i<6;i++)
			{
				IMU_2[i] = IMU_Read[i+16];
			}

			for (i=0;i<9;i++)
			{
				MINS[i] = IMU_Read[i];  //主惯导姿态、速度、位置
			}

		/**********姿态以及杆臂效应修正***************/
			Convert_Att_Cbn(MINS, Cmn[0]);
			Matrix33_Mult(Cmn[0], Csm[0], tmp[0]);
			Convert_Cbn_To_Att(tmp[0], MINS);

			for (i=0;i<3;i++)
			{
				wib_b[i] = IMU_2[i]/400;
			}
			Convert_Att_Cbn(MINS, Cbn[0]);
			Cross(wib_b, dL, res_cross);
			Matrix31_Mult(Cbn[0], res_cross, dV);
			for (i=0;i<3;i++)
			{
				MINS[i+3] += dV[i];
			}
		/*********************************************/
		}
		
		IMU_Count += 2;

/*	  printf("%lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n",IMU_1[0],IMU_1[1],IMU_1[2],IMU_1[3],IMU_1[4],IMU_1[5],
                IMU_2[0],IMU_2[1],IMU_2[2],IMU_2[3],IMU_2[4],IMU_2[5]);
	break;   */

/*			printf("%lf %lf %lf\n %lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf\n",
			 
			 SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],
			 MINS[3],MINS[4],MINS[5],
			 SystemData.OutData.Att[0],SystemData.OutData.Att[1],SystemData.OutData.Att[2],
			 SystemData.OutData.Pos[0],SystemData.OutData.Pos[1],SystemData.OutData.Pos[2],
             MINS[6],MINS[7],MINS[8],
			 MINS[0],MINS[1],MINS[2]);
			   break;   */

					   Navigation(&OutputData,IMU_1, IMU_2);
					   
/*			printf("%lf %lf %lf\n %lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf\n",
			 
			 SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],
			 MINS[3],MINS[4],MINS[5],
			 SystemData.OutData.Att[0],SystemData.OutData.Att[1],SystemData.OutData.Att[2],
			 SystemData.OutData.Pos[0],SystemData.OutData.Pos[1],SystemData.OutData.Pos[2],
             MINS[6],MINS[7],MINS[8],
			 MINS[0],MINS[1],MINS[2]);
			   break;*/

					   SystemData.T100MS_Flag = (IMU_Count % 320 == 0);   //0.8s
					   
					 //  printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",IMU_1[0],IMU_1[1],IMU_1[2],IMU_1[3],IMU_1[4],IMU_1[5],
                     //     IMU_2[0],IMU_2[1],IMU_2[2],IMU_2[3],IMU_2[4],IMU_2[5]);
					//	   break;


                       if (IMU_Count % 320 == 0)  //0.8s
					   {	        

                           SystemData.MesFlag[0] = 1;   //量测标志位

						   MINSAtt[0] = MINS[0];
                           MINSAtt[1] = MINS[1];
						   MINSAtt[2] = MINS[2];

						   Convert_Att_To_Q(MINSAtt,MINSQ);
                           QQ2Phi(SystemData.Q,MINSQ,phi);
						   
						   SystemData.Z1[0] = phi[0];
                           SystemData.Z1[1] = phi[1];
						   SystemData.Z1[2] = phi[2];

			               SystemData.Z1[3] = SystemData.OutData.Vg[0] - MINS[3];  //量测标志有效，计算量测量
						   SystemData.Z1[4] = SystemData.OutData.Vg[1] - MINS[4];
                           SystemData.Z1[5] = SystemData.OutData.Vg[2] - MINS[5];
                  
						   //注意：加入到烧写代码中
						   Convert_Att_Cbn(MINSAtt,Cbm2n[0]);
						   SystemData.H1[0][12] = -Cbm2n[0][0];
						   SystemData.H1[0][13] = -Cbm2n[0][1];
						   SystemData.H1[0][14] = -Cbm2n[0][2];

						   SystemData.H1[1][12] = -Cbm2n[1][0];
						   SystemData.H1[1][13] = -Cbm2n[1][1];
						   SystemData.H1[1][14] = -Cbm2n[1][2];

						   SystemData.H1[2][12] = -Cbm2n[2][0];
						   SystemData.H1[2][13] = -Cbm2n[2][1];
						   SystemData.H1[2][14] = -Cbm2n[2][2];

						   //


						   
						 //  SystemData.Z1[8] = SystemData.OutData.Pos[2] - MINS[8];


  /*     printf("%d %d\n",
			      SystemData.T100MS_Flag,SystemData.MesFlag[0]);

         break;*/
 /*		  printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
			      "%lf %lf %lf\n"
				  "%lf %lf %lf\n"
				  "%lf %lf %lf\n",
			 
			 SystemData.Q[0],SystemData.Q[1],SystemData.Q[2],SystemData.Q[3],
			 MINSQ[0],MINSQ[1],MINSQ[2],MINSQ[3],
			 SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],
			 MINS[3],MINS[4],MINS[5],

			  SystemData.Z1[0],SystemData.Z1[1],SystemData.Z1[2],SystemData.Z1[3],SystemData.Z1[4],SystemData.Z1[5],
             SystemData.Z1[6],SystemData.Z1[7],SystemData.Z1[8],
			 SystemData.H1[0][0],SystemData.H1[0][1],SystemData.H1[0][2],
			 SystemData.H1[0][12],SystemData.H1[0][13],SystemData.H1[0][14],
			 SystemData.H1[1][0],SystemData.H1[1][1],SystemData.H1[1][2]);
			   break;   */

						   
					   }
	
		               //滤波解算
		               Kalman_Filter();
					   
/*if (IMU_Count % 320 == 0) 
{
	break;
}*/

					  if (IMU_Count % 400 == 0)   //1s到了，存储导航结果
					  {
									          // 计算结果保存
			             fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f "
					                 "%.16f  %.16f %.16f %.16f  %.16f %.16f %.16f  %.16f %.16f %.16f  %.16f %.16f\n",
				                OutputData.Att[0]*Rad_Deg,OutputData.Att[1]*Rad_Deg,OutputData.Att[2]*Rad_Deg,
				                OutputData.Vg[0],OutputData.Vg[1],OutputData.Vg[2],
					            OutputData.Pos[0],OutputData.Pos[1],OutputData.Pos[2],
								SystemData.Xk[6],SystemData.Xk[7],SystemData.Xk[8],
								SystemData.Xk[9],SystemData.Xk[10],SystemData.Xk[11],
								SystemData.Xk[12],SystemData.Xk[13],SystemData.Xk[14]);
                       
/*			printf("%lf %lf %lf\n %lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf \n%lf %lf %lf\n",
			 
			 OutputData.Vg[0],OutputData.Vg[1],OutputData.Vg[2],
			 MINS[3],MINS[4],MINS[5],
			 OutputData.Att[0],OutputData.Att[1],OutputData.Att[2],
			 OutputData.Pos[0],OutputData.Pos[1],OutputData.Pos[2],
             MINS[6],MINS[7],MINS[8],
			 MINS[0],MINS[1],MINS[2]);
			   break;*/
						//  break;
		 


						        // 程序进程显示
				                printf("\b\b\b\b\b\b\b%5.0f s",IMU_Count*SAMPLE_TIME_S);

								if (IMU_Count == 270000)
								{
			 	                          break;
								}
					  }
      }

	fclose(fileIMU);
	//fclose(fileGPS);
	fclose(fileResult);
    fclose(fileINSGPSResult);
}

/////////////////////////////////////////////////////////////////////////
/********* 标定子函数定义 **********************************************/
/*void Calibration(DBL IMU1[], DBL IMU2[])
{
	DBL Temp1[6] = {0.0}, Temp2[6] = {0.0};
	int i = 0;
	for (i=0;i<3;i++)
	{
		Temp1[i] = IMU1[i] + 100*Dph_Rps*0.01 + (double)rand()*sqrt(0.01)*0.1*dp05h;
        Temp1[i+3] = IMU1[i+3] + 10*mg*0.01 + (double)rand()*sqrt(0.01)*1*mg;
		
		Temp2[i] = IMU2[i] + 100*Dph_Rps*0.01 + (double)rand()*sqrt(0.01)*0.1*dp05h;
        Temp2[i+3] = IMU2[i+3] + 10*mg*0.01 + (double)rand()*sqrt(0.01)*1*mg;
	}
	for (i=0;i<6;i++)
	{
		IMU1[i] = Temp1[i];
        IMU2[i] = Temp2[i];
	}

}*/

/********* 标定子函数定义 **********************************************/


/********* 对准结构体初始化函数定义 ************************************/
void Initial(ALIGN_DATA *Align, DBL Pos_0[6])
{
	int i = 0, j = 0;

	for (i=0;i<3;i++)
	{
		Align->Qib2b[i+1] = 0.0;
		Align->Qin2n[i+1] = 0.0;
		Align->Qin2ib[i+1] = 0.0;
		Align->Qn2b[i+1] = 0.0;

		Align->Euler_c[i] = 0.0;
		Align->Euler_t[i] = 0.0;

		Align->Vib_f[i] = 0.0;
		Align->Vin_g[i] = 0.0;

		Align->mk = 0.0;
		Align->Vx[i] = 0.0;
		Align->Vy[i] = 0.0;
		Align->Vz[i] = 0.0;

		Align->VPos_0[i] = Pos_0[i];
		Align->VPos_0[i+3] = Pos_0[i+3];

		for (j=0;j<3;j++)
		{
			Align->Vb[i][j] = 0.0;
			Align->Vn[i][j] = 0.0;
		}

	}

	Align->Qib2b[0] = 1.0;
	Align->Qin2n[0] = 1.0;
	Align->Qin2ib[0] = 1.0;
	Align->Qn2b[0] = 1.0;

	Align->Align_Count = 0;
	Align->Ref_T = 1.0;

	Align->cL0 = cos(Align->VPos_0[3]);
	Align->sL0 = sin(Align->VPos_0[3]);

	Align->winie[0] = 0.0;
	Align->winie[1] = Wie*Align->cL0;
	Align->winie[2] = Wie*Align->sL0;

}

/********* 对准结构体初始化函数定义 ************************************/

/********* 对准结构体更新函数声明 **************************************/

void Align_IMU_Update(ALIGN_DATA *Align, DBL IMU1[6], DBL IMU2[6])
{
	DBL Dv1[3] = {IMU1[3],IMU1[4],IMU1[5]}, Dv2[3] = {IMU2[3],IMU2[4],IMU2[5]};
	DBL Dw1[3] = {IMU1[0],IMU1[1],IMU1[2]}, Dw2[3] = {IMU2[0],IMU2[1],IMU2[2]};

	DBL Dv[3] = {Dv1[0]+Dv2[0], Dv1[1]+Dv2[1], Dv1[2]+Dv2[2]};
	DBL Dw[3] = {Dw1[0]+Dw2[0], Dw1[1]+Dw2[1], Dw1[2]+Dw2[2]};

	DBL Temp1[3] = {0.0}, Temp2[3] = {0.0}, Temp3[3] = {0.0}, Temp4[3] = {0.0};

	DBL Cb2ib[3][3] = {0.0}, dvf_b[3] = {0.0}, dvf_ib[3] = {0.0};

	DBL Phi[3] = {0.0}, Norm_Phi = 0.0, Cof = 0.0, Qtemp[4] = {1.0, 0.0, 0.0, 0.0};

	int i = 0, j = 0;

// 惯性系比力积分及姿态跟踪二子样更新计算
	
	Cross(Dw,Dv,Temp1);
	Cross(Dv1,Dw2,Temp2);
	Cross(Dw1,Dv2,Temp3);
	Cross(Dw1,Dw2,Temp4);

	for (i=0;i<3;i++)
	{
		dvf_b[i] = Dv[i] + 0.5*Temp1[i] + 2.0/3.0 * ( Temp2[i] + Temp3[i] );

		Phi[i] = Dw[i] + 2.0/3.0 * Temp4[i];
	}

	Convert_Q_To_Cbn(Align->Qib2b,Cb2ib[0]);

	Matrix31_Mult(Cb2ib[0],dvf_b,dvf_ib);

	Align->Vib_f[0] += dvf_ib[0];
	Align->Vib_f[1] += dvf_ib[1];
	Align->Vib_f[2] += dvf_ib[2];

	Refresh_Q(Phi, Align->Qib2b);

}

/********* 对准结构体更新函数声明 **************************************/

/********* 对准结构体更新函数声明 **************************************/
void Align_GPS_Update(ALIGN_DATA *Align, DBL GPS_1[6], DBL GPS_2[6])
{
	DBL dlm = 0.0, s_dlm = 0.0, c_dlm = 0.0;
	DBL sLt_1 = sin(GPS_1[3]), cLt_1 = cos(GPS_1[3]);
	DBL sLt_2 = sin(GPS_2[3]), cLt_2 = cos(GPS_2[3]);
	
	DBL Cnt2in_t[3][3] = {0.0}, Cnt2in_0[3][3] = {0.0};

	DBL dvrk1[3] = {0.0}, dvrk2[3] = {0.0}, dvrk3[3] = {0.0};

	DBL Temp_1[3] = {0.0}, Temp_2[3] = {0.0}, Temp_3[3] = {0.0};

	DBL vn_1[3] = {GPS_1[0],GPS_1[1],GPS_1[2]};
	DBL vn_2[3] = {GPS_2[0],GPS_2[1],GPS_2[2]};

	DBL wnin[3] = { -(vn_1[1]+vn_2[1])/2.0/Re, 
					(cLt_1+cLt_2)/2.0*Wie + (vn_1[0]+vn_2[0])/2.0/Re, 
					(sLt_1+sLt_2)/2.0*Wie + (vn_1[0]+vn_2[0])/2.0/Re*sLt_2/cLt_2 };

	DBL gn[3] = {0.0,0.0,0.0}, Qtemp[4] = {1.0, 0.0, 0.0, 0.0}, Qtemp_1[4] = {1.0, 0.0, 0.0, 0.0};

	DBL GWM[3][3] = {0.0}, ak = 0.0, Norm_Vf = 0.0, Norm_Vg = 0.0;

	DBL M1_1[3][3] = {0.0}, M1_2[3][3] = {0.0};
	DBL C1[3][3] = {{1,0,0},{0,1,0},{0,0,1}}, C2[3][3] = {{-1,0,0},{0,-1,0},{0,0,1}};

	int i = 0, j = 0;


	gn[2] = -( g0*(1 - 2*GPS_1[5]/Re)*(1 + 5.2884e-3*sLt_1*sLt_1 - 5.9e-6*sin(2*GPS_1[3])*sin(2*GPS_1[3])) + 
			   g0*(1 - 2*GPS_2[5]/Re)*(1 + 5.2884e-3*sLt_2*sLt_2 - 5.9e-6*sin(2*GPS_2[3])*sin(2*GPS_2[3])) ) / 2.0;
	

	Align->Align_Count++;
	
	dlm = GPS_2[4] - Align->VPos_0[4] + Wie * Align->Align_Count * Align->Ref_T;
	s_dlm = sin(dlm);	c_dlm = cos(dlm);

	Cnt2in_t[0][0] =  c_dlm; 
	Cnt2in_t[0][1] = -sLt_2*s_dlm; 
	Cnt2in_t[0][2] =  cLt_2*s_dlm; 
	Cnt2in_t[1][0] =  Align->sL0*s_dlm;
	Cnt2in_t[1][1] =  Align->sL0*sLt_2*c_dlm + Align->cL0*cLt_2;
	Cnt2in_t[1][2] = -Align->sL0*cLt_2*c_dlm + Align->cL0*sLt_2;
	Cnt2in_t[2][0] = -Align->cL0*s_dlm;
	Cnt2in_t[2][1] = -Align->cL0*sLt_2*c_dlm + Align->sL0*cLt_2;
	Cnt2in_t[2][2] =  Align->cL0*cLt_2*c_dlm + Align->sL0*sLt_2;


	Convert_Q_To_Cbn(Align->Qin2n,Cnt2in_0[0]);


	Matrix31_Mult(Cnt2in_0[0],vn_1,Temp_1);
	Matrix31_Mult(Cnt2in_t[0],vn_2,Temp_2);
	for (i=0;i<3;i++)
	{
		dvrk1[i] = Temp_2[i] - Temp_1[i];
	}
	
	
	LCross(wnin,GWM[0]);
	Matrix31_Mult(GWM[0],vn_1,Temp_1);
	Matrix31_Mult(GWM[0],vn_2,Temp_2);
	for (i=0;i<3;i++)
	{
		Temp_3[i] = (vn_1[i]+vn_2[i])*Align->Ref_T/2 + (Temp_1[i]/6 + Temp_2[i]/3)*Align->Ref_T*Align->Ref_T;
	}
	Matrix31_Mult(Cnt2in_0[0],Temp_3,dvrk2);

	
	Matrix31_Mult(GWM[0],gn,Temp_1);
	for (i=0;i<3;i++)
	{
		Temp_2[i] = -gn[i]*Align->Ref_T - Temp_1[i]*Align->Ref_T*Align->Ref_T/2;
	}
	Matrix31_Mult(Cnt2in_0[0],Temp_2,dvrk3);
	

	Cross(Align->winie,dvrk2,Temp_3);
	
	for (i=0;i<3;i++)
	{
		Align->Vin_g[i] += dvrk1[i] + Temp_3[i] + dvrk3[i];
	}
//	printf("\n%.16f %.16f %.16f\n",Align->Vin_g[0],Align->Vin_g[1],Align->Vin_g[2]);


	Convert_Cbn_To_Att(Cnt2in_t[0],Temp_1);
	Temp_1[2] = Limit_Angle(Temp_1[2],1);
//	printf("\n%.16f %.16f %.16f\n",Temp_1[0],Temp_1[1],Temp_1[2]-2*PI);

	Convert_Att_To_Q(Temp_1,Align->Qin2n);

//////////////////////////////////////////////////////////////////////////////////////////////////////

// 以上比力矢量在导航惯性系内积分结果计算完毕，下面完成对准所需递推更新计算
/**/	
	ak = Align->Vib_f[0]*Align->Vib_f[0] + Align->Vib_f[1]*Align->Vib_f[1] + 
		 Align->Vib_f[2]*Align->Vib_f[2];   //V_b0模值
	
	Norm_Vf = sqrt(ak);
	Norm_Vg = sqrt( Align->Vin_g[0]*Align->Vin_g[0] + Align->Vin_g[1]*Align->Vin_g[1] + 
			        Align->Vin_g[2]*Align->Vin_g[2] );

	dvrk2[0] = Align->Vin_g[0] / Norm_Vg;   //V_no单位化  
	dvrk2[1] = Align->Vin_g[1] / Norm_Vg;
	dvrk2[2] = Align->Vin_g[2] / Norm_Vg;

	for (i=0;i<3;i++)
	{
		dvrk1[i] = Align->Vib_f[i] / Norm_Vf;   //V_b0单位化		

		Temp_1[i] = Align->Vx[i] * Align->mk / ( Align->mk + ak ) + 
					dvrk2[0] * dvrk1[i] * ak / ( Align->mk + ak );

		Temp_2[i] = Align->Vy[i] * Align->mk / ( Align->mk + ak ) + 
					dvrk2[1] * dvrk1[i] * ak / ( Align->mk + ak );

		Temp_3[i] = Align->Vz[i] * Align->mk / ( Align->mk + ak ) + 
					dvrk2[2] * dvrk1[i] * ak / ( Align->mk + ak );

		Align->Vx[i] = Temp_1[i];
		Align->Vy[i] = Temp_2[i];
		Align->Vz[i] = Temp_3[i];  //SVV
	}


	LCross(dvrk1,Cnt2in_0[0]);  //Cnt2in_0 = LCross(vb0k)
	Matrix33_Mult(Cnt2in_0[0],Cnt2in_0[0],GWM[0]);   //GWM = (LCross(vb0k))^2

	LCross(dvrk2,Cnt2in_t[0]);  //Cnt2in_t = LCross(vn0k)
	Matrix33_Mult(Cnt2in_t[0],Cnt2in_t[0],Cnt2in_0[0]);   //Cnt2in_0 = (LCross(vn0k))^2

	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)   //mk/(mk + atk)*SVb0 + atk/(mk + atk)*(LCross(vb0k))^2
		{
			Cnt2in_t[i][j] = Align->Vb[i][j] * Align->mk / ( Align->mk + ak ) +    //SVb0 = Align->Vb
							 GWM[i][j] * ak / ( Align->mk + ak );
			Align->Vb[i][j] = Cnt2in_t[i][j];  //保存为下一时刻

			Cnt2in_t[i][j] = Align->Vn[i][j] * Align->mk / ( Align->mk + ak ) + 
							 Cnt2in_0[i][j] * ak / ( Align->mk + ak );   //mk/(mk + atk)*SVn0 + atk/(mk + atk)*(LCross(vn0k))^2

			Align->Vn[i][j] = Cnt2in_t[i][j];
		}
	}


	Align->mk += ak;

/////////////////////////////////////////////////////////////////////////////////////////

// 下面根据递推更新结果完成对准，计算姿态信息
/**/
	if (Align->Align_Count >= 2)   //M1  
	{
		M1_1[0][0] = -Align->Vz[2] - Align->Vy[1];
		M1_1[0][1] =  Align->Vx[1]; M1_1[0][2] = Align->Vx[2];   //SVV = Align->Vx  Align->Vy   Align->Vz

		M1_1[1][0] =  Align->Vy[0]; M1_1[1][2] = Align->Vy[2];
		M1_1[1][1] = -Align->Vz[2] - Align->Vx[0];

		M1_1[2][0] =  Align->Vz[0]; M1_1[2][1] = Align->Vz[1];
		M1_1[2][2] = -Align->Vx[0] - Align->Vy[1];


		M1_2[0][0] = -Align->Vz[2] + Align->Vy[1];
		M1_2[0][1] = -Align->Vx[1];
		M1_2[0][2] =  Align->Vx[2];

		M1_2[1][0] = -Align->Vy[0];
		M1_2[1][1] = -Align->Vz[2] + Align->Vx[0];
		M1_2[1][2] =  Align->Vy[2];

		M1_2[2][0] = -Align->Vz[0];
		M1_2[2][1] = -Align->Vz[1];
		M1_2[2][2] =  Align->Vx[0] + Align->Vy[1];


		
		Matrix33_Mult(C2[0],Align->Vb[0],Cnt2in_0[0]);   //SVb0 = Align->Vb
		Matrix33_Mult(Cnt2in_0[0],C2[0],GWM[0]);

		for (i=0;i<3;i++)
		{
			for (j=0;j<3;j++)
			{
				Cnt2in_0[i][j] = Align->Vb[i][j] + Align->Vn[i][j] + M1_1[i][j] + M1_1[j][i];//Cnt2in_0 = SVb0 + SVn0 + M1 + M1'
				Cnt2in_t[i][j] = GWM[i][j] + Align->Vn[i][j] + M1_2[i][j] + M1_2[j][i];
			}
		}

/*
		printf("%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n",
			   Cnt2in_0[0][0],Cnt2in_0[0][1],Cnt2in_0[0][2],
			   Cnt2in_0[1][0],Cnt2in_0[1][1],Cnt2in_0[1][2],
			   Cnt2in_0[2][0],Cnt2in_0[2][1],Cnt2in_0[2][2] );

		printf("%d %16.16f %16.16f\n",Align->Align_Count,sdet_33(Cnt2in_0[0]),sdet_33(Cnt2in_t[0]));
*/

		if (sdet_33(Cnt2in_0[0])*sdet_33(Cnt2in_0[0]) >= sdet_33(Cnt2in_t[0])*sdet_33(Cnt2in_t[0])) 
		{
			Temp_1[0] = -(-Align->Vy[2] + Align->Vz[1]) * 2.0;
			Temp_1[1] = -(-Align->Vz[0] + Align->Vx[2]) * 2.0;
			Temp_1[2] = -(-Align->Vx[1] + Align->Vy[0]) * 2.0;

			Matrix_Inverse(Cnt2in_0[0]);

			Matrix31_Mult(Cnt2in_0[0],Temp_1,gn);

		}
		else
		{
			Temp_1[0] = -(-Align->Vy[2] - Align->Vz[1]) * 2.0;
			Temp_1[1] = -( Align->Vz[0] + Align->Vx[2]) * 2.0;
			Temp_1[2] = -( Align->Vx[1] - Align->Vy[0]) * 2.0;    //Temp_1 = Dk

			Matrix_Inverse(Cnt2in_t[0]);   //Cnt2in_t = Sk^-1

			Matrix31_Mult(Cnt2in_t[0],Temp_1,gn);   //lk = gn = Sk^-1*Dk

			C1[0][0] = -1.0; C1[1][1] = -1.0;
		}

		//printf("%16.16f %16.16f %16.16f \n",gn[0],gn[1],gn[2]);

		LCross(gn,GWM[0]);   //GWM = (lX)
		GWM[0][0] += 1.0; GWM[1][1] += 1.0; GWM[2][2] += 1.0;    //GWM = (I + (lx))
		Matrix_Inverse(GWM[0]);   //GWM = (I + (lx))^-1
		Temp_1[0] = gn[0] * (-1.0);    
		Temp_1[1] = gn[1] * (-1.0);
		Temp_1[2] = gn[2] * (-1.0);   //Temp_1 = -l
		LCross(Temp_1,C2[0]);    //C2 = (-lX) = -(lx)
		C2[0][0] += 1.0; C2[1][1] += 1.0; C2[2][2] += 1.0;   //c2 = (I -(lx))

		Matrix33_Mult(GWM[0],C2[0],Cnt2in_0[0]);   //Cnt2in_0 = (I + (lx))^-1*(I -(lx))
		Matrix33_Mult(C1[0],Cnt2in_0[0],Cnt2in_t[0]);

/*
		printf("%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n"
			   "\n%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n"
			   "%12.12f %12.12f %12.12f \n",
			   Cnt2in_0[0][0],Cnt2in_0[0][1],Cnt2in_0[0][2],
			   Cnt2in_0[1][0],Cnt2in_0[1][1],Cnt2in_0[1][2],
			   Cnt2in_0[2][0],Cnt2in_0[2][1],Cnt2in_0[2][2],
			   Cnt2in_t[0][0],Cnt2in_t[0][1],Cnt2in_t[0][2],
			   Cnt2in_t[1][0],Cnt2in_t[1][1],Cnt2in_t[1][2],
			   Cnt2in_t[2][0],Cnt2in_t[2][1],Cnt2in_t[2][2]);
*/
		// Cnt2in_t即是所估计的常值姿态阵Cin2ib
		Matrix33_Tran(Cnt2in_t[0]);
		Convert_Cbn_To_Att(Cnt2in_t[0],Align->Euler_c);
		Align->Euler_c[2] = Limit_Angle(Align->Euler_c[2],1);
		Convert_Att_To_Q(Align->Euler_c,Align->Qin2ib);

		
		
		Qtemp[0] = Align->Qin2n[0];
		Qtemp[1] = -Align->Qin2n[1];
		Qtemp[2] = -Align->Qin2n[2];
		Qtemp[3] = -Align->Qin2n[3];

		Quat_Mul(Align->Qin2ib,Qtemp,Qtemp_1);
		Quat_Mul(Align->Qib2b,Qtemp_1,Align->Qn2b);

		Convert_Q_To_Cbn(Align->Qn2b,Cnt2in_t[0]);
		Convert_Cbn_To_Att(Cnt2in_t[0],Align->Euler_t);
		Align->Euler_t[2] = Limit_Angle(Align->Euler_t[2],1);

/*
		printf("%16.16f %16.16f %16.16f \n"
			   "%16.16f %16.16f %16.16f %16.16f \n",
			   Align->Euler_t[0],Align->Euler_t[1],Align->Euler_t[2],
			   Align->Qn2b[0],Align->Qn2b[1],Align->Qn2b[2],Align->Qn2b[3]);
*/

	}


}
/********* 对准结构体更新函数声明 **************************************/