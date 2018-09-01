/************************************************************************/
/* 滤波函数测试正常                                                      */
/************************************************************************/
#include <stdio.h>
#include "matrix.h"
#include "Navigation.h"
#include "myKalmanFilter.h"

KF_Data kfData;
ALIGN_DATA Align;
extern  SYSTEM_DATA SystemData;
extern int matrix_ErrorNum;
NAVI_OUT   OutputData;

double MINSAtt[3];
double MINSQ[4];
double phi[3];
double Cbm2n[3][3];

int main()
{
	int t0 = clock();
	FILE *fileIMU, *fileGPS, *fileINSGPSResult;
	int IMU_IDX = 0, IMU_cnt = 0,GPS_IDX_1 = 0, GPS_cnt = 0,  GPS_IDX_2 = 0, IMU_Count = 0, resFileRead = 0;
	int i = 0, j = 0, k = 0, GPS_Lost = 0, GPS_num;
	double IMU_1[6] = {0.0}, IMU_2[6] = {0.0}, GPS_1[6] = {0.0}, GPS_2[6] = {0.0}, GPS_Save[6] = {0.0};
	double IMU_Read[22] = {0.0}, Head_c = 0.0, Head_t = 0.0, GPS_second = 0,GPS_GDOP = 0.0;
	double MINS[9] = {0,0,0, 0,10,0, 0.6108652539,1.8849555922,400};		//初始信息
	double imu_avp[15] = {0.0};

	fileIMU = fopen("vcdata.txt","r");
	fileINSGPSResult = fopen("INSResult.dat","w");
	FILE *fp_KFRes = fopen("out.dat", "w");

	//初始化
	Navigation_Init(MINS);   //导航初始化
	InitKalmanFilter();   //滤波器初始化

	while ((resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
		&IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
		&IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
		&IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
		&IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21])) != EOF)
	{
		// 循环判断中读取第一帧数据，有效则赋值；然后接着读取第二帧数据，有效则赋值。更新IMU数据读取个数计数器
		for (i=0;i<6;i++)
			IMU_1[i] = IMU_Read[i];

		resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
			&IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
			&IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
			&IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
			&IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21]);

		if (EOF == resFileRead)
		{
			break;   //数据读取完毕，程序结束
		}
		else
		{
			for (i=0;i<6;i++)
				IMU_2[i] = IMU_Read[i];
			
		}
		IMU_Count += 2;

		Navigation(&OutputData,IMU_1, IMU_2);
		for (int ii=0; ii<3; ii++)
		{
			imu_avp[ii]    = IMU_2[ii];
			imu_avp[ii+3]  = IMU_2[ii+3];
			imu_avp[ii+6]  = OutputData.Att[ii];
			imu_avp[ii+9]  = OutputData.Vg[ii];
			imu_avp[ii+12] = OutputData.Pos[ii];
		}
		KalmanFilterSetPhi(imu_avp);

		if (IMU_Count % 100 == 0)  //1s
		{
/*			double Fai[3] = {0.0};*/
			for (i=0;i<6;i++)
				MINS[i+3] = IMU_Read[i+6];  //GPS数据更新
			matrix_set(kfData.Yk, 0,0, SystemData.OutData.Vg[0] - MINS[3]);
			matrix_set(kfData.Yk, 1,0, SystemData.OutData.Vg[1] - MINS[4]);
			matrix_set(kfData.Yk, 2,0, SystemData.OutData.Vg[2] - MINS[5]);
			matrix_set(kfData.Yk, 3,0, SystemData.OutData.Pos[0] - MINS[6]);
			matrix_set(kfData.Yk, 4,0, SystemData.OutData.Pos[1] - MINS[7]);
			matrix_set(kfData.Yk, 5,0, SystemData.OutData.Pos[2] - MINS[8]);
			kfData.MeasureEnable  = 1;
		}

		//滤波解算
		KalmanFilterStd();

		//判断滤波过程中是否有计算异常情况并处理
		if (matrix_ErrorNum != 0)
		{
			switch (matrix_ErrorNum)
			{
			case 1:
				printf("矩阵运算异常 #1：矩阵转置运算错误。\n");
				break;
			case 2:
				printf("矩阵运算异常 #2：索引超出矩阵范围。\n");
				break;
			case 3:
				printf("矩阵运算异常 #3：形状不同的矩阵不能赋值（赋值、加减运算中）。\n");
				break;
			case 4:
				printf("矩阵运算异常 #4：矩阵不满足乘法规则。\n");
				break;
			case 5:
				printf("矩阵运算异常 #5：矩阵不为方阵不能求逆。\n");
				break;
			default:
				printf("矩阵运算异常 #9：未知错误。\n");
				break;
			}
			printf("按任意键退出……\n");
			getchar();
			break;
		}

		//量测更新后状态反馈
		if (kfData.MeasureUdtFinished == 1)
		{
			kfData.MeasureUdtFinished = 0;
			for(i=0; i<3; i++)
			{
				SystemData.OutData.Vg[i]  -= matrix_get(kfData.Xk,i+3,0);
				SystemData.OutData.Pos[i] -= matrix_get(kfData.Xk,i+6,0);
				matrix_set(kfData.Xk,i+3,0, 0);
				matrix_set(kfData.Xk,i+6,0, 0);
			}
		}

		if (IMU_Count % 100 == 0)   //1s到了，存储导航结果
		{
			// 计算结果保存
			fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f  %.16f %.16f\n",
				SystemData.OutData.Att[0]*Rad_Deg,SystemData.OutData.Att[1]*Rad_Deg,SystemData.OutData.Att[2]*Rad_Deg,
				SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],
				SystemData.OutData.Pos[0],SystemData.OutData.Pos[1],SystemData.OutData.Pos[2]);
			KalmanFilterOut(fp_KFRes);

			// 程序进程显示
			printf("\b\b\b\b\b\b\b%5.0f s",IMU_Count*SAMPLE_TIME_S);

			if (IMU_Count == 30000)
			{
				break;
			}
		}
	}

	int t1 = clock();
	printf("running time： %d ms", t1-t0);

	fclose(fileIMU);
	fclose(fileINSGPSResult);
	fclose(fp_KFRes);
	getchar();
	return 0;
}
