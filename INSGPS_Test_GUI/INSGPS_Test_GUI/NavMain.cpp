/************************************************************************/
/* 滤波函数测试正常                                                      */
/************************************************************************/
#include <stdio.h>
#include "matrix.h"
#include "myNavigation.h"
#include "myKalmanFilter.h"


int navmain()
{
	FILE *fileIMU, *fileINSGPSResult;
	int IMU_Count = 0, resFileRead = 0;
	int i = 0, j = 0, k = 0;
	double IMU_1[6] = {0.0}, IMU_2[6] = {0.0};
	double IMU_Read[22] = {0.0};
	double imu_avp[15] = {0.0};
	extern insdata ins;
	extern avpdata avp;
	extern KF_Data kfData;
	extern int matrix_ErrorNum;
	matrix *wm = matrix_calloc(3,nSimple);
	matrix *vm = matrix_calloc(3,nSimple);
	int ii = 0;

	fileIMU = fopen("vcdata.txt","r");
	fileINSGPSResult = fopen("INSResult.dat","w");
	FILE *fp_KFRes = fopen("out.dat", "w");
	FILE *fp_debug = fopen("DebugPhi.dat", "w");

	//惯导解算初始化
	insinitial(avp0, nSimple);		//导航初始化
	InitKalmanFilter();				//滤波器初始化

	while ((resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
		&IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
		&IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
		&IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
		&IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21])) != EOF)
	{
		ii++;
		for (i=0;i<6;i++)
			IMU_1[i] = IMU_Read[i];

		IMU_Count += 1;

		matrix_set(wm,0,ii-1, IMU_1[0]);
		matrix_set(wm,1,ii-1, IMU_1[1]);
		matrix_set(wm,2,ii-1, IMU_1[2]);
		matrix_set(vm,0,ii-1, IMU_1[3]);
		matrix_set(vm,1,ii-1, IMU_1[4]);
		matrix_set(vm,2,ii-1, IMU_1[5]);

		//nSimple次数据进行一次惯导解算
		if (ii<nSimple)
			continue;
		else
			ii = 0;

		insupdate(wm,vm);
		for (int ii=0; ii<3; ii++)
		{
			imu_avp[ii]    = IMU_1[ii];
			imu_avp[ii+3]  = IMU_1[ii+3];
			imu_avp[ii+6]  = avp.att[ii];
			imu_avp[ii+9]  = avp.vn[ii];
			imu_avp[ii+12] = avp.pos[ii];
		}
		KalmanFilterSetPhi(imu_avp);

		if (IMU_Count % 100 == 0)  //1s
		{
			matrix* Vb = matrix_calloc(3,1);
			matrix* Vn = matrix_calloc(3,1);

			matrix_set(kfData.Yk, 0,0, matrix_get(ins.vn,0,0) - IMU_Read[6]);
			matrix_set(kfData.Yk, 1,0, matrix_get(ins.vn,1,0) - IMU_Read[7]);
			matrix_set(kfData.Yk, 2,0, matrix_get(ins.vn,2,0) - IMU_Read[8]);
			matrix_set(kfData.Yk, 3,0, matrix_get(ins.pos,0,0) - IMU_Read[9]);
			matrix_set(kfData.Yk, 4,0, matrix_get(ins.pos,1,0) - IMU_Read[10]);
			matrix_set(kfData.Yk, 5,0, matrix_get(ins.pos,2,0) - IMU_Read[11]);
			matrix_free(Vb);
			matrix_free(Vn);
			kfData.MeasureEnable  = 1;
		}

		//滤波解算
		KalmanFilterStd();
		KalmanFilterOut(fp_KFRes);

		//量测更新后状态反馈
		if (kfData.MeasureUdtFinished == 1)
		{
			kfData.MeasureUdtFinished = 0;
			for(i=0; i<3; i++)
			{
				matrix_sub(ins.vn, kfData.dvn);
				matrix_sub(ins.pos, kfData.dpos);
				matrix_set(kfData.Xk,i+3,0, 0);
				matrix_set(kfData.Xk,i+6,0, 0);
			}
		}

		ins2avp();

		if (IMU_Count % 100 == 0)   //1s到了，存储导航结果
		{
			// 计算结果保存
			fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f  %.16f %.16f\n",
				avp.att[0]*Rad_Deg,avp.att[1]*Rad_Deg,avp.att[2]*Rad_Deg,
				avp.vn[0],avp.vn[1],avp.vn[2],
				avp.pos[0],avp.pos[1],avp.pos[2]);
			/*KalmanFilterOut(fp_KFRes);*/

			// 程序进程显示
			printf("\b\b\b\b\b\b\b%5.0f s",IMU_Count*ins.ts);

			if (IMU_Count == 30000)
			{
				break;
			}
		}
	}

	fclose(fileIMU);
	fclose(fileINSGPSResult);
	fclose(fp_KFRes);
	fclose(fp_debug);
	getchar();
	return 0;
}
