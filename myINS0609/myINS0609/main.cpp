/************************************************************************/
/* �˲�������������                                                      */
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
	double MINS[9] = {0,0,0, 0,10,0, 0.6108652539,1.8849555922,400};		//��ʼ��Ϣ
	double imu_avp[15] = {0.0};

	fileIMU = fopen("vcdata.txt","r");
	fileINSGPSResult = fopen("INSResult.dat","w");
	FILE *fp_KFRes = fopen("out.dat", "w");

	//��ʼ��
	Navigation_Init(MINS);   //������ʼ��
	InitKalmanFilter();   //�˲�����ʼ��

	while ((resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
		&IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
		&IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
		&IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
		&IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21])) != EOF)
	{
		// ѭ���ж��ж�ȡ��һ֡���ݣ���Ч��ֵ��Ȼ����Ŷ�ȡ�ڶ�֡���ݣ���Ч��ֵ������IMU���ݶ�ȡ����������
		for (i=0;i<6;i++)
			IMU_1[i] = IMU_Read[i];

		resFileRead = fscanf(fileIMU,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
			&IMU_Read[0],&IMU_Read[1],&IMU_Read[2],&IMU_Read[3],&IMU_Read[4],&IMU_Read[5],
			&IMU_Read[6],&IMU_Read[7],&IMU_Read[8],&IMU_Read[9],&IMU_Read[10],&IMU_Read[11],
			&IMU_Read[12],&IMU_Read[13],&IMU_Read[14],&IMU_Read[15],&IMU_Read[16],&IMU_Read[17],
			&IMU_Read[18],&IMU_Read[19],&IMU_Read[20],&IMU_Read[21]);

		if (EOF == resFileRead)
		{
			break;   //���ݶ�ȡ��ϣ��������
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
				MINS[i+3] = IMU_Read[i+6];  //GPS���ݸ���
			matrix_set(kfData.Yk, 0,0, SystemData.OutData.Vg[0] - MINS[3]);
			matrix_set(kfData.Yk, 1,0, SystemData.OutData.Vg[1] - MINS[4]);
			matrix_set(kfData.Yk, 2,0, SystemData.OutData.Vg[2] - MINS[5]);
			matrix_set(kfData.Yk, 3,0, SystemData.OutData.Pos[0] - MINS[6]);
			matrix_set(kfData.Yk, 4,0, SystemData.OutData.Pos[1] - MINS[7]);
			matrix_set(kfData.Yk, 5,0, SystemData.OutData.Pos[2] - MINS[8]);
			kfData.MeasureEnable  = 1;
		}

		//�˲�����
		KalmanFilterStd();

		//�ж��˲��������Ƿ��м����쳣���������
		if (matrix_ErrorNum != 0)
		{
			switch (matrix_ErrorNum)
			{
			case 1:
				printf("���������쳣 #1������ת���������\n");
				break;
			case 2:
				printf("���������쳣 #2��������������Χ��\n");
				break;
			case 3:
				printf("���������쳣 #3����״��ͬ�ľ����ܸ�ֵ����ֵ���Ӽ������У���\n");
				break;
			case 4:
				printf("���������쳣 #4����������˷�����\n");
				break;
			case 5:
				printf("���������쳣 #5������Ϊ���������档\n");
				break;
			default:
				printf("���������쳣 #9��δ֪����\n");
				break;
			}
			printf("��������˳�����\n");
			getchar();
			break;
		}

		//������º�״̬����
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

		if (IMU_Count % 100 == 0)   //1s���ˣ��洢�������
		{
			// ����������
			fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f  %.16f %.16f\n",
				SystemData.OutData.Att[0]*Rad_Deg,SystemData.OutData.Att[1]*Rad_Deg,SystemData.OutData.Att[2]*Rad_Deg,
				SystemData.OutData.Vg[0],SystemData.OutData.Vg[1],SystemData.OutData.Vg[2],
				SystemData.OutData.Pos[0],SystemData.OutData.Pos[1],SystemData.OutData.Pos[2]);
			KalmanFilterOut(fp_KFRes);

			// ���������ʾ
			printf("\b\b\b\b\b\b\b%5.0f s",IMU_Count*SAMPLE_TIME_S);

			if (IMU_Count == 30000)
			{
				break;
			}
		}
	}

	int t1 = clock();
	printf("running time�� %d ms", t1-t0);

	fclose(fileIMU);
	fclose(fileINSGPSResult);
	fclose(fp_KFRes);
	getchar();
	return 0;
}
