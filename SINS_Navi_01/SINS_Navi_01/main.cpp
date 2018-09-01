/************************************************************************/
/* 惯导解算测试程序                                                      */
/************************************************************************/
#include <stdio.h>
#include "matrix.h"
#include "myNavigation.h"
#include "setting.h"

int main()
{
  	FILE *fileIMU, *fileINSGPSResult;
  	int IMU_Count = 0, resFileRead = 0;
  	int i = 0, j = 0, k = 0;
  	double IMU_1[6] = {0.0}, IMU_2[6] = {0.0};
  	double IMU_Read[22] = {0.0};
	extern avpdata avp;
	matrix *wm = matrix_calloc(3,nSimple);
	matrix *vm = matrix_calloc(3,nSimple);
	int ii = 0;
  
  	fileIMU = fopen("vcdata.txt","r");
  	fileINSGPSResult = fopen("INSResult.dat","w");
    
  	//惯导解算初始化
  	insinitial(avp0, nSimple);   //导航初始化
  
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

  		if (IMU_Count % 100 == 0)   //1s到了，存储导航结果
  		{
			ins2avp();
  			// 计算结果保存
   			fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f  %.16f %.16f\n",
   				avp.att[0]*Rad_Deg,avp.att[1]*Rad_Deg,avp.att[2]*Rad_Deg,avp.vn[0],avp.vn[1],avp.vn[2],avp.pos[0],avp.pos[1],avp.pos[2]);
   			// 程序进程显示
   			printf("\b\b\b\b\b\b\b%5.0f s",IMU_Count*avp0[9]);
   			if (IMU_Count == 30000)
   			{
   				break;
   			}
  		}
  	}
  
  	fclose(fileIMU);
  	fclose(fileINSGPSResult);
	printf("\n解算完毕，按任意键退出……");
  	getchar();
	return 0;
}


