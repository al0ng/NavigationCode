#include <stdio.h>
#include <Windows.h>
#include "UART.h"

int main(){
	int i = 1000;
	unsigned char buffer[32];
	/*HANDLE handle = UART_init(L"COM11");*/
	FILE* fp_data=	fopen("E:\\WorkSpace\\Visual Studio 2010\\Projects\\MPU6050READ\\MPU6050READ\\data.txt", "r");
	FILE* fp_imu=	fopen("imu.dat", "w");
	while(i--){
		fscanf(fp_data, "%X %X %X %X %X %X %X %X %X %X %X %X ", buffer, buffer+1, buffer+2, buffer+3, buffer+4, buffer+5,
			buffer+6, buffer+7, buffer+8, buffer+9, buffer+10, buffer+11);
		short acc[3], gro[3];
		acc[0] = *(short *)(buffer+0);
		acc[1] = *(short *)(buffer+2);
		acc[2] = *(short *)(buffer+4);
		gro[0] = *(short *)(buffer+6);
		gro[1] = *(short *)(buffer+8);
		gro[2] = *(short *)(buffer+10);
		fprintf(fp_imu, "%d %d %d    %d %d %d\n", acc[0], acc[1], acc[2], gro[0], gro[1], gro[2]);
	}
	
// 	while (i--)
// 	{
// 		/*readUART(handle, buffer, 6);*/
// 		printf("%d\n", i);
// 		printf("%d %d %d   %d %d %d", *(short*)buffer, *(short*)(buffer+2), *(short*)(buffer+4), *(short*)(buffer+6), *(short*)(buffer+8), *(short*)(buffer+10));
// 	}

	return 0;
}