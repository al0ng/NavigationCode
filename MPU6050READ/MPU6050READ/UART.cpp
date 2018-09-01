//======================================================================
//
//  UART.cpp  串口操作
//
//======================================================================
#include "UART.h"

//======================================================================
//
//  名称    UART_init
//  描述    串口初始化函数（9600-8-1）
//  输入		LPCTSTR COMx 初始化的串口
//  输出    HANDLE handle  串口操作句柄
//  版本    V1.00
//
HANDLE UART_init(LPCTSTR COMx)
{
	int tag_getState;
	DCB DCB_SPC;
	HANDLE handle;

	handle = CreateFile(COMx, GENERIC_READ|GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
	if (handle == INVALID_HANDLE_VALUE)
	{
		printf("Open COM Failed!!\n\n按任意键退出……");
		getchar();
		exit(0);
	}

	tag_getState = GetCommState(handle, &DCB_SPC);
	if (tag_getState != 0)
	{
		DCB_SPC.BaudRate = CBR_9600;
		//DCB_SPC.fParity = 0;
		DCB_SPC.ByteSize = 8;
		DCB_SPC.Parity = NOPARITY;
		DCB_SPC.StopBits = ONESTOPBIT;
	}
	SetCommState(handle, &DCB_SPC);

	SetupComm(handle, 2048, 1024);				
	PurgeComm(handle, PURGE_TXCLEAR|PURGE_RXCLEAR);

	return handle;
}

//======================================================================
//
//  名称    readTag
//  描述    读取开头标志
//  输入		HANDLE handle 串口操作句柄
//  输出    int 读取完成标志
//  版本    V1.00
//
int readTag(HANDLE handle)
{
	static unsigned char tag_string[4] = {0};
	unsigned char tag;
	DWORD wCount;
	BOOL bReadState;
	bReadState = ReadFile(handle, &tag, 1, &wCount, NULL);
	tag_string[0] = tag_string[1];
	tag_string[1] = tag_string[2];
	tag_string[2] = tag_string[3];
	tag_string[3] = tag;
	if (!bReadState)
	{
		printf("Read UART Failed!\n");
		return 0;
	}
	else if (tag_string[0]==250 && tag_string[1]==255 && tag_string[2]==0x36 && tag_string[3]==0x39)
		return 1;
	else
		return 0;
}

//======================================================================
//
//  名称    readUART
//  描述    读取一帧数据
//  输入		HANDLE handle 串口操作句柄
//			unsigned char* data  缓冲区指针
//			int num  读取数据长度
//  输出    
//  版本    V1.00
//
void readUART(HANDLE handle, unsigned char* data, int num)
{
	DWORD wCount;
	BOOL bReadState;
	bReadState = ReadFile(handle, data, num, &wCount, NULL);
	if (!bReadState)
	{
		printf("Read UART Failed!\n");
	}
}

//======================================================================
//
//  名称    writeUART
//  描述    串口写操作函数
//  输入		HANDLE handle 串口操作句柄
//			unsigned char* data  写入数据指针
//  输出    int 写操作结果
//  版本    V1.00
//
int writeUART(HANDLE handle, unsigned char* data)
{
	DWORD wCount;
	BOOL bWriteState;
	bWriteState = WriteFile(handle, data, sizeof(data), &wCount, NULL);
	if(bWriteState)
		return 1;
	else
		return 0;
}

//======================================================================
//
//  名称    CheckSum
//  描述    数据求和检验
//  输入		unsigned char* data  待检验数据指针
//  输出    unsigned char sum  求和检验结果
//  版本    V1.00
//
unsigned char CheckSum(unsigned char *data)
{
	unsigned char sum = 0x36+0x39;
	for (int i=0; i<57; i++)
	{
		sum += data[i];
	}
// 	if (sum > 255)
// 	{
// 		sum=~sum;
// 		sum+=1;
// 	}
	return sum;
}