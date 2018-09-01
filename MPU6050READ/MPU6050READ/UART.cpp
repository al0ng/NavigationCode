//======================================================================
//
//  UART.cpp  ���ڲ���
//
//======================================================================
#include "UART.h"

//======================================================================
//
//  ����    UART_init
//  ����    ���ڳ�ʼ��������9600-8-1��
//  ����		LPCTSTR COMx ��ʼ���Ĵ���
//  ���    HANDLE handle  ���ڲ������
//  �汾    V1.00
//
HANDLE UART_init(LPCTSTR COMx)
{
	int tag_getState;
	DCB DCB_SPC;
	HANDLE handle;

	handle = CreateFile(COMx, GENERIC_READ|GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
	if (handle == INVALID_HANDLE_VALUE)
	{
		printf("Open COM Failed!!\n\n��������˳�����");
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
//  ����    readTag
//  ����    ��ȡ��ͷ��־
//  ����		HANDLE handle ���ڲ������
//  ���    int ��ȡ��ɱ�־
//  �汾    V1.00
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
//  ����    readUART
//  ����    ��ȡһ֡����
//  ����		HANDLE handle ���ڲ������
//			unsigned char* data  ������ָ��
//			int num  ��ȡ���ݳ���
//  ���    
//  �汾    V1.00
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
//  ����    writeUART
//  ����    ����д��������
//  ����		HANDLE handle ���ڲ������
//			unsigned char* data  д������ָ��
//  ���    int д�������
//  �汾    V1.00
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
//  ����    CheckSum
//  ����    ������ͼ���
//  ����		unsigned char* data  ����������ָ��
//  ���    unsigned char sum  ��ͼ�����
//  �汾    V1.00
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