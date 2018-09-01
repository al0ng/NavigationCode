#include <Windows.h>
#include <stdio.h>

HANDLE UART_init(LPCTSTR COMx);
int readTag(HANDLE);
void readUART(HANDLE, unsigned char*, int);
int writeUART(HANDLE, unsigned char*);
unsigned char CheckSum(unsigned char *);