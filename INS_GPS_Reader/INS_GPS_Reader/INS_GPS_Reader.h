
// INS_GPS_Reader.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CINS_GPS_ReaderApp:
// �йش����ʵ�֣������ INS_GPS_Reader.cpp
//

class CINS_GPS_ReaderApp : public CWinApp
{
public:
	CINS_GPS_ReaderApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CINS_GPS_ReaderApp theApp;