#include "MainDlg.h"
#include <string.h>
#include <stdio.h>
#include <Windows.h>
#include <commctrl.h>
#include "NavMain.h"
#include "myKalmanFilter.h"

HANDLE thread;
HWND hListview;
DWORD WINAPI Thread_INSGPS(LPVOID lpvoid);

BOOL WINAPI Main_Proc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch(uMsg)
	{
		HANDLE_MSG(hWnd, WM_INITDIALOG, Main_OnInitDialog);
		HANDLE_MSG(hWnd, WM_COMMAND, Main_OnCommand);
		HANDLE_MSG(hWnd, WM_CLOSE, Main_OnClose);
	}

	return FALSE;
}

BOOL Main_OnInitDialog(HWND hwnd, HWND hwndFocus, LPARAM lParam)
{
	hListview = GetDlgItem(hwnd, IDC_LIST);  
	LVCOLUMN vcl;  
	vcl.mask = LVCF_TEXT | LVCF_WIDTH | LVCF_SUBITEM| LVCF_FMT;  
	vcl.fmt = LVCFMT_RIGHT;
	vcl.pszText = L"num";//列标题  
	vcl.cx = 50;//列宽  
	vcl.iSubItem = 0;//子项索引，第一列无子项  
	ListView_InsertColumn(hListview, 0, &vcl);   
	// 第一列  
	vcl.fmt = LVCFMT_RIGHT;
	vcl.pszText = L"phi_x           ";//列标题  
	vcl.cx = 150;//列宽  
	vcl.iSubItem = 1;//子项索引，第一列无子项  
	ListView_InsertColumn(hListview, 1, &vcl);  
	// 第二列  
	vcl.fmt = LVCFMT_RIGHT;
	vcl.pszText = L"phi_y           ";  
	vcl.cx = 150;  
	vcl.iSubItem = 2;//子项索引  
	ListView_InsertColumn(hListview, 2, &vcl);  
	// 第三列  
	vcl.fmt = LVCFMT_RIGHT;
	vcl.pszText = L"phi_z           ";  
	vcl.cx = 150;  
	vcl.iSubItem = 3;  
	ListView_InsertColumn(hListview, 3, &vcl);  
	return TRUE;
}

void Main_OnCommand(HWND hwnd, int id, HWND hwndCtl, UINT codeNotify)
{
	switch(id)
	{
	case IDC_OK:
		{
			HWND shwnd = hwnd;
			thread = CreateThread(NULL,0,Thread_INSGPS, &shwnd, 0,0);
		}
		break;
	case IDC_STOP:
		{
			if (thread)
			{
				TerminateThread(thread,NULL);
			}
		}
		break;
	case ID_MENU_OPEN:
		{
			if (thread)
			{
				TerminateThread(thread,NULL);
			}
		}
		break;
	default:
		break;
	}
}

void Main_OnClose(HWND hwnd)
{
	CloseHandle(thread);
	EndDialog(hwnd, 0);
}

DWORD WINAPI Thread_INSGPS(LPVOID lpvoid)
{
	HWND hwnd = *(HWND *)lpvoid;
	HWND handle_processBar = GetDlgItem(hwnd, IDC_PROGRESS1);
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
	WCHAR is[15][1024];

	fileIMU = fopen("vcdata.txt","r");
	fileINSGPSResult = fopen("INSResult.dat","w");
	FILE *fp_KFRes = fopen("out.dat", "w");
	FILE *fp_debug = fopen("DebugPhi.dat", "w");

	SendMessage(handle_processBar, PBM_SETRANGE, 0, MAKELPARAM(0, 30000/2));
	SendMessage(handle_processBar, PBM_SETSTEP, (WPARAM)1,0);
	SendMessage(handle_processBar, PBM_SETPOS, (WPARAM)0,0);
	//惯导解算初始化
	insinitial(avp0, nSimple);		//导航初始化
	ins2avp();
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
		SendMessage(handle_processBar, PBM_STEPIT, 0, 0);

		if (IMU_Count % 100 == 0)   //1s到了，存储导航结果
		{
			LVITEM vitem; 
			WCHAR listtext[512];
			vitem.mask = LVIF_TEXT;
			vitem.iItem = IMU_Count/100-1;

			swprintf(listtext, L"%5d", IMU_Count/100);
			vitem.pszText = listtext;
			vitem.iSubItem = 0;
			ListView_InsertItem(hListview,&vitem);

			swprintf(listtext, L"%.16lf", matrix_get(kfData.Xk,0,0));
			vitem.pszText = listtext;
			vitem.iSubItem = 1;
			ListView_SetItem(hListview,&vitem);

			swprintf(listtext, L"%.16lf", matrix_get(kfData.Xk,1,0));
			vitem.pszText = listtext;
			vitem.iSubItem = 2;
			ListView_SetItem(hListview,&vitem);

			swprintf(listtext, L"%.16lf", matrix_get(kfData.Xk,2,0));
			vitem.pszText = listtext;
			vitem.iSubItem = 3;
			ListView_SetItem(hListview,&vitem);
			// 计算结果保存
			fprintf(fileINSGPSResult,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f  %.16f %.16f %.3f\n",
				avp.att[0]*Rad_Deg,avp.att[1]*Rad_Deg,avp.att[2]*Rad_Deg,
				avp.vn[0],avp.vn[1],avp.vn[2],
				avp.pos[0],avp.pos[1],avp.pos[2], IMU_Count*ins.ts);
			/*KalmanFilterOut(fp_KFRes);*/
			for (int i=0; i<15; i++)
			{
				swprintf(is[i], L"%.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf ",
					matrix_get(kfData.Phi,0,i),matrix_get(kfData.Phi,1,i),matrix_get(kfData.Phi,2,i),matrix_get(kfData.Phi,3,i),matrix_get(kfData.Phi,4,i),
					matrix_get(kfData.Phi,5,i),matrix_get(kfData.Phi,6,i),matrix_get(kfData.Phi,7,i),matrix_get(kfData.Phi,8,i),matrix_get(kfData.Phi,9,i),
					matrix_get(kfData.Phi,10,i),matrix_get(kfData.Phi,11,i),matrix_get(kfData.Phi,12,i),matrix_get(kfData.Phi,13,i),matrix_get(kfData.Phi,14,i));
				SetDlgItemText(hwnd, IDC_STATIC0+i, (is[i]));
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