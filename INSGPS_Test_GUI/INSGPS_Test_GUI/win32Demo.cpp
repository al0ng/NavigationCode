#include "win32Demo.h"
#include "MainDlg.h"
#include <CommCtrl.h>
#pragma comment(lib,"ComCtl32.lib")

int WINAPI wWinMain(HINSTANCE hInstance,
                     HINSTANCE hPrevInstance,
                     LPTSTR    lpCmdLine,
                     int       nCmdShow)
{
	InitCommonControls();
	DialogBox(hInstance, MAKEINTRESOURCE(IDD_MAIN), NULL, Main_Proc);
	return 0;
}
