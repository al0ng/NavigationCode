
// INS_GPS_ReaderDlg.h : 头文件
//

#pragma once


// CINS_GPS_ReaderDlg 对话框
class CINS_GPS_ReaderDlg : public CDialogEx
{
// 构造
public:
	CINS_GPS_ReaderDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_INS_GPS_READER_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
};
