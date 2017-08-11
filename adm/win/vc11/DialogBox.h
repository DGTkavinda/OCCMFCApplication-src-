#pragma once


// CDialogBox dialog

class CDialogBox : public CDialog
{
	DECLARE_DYNAMIC(CDialogBox)

public:
	CDialogBox(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDialogBox();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
