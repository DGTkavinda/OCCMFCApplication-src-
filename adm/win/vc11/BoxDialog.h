#pragma once
//#include "resource.h"

// CBoxDialog dialog

class CBoxDialog : public CDialog
{
	DECLARE_DYNAMIC(CBoxDialog)

public:
	CBoxDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~CBoxDialog();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
