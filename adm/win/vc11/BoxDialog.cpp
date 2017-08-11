// BoxDialog.cpp : implementation file
//

#include "stdafx.h"
#include "BoxDialog.h"
#include "afxdialog.h"


// CBoxDialog dialog

IMPLEMENT_DYNAMIC(CBoxDialog, CDialog)

CBoxDialog::CBoxDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CBoxDialog::IDD, pParent)
{

}

CBoxDialog::~CBoxDialog()
{
}

void CBoxDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CBoxDialog, CDialog)
END_MESSAGE_MAP()
	

// CBoxDialog message handlers
