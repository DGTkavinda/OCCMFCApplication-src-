// DialogBox.cpp : implementation file
//

#include "stdafx.h"
#include "DialogBox.h"
#include "afxdialogex.h"


// CDialogBox dialog

IMPLEMENT_DYNAMIC(CDialogBox, CDialog)

CDialogBox::CDialogBox(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogBox::IDD, pParent)
{

}

CDialogBox::~CDialogBox()
{
}

void CDialogBox::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CDialogBox, CDialog)
END_MESSAGE_MAP()


// CDialogBox message handlers
