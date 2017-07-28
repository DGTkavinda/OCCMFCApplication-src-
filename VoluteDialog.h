#pragma once

#include "res/resource.h"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "Geom_BezierCurve.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepGprop.hxx"
#include "GProp_GProps.hxx"
#include <ColoredShapes.h>
#include <OCC_3dDoc.h>
#include "OCC_BaseDoc.h"
#include <BRepOffsetAPI_ThruSections.hxx>
#include "BRepBuilderAPI_Transform.hxx"
#include <AIS_InteractiveContext.hxx>
#include "gp_Circ.hxx"
// CVoluteDialog dialog

class CVoluteDialog : public CDialog
{
	DECLARE_DYNAMIC(CVoluteDialog)

public:
	CVoluteDialog(Handle_AIS_InteractiveContext,CWnd* pParent = NULL);   // standard constructor
	virtual ~CVoluteDialog();
	TopoDS_Wire createNewShapeAccordingToR1height(double,double,double,double,double);
	TopoDS_Wire createNewShapeWithRightArea(double,double,double,double,double,double);
	double getSurfaceArea(TopoDS_Wire);
	double getTrapezuimHeight(double,double,double,double);
	CEdit m_edtHeight;
	CEdit m_edtWidth;
	CEdit m_edtAngle1;
	CEdit m_edtAngle2;
	CEdit m_edtPercentage;
	CEdit m_edtFraction;
	CEdit m_edtLength;
	CEdit m_edtInlet;
// Dialog Data
	enum { IDD = IDD_DIALOG_Volute };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	CColoredShapes* m_pcoloredshapeList;
	  Handle(AIS_InteractiveContext) myAISContext;

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClicked1();
};
