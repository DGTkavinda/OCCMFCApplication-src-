#pragma once
#include "res/resource.h"
#include <AIS_InteractiveContext.hxx>
#include <ColoredShapes.h>

#include "res/resource.h"
#include "Geom_BezierCurve.hxx"
#include <BRepGProp.hxx>
#include "GProp_GProps.hxx"
#include <BRepOffsetAPI_ThruSections.hxx>
#include "gp_Trsf.hxx"
#include "VoluteDialog.h"
#include "GeomAbs_JoinType.hxx"
#include "BRepOffsetAPI_MakeOffsetShape.hxx"
#include "BRepOffset_Mode.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include <fstream>  
#include <iostream>  
#include <string>  
#include <Geom_BSplineCurve.hxx>
#include <BRepFill_Filling.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepExtrema_DistanceSS.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomLib_Tool.hxx>
#include "GC_MakeArcOfEllipse.hxx"
#include <BRepFilletAPI_MakeFillet2d.hxx>
#include <ChFi2d_FilletAPI.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <ChFi2d_AnaFilletAlgo.hxx>
#include <ChFi2d_FilletAlgo.hxx>
#include <BRepExtrema_ExtCC.hxx>
#include <BRepExtrema_ExtPC.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepLib_FuseEdges.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <AIS_InteractiveContext.hxx>
// CZXDialog dialog

class CZXDialog : public CDialog
{
	DECLARE_DYNAMIC(CZXDialog)

public:
	CZXDialog(Handle_AIS_InteractiveContext myContext,CWnd* pParent = NULL);   // standard constructor
	virtual ~CZXDialog();

	// Dialog Data
	enum { IDD = IDD_DIALOG_ZX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	CColoredShapes* m_pcoloredshapeList;
	Handle(AIS_InteractiveContext) myAISContext;

	DECLARE_MESSAGE_MAP()

public:
	void OnBearingVolute(double Z1,double X1,double Z2,double X2, double width,double exhaustFlankHeight,double bearingFlankHeight,double bearingSideAngle,double exhaustSideAngle,double wholeVoluteArea,
		double tipRadius,double dividerWallHeight,double dividerAngle,double exhaustThickness,
		double bearingThickness,double transitionPartLength,double exitPipeLength,double exitDividerAngle,double voluteRadius,double exitDividerWallWidth,double exitPipeRadius,
		double toungAreaPercentage);

	TopoDS_Wire getDualVoluteCrossSection(double Z1,double X1,double Z2,double X2,double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle, double exhaustSideAngle,
		double expectedArea, double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle,
		double exhaustThickness, double bearingThickness);

	gp_Pnt getMinimumDistancePoint(TopoDS_Edge,TopoDS_Edge);
	gp_Pnt getMinimumDistancePoint(TopoDS_Edge,TopoDS_Vertex);


	void makeTwoLobeVolute();
	void displayPoint(gp_Pnt);
	void displayPoint( Quantity_NameOfColor,gp_Pnt);
	void displayEdge(TopoDS_Edge);
	 double getSurfaceArea(TopoDS_Wire);

	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedButtonVolute();
};
