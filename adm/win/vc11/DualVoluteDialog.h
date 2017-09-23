#pragma once
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

// CDualVoluteDialog dialog

class CDualVoluteDialog : public CDialog
{
	DECLARE_DYNAMIC(CDualVoluteDialog)

public:
	CDualVoluteDialog(Handle_AIS_InteractiveContext,CWnd* pParent = NULL);   // standard constructor
	virtual ~CDualVoluteDialog();

// Dialog Data
	enum { IDD = IDD_DIALOG3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	CColoredShapes* m_pcoloredshapeList;
	  Handle(AIS_InteractiveContext) myAISContext;
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();

public:
	//initiating program,create compund shape
	void OnBearingVolute();

	//creates cross section for a given area
	TopoDS_Wire getDualVoluteCrossSection(double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle, double exhaustSideAngle,
		double expectedArea, double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle,
		double exhaustThickness, double bearingThickness);

	//create volute cross section for transition part
	void createTrasitionVoluteCrossSection(TopoDS_Wire& w1,TopoDS_Wire& w2,TopoDS_Edge& p2q2e1,gp_Pnt& leftTopPoint,gp_Pnt& rightTopPoint,gp_Vec& leftBottomVec,gp_Vec& dividerLineVec,
		double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle, double exhaustSideAngle, double wholeVoluteArea,
		double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle, double exhaustThickness, double bearingThickness);

	//create exit pipe cross section for trasition part
	void getExitPipeTransitionCrossSection(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,double exitDividerWallWidth,
		double exitPipeRadius,gp_Vec divideLineVec);

	//create volute-exit pipe transition part
	void createTransitionExitPipePart(TopoDS_Shape& s1,TopoDS_Shape& s2,TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt& centre,TopoDS_Wire crossSection27,TopoDS_Wire leftVoluteWire,
		TopoDS_Wire rightVoluteWire ,TopoDS_Edge p2q2Edge,gp_Vec leftBottomVector ,gp_Vec divideLineVec , double exitPipeRadius,
		double exitDividerAngle,double maxWidthOfDividerWall,double transitionPartLength, double voluteRadius );

	// create exit pipe ending cross-section
	void getExitPipeEnding(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,
		double exitDividerWallWidth,double exitPipeRadius,gp_Vec divideLineVec);

	//create exit pipe
	void createExitPipeEndigPart(TopoDS_Shape& s1,TopoDS_Shape& s2,TopoDS_Wire leftTransition,TopoDS_Wire rightTransition,gp_Pnt centrePointOfCircle,
		gp_Vec leftBottomVector,gp_Vec divideLineVec ,double exitDividerWallWidth,double exitPipeRadius,double exitDividerAngle, 
		double transitionPartLength,double exitPipeLength);


	//calculate maximum height of a volute cross-section
	double getMaximumHeight(TopoDS_Wire);

	double getTrapezuimHeight(double,double,double,double);
	
	gp_Pnt getMinimumDistancePoint(TopoDS_Edge,TopoDS_Edge);
	gp_Pnt getMinimumDistancePoint(TopoDS_Edge,TopoDS_Vertex);
	TopoDS_Edge convertTrimmToBezier(Handle_Geom_Curve,gp_Vec,gp_Vec,double);
	double getDividerWallMaximumWidth(TopoDS_Wire);
	void MergeEdges(TopoDS_Edge edge1, TopoDS_Edge edge2);

	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedButton1();
};
