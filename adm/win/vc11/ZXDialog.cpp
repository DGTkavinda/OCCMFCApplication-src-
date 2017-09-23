// ZXDialog.cpp : implementation file
//

#include "stdafx.h"
#include "ZXDialog.h"
#include "afxdialogex.h"
#include "ImportExportDoc.h"
#define PI 3.14159265

// CZXDialog dialog

IMPLEMENT_DYNAMIC(CZXDialog, CDialog)

	CZXDialog::CZXDialog(Handle_AIS_InteractiveContext myContext,CWnd* pParent /*=NULL*/)
	: CDialog(CZXDialog::IDD, pParent)
{
	m_pcoloredshapeList = new CColoredShapes();
	myAISContext = myContext;
	myAISContext->SetDisplayMode(AIS_Shaded);
}

CZXDialog::~CZXDialog()
{
}

void CZXDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

}


BEGIN_MESSAGE_MAP(CZXDialog, CDialog)
	ON_BN_CLICKED(IDOK, &CZXDialog::OnBnClickedOk)
	ON_BN_CLICKED(IDC_BUTTON_Volute, &CZXDialog::OnBnClickedButtonVolute)
END_MESSAGE_MAP()


// CZXDialog message handlers


void CZXDialog::OnBnClickedOk()
{


	// TODO: Add your control notification handler code here
	CDialog::OnOK();
}




void CZXDialog::makeTwoLobeVolute()
{

	double width;
	double exhaustFlankHeight;
	double bearingFlankHeight;
	double bearingSideAngle;
	double exhaustSideAngle;
	double wholeVoluteArea;	
	double tipRadius;
	double dividerWallHeight;
	double dividerAngle;
	double exhaustThickness;
	double bearingThickness;
	double toungAreaPercentage;

	double exitPipeRadius;
	double transitionPartLength;
	double exitPipeLength;
	double exitDividerAngle;
	double voluteRadius;
	double exitDividerWallWidth;
	double inletZ1;
	double inletX1;
	double inletZ2;
	double inletX2;



	double areaRatio;

	width=10*2;
	bearingFlankHeight=13.526;
	exhaustFlankHeight=10;

	bearingSideAngle=41;
	exhaustSideAngle=30;
	wholeVoluteArea=2000;

	tipRadius=2;
	dividerWallHeight=15;

	dividerAngle=20;
	exhaustThickness=15;
	bearingThickness=15;
	toungAreaPercentage=0.01;
	transitionPartLength=-150;
	exitPipeLength=200;
	exitDividerAngle=90;
	exitDividerWallWidth=15;
	exitPipeRadius=width*3;
	voluteRadius=width*5;

	inletZ1=10;
	inletX1=10;
	inletZ2=10+20; 
	inletX2=10;

	//displayPoint(gp_Pnt(0,0,0));

	OnBearingVolute(inletZ1,inletX1,inletZ2,inletX2,width, exhaustFlankHeight, bearingFlankHeight, bearingSideAngle, exhaustSideAngle, wholeVoluteArea,
		tipRadius, dividerWallHeight, dividerAngle, exhaustThickness,bearingThickness, transitionPartLength, exitPipeLength,
		exitDividerAngle, voluteRadius, exitDividerWallWidth, exitPipeRadius,toungAreaPercentage);


	return;
}

void CZXDialog::OnBearingVolute(double Z1,double X1,double Z2,double X2,double width,double exhaustFlankHeight,double bearingFlankHeight,double bearingSideAngle,double exhaustSideAngle,double wholeVoluteArea,
								double tipRadius,double dividerWallHeight,double dividerAngle,double exhaustThickness,
								double bearingThickness,double transitionPartLength,double exitPipeLength,double exitDividerAngle,double voluteRadius,double exitDividerWallWidth,
								double exitPipeRadius,double toungAreaPercentage)
{

	double bearingFlankHeightGap=abs(bearingFlankHeight-exhaustFlankHeight);
	double expectedArea=1000;
	width=Z2-Z1;

	TopoDS_Wire newCrossSection;


	newCrossSection=getDualVoluteCrossSection(Z1,X1,Z2,X2,width,bearingFlankHeightGap, bearingFlankHeight, bearingSideAngle, exhaustSideAngle,expectedArea, tipRadius, dividerWallHeight, exhaustFlankHeight,dividerAngle, exhaustThickness,bearingThickness);



}


TopoDS_Wire CZXDialog::getDualVoluteCrossSection(double Z1,double X1,double Z2,double X2,double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle, double exhaustSideAngle, 
												 double expectedArea,double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle,
												 double exhaustThickness,double bearingThickness)
{



	double baseAngle;
	double baseWidth;

	baseAngle=atan2(bearingFlankHeightGap,width)*180/PI;
	baseWidth=sqrt(width*width+bearingFlankHeightGap*bearingFlankHeightGap);

	double bearingSideAngTan=tan((bearingSideAngle-baseAngle)*PI/180);
	double exhaustSideAngTan=tan((exhaustSideAngle-baseAngle)*PI/180);

	double bearingSideInitialScaler=bearingThickness;
	double exhaustSideInitialScaler=exhaustThickness;



	gp_Pnt q0(X1,0,Z1);
	gp_Pnt p0(X2,0,Z2);
	gp_Pnt q1(X1+bearingFlankHeight,0,Z1);
	gp_Pnt p1(X2+exhaustFlankHeight,0,Z2);

	displayPoint(Quantity_NOC_IVORY,gp_Pnt(0,0,0));

	displayPoint(Quantity_NOC_YELLOW,q0);
	


	TopoDS_Edge q0q1Edge=BRepBuilderAPI_MakeEdge(q0,q1);
	TopoDS_Edge p0p1Edge=BRepBuilderAPI_MakeEdge(p0,p1);

	gp_Vec leftHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,-10));
	leftHorizontalVec.Normalize();
	gp_Ax1 q1YAxis(q1,gp_Dir(0,1,0));
	gp_Vec q1q2Vec=leftHorizontalVec.Rotated(q1YAxis,(360-bearingSideAngle)*PI/180);
	q1q2Vec.Multiply(bearingSideInitialScaler);
	gp_Pnt q2=q1.Translated(q1q2Vec);

	gp_Vec rightHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,10));
	rightHorizontalVec.Normalize();
	gp_Ax1 p1YAxis(p1,gp_Dir(0,1,0));
	gp_Vec p1p2Vec=rightHorizontalVec.Rotated(p1YAxis,(exhaustSideAngle)*PI/180);
	p1p2Vec.Multiply(exhaustSideInitialScaler);
	gp_Pnt p2=p1.Translated(p1p2Vec);

	TopoDS_Edge q1q2Edge=BRepBuilderAPI_MakeEdge(q1,q2);
	TopoDS_Edge p1p2Edge=BRepBuilderAPI_MakeEdge(p1,p2);
	TopoDS_Edge p1q1Edge=BRepBuilderAPI_MakeEdge(p1,q1);


	gp_Vec q1q2Perpendiculer=q1q2Vec.Rotated(q1YAxis,(270)*PI/180);
	q1q2Perpendiculer.Normalize();
	q1q2Perpendiculer.Multiply(bearingThickness+tipRadius);
	gp_Pnt q1q2ParallelPoint=q1.Translated(q1q2Perpendiculer);
	gp_Dir q1q2Dir(q1q2Vec);
	gp_Lin q1q2ParallelLine(q1q2ParallelPoint,q1q2Dir); 
	TopoDS_Edge q1q2ParallelEdge=BRepBuilderAPI_MakeEdge(q1q2ParallelLine);




	gp_Vec p1p2Perpendiculer=p1p2Vec.Rotated(p1YAxis,90*PI/180);
	p1p2Perpendiculer.Normalize();
	p1p2Perpendiculer.Multiply(exhaustThickness+tipRadius);
	gp_Pnt p1p2ParallePoint=p1.Translated(p1p2Perpendiculer);
	gp_Dir p1p2Dir(p1p2Vec);
	gp_Lin p1p2ParallelLine(p1p2ParallePoint,p1p2Dir);
	TopoDS_Edge p1p2ParallelEdge=BRepBuilderAPI_MakeEdge(p1p2ParallelLine);

	gp_Pnt parallelEdgesIntersectionPoint=getMinimumDistancePoint(q1q2ParallelEdge,p1p2ParallelEdge);
	TopoDS_Vertex intersectionVertex=BRepBuilderAPI_MakeVertex(parallelEdgesIntersectionPoint);

	gp_Pnt q1q2minimumDistPnt=getMinimumDistancePoint(q1q2Edge,intersectionVertex);
	gp_Pnt p1p2minimumDistPnt=getMinimumDistancePoint(p1p2Edge,intersectionVertex);
	gp_Pnt p1q1minimumDistPnt=getMinimumDistancePoint(p1q1Edge,intersectionVertex);
	gp_Pnt r1=p1q1minimumDistPnt;
	
	gp_Vec p1q1Vec(p1,q1);
	
	p1q1Vec.Normalize();
	p1q1Vec.Multiply(tipRadius);
	gp_Pnt leftBottomPointOfDividerWall=parallelEdgesIntersectionPoint.Translated(p1q1Vec);
	p1q1Vec.Reverse();
	gp_Pnt rightBottomPointOfDividerWall=parallelEdgesIntersectionPoint.Translated(p1q1Vec);

	double leftSideLengthOfDividerWall=q1q2minimumDistPnt.Distance(q2);
	double rightSideLengthOfDividerWall=p1p2minimumDistPnt.Distance(p2);

	gp_Vec dividerLineVec(p1q1minimumDistPnt,parallelEdgesIntersectionPoint);
	dividerLineVec.Normalize();
	
	gp_Ax1 leftBottomPointAxis(leftBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec leftDividerWallScaleVec=dividerLineVec.Rotated(leftBottomPointAxis,((dividerAngle/2))*PI/180);
	leftDividerWallScaleVec.Normalize();
	leftDividerWallScaleVec.Multiply(leftSideLengthOfDividerWall);
	gp_Pnt leftTopPointOfDividerWall=leftBottomPointOfDividerWall.Translated(leftDividerWallScaleVec);
	TopoDS_Edge leftDividerWallEdge=BRepBuilderAPI_MakeEdge(leftTopPointOfDividerWall,leftBottomPointOfDividerWall);
	
	gp_Ax1 rightBottomPointAxis(rightBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec rightDividerWallScaleVec=dividerLineVec.Rotated(rightBottomPointAxis,(360-(dividerAngle/2))*PI/180);
	rightDividerWallScaleVec.Normalize();
	rightDividerWallScaleVec.Multiply(rightSideLengthOfDividerWall);
	gp_Pnt rightTopPointOfDividerWall=rightBottomPointOfDividerWall.Translated(rightDividerWallScaleVec);
	TopoDS_Edge rightDividerWallEdge=BRepBuilderAPI_MakeEdge(rightTopPointOfDividerWall,rightBottomPointOfDividerWall);


	Standard_Real U1;
	Standard_Real U2;

	gp_Pnt q3=leftTopPointOfDividerWall;
	gp_Pnt q4=leftBottomPointOfDividerWall;
	
	gp_Pnt p3=rightTopPointOfDividerWall;
	

	gp_Vec q2q3Vec(q2,q3);
	TopoDS_Edge q2q3Edge=BRepBuilderAPI_MakeEdge(q2,q3);
	Handle_Geom_Curve q2q3LineCurve=BRep_Tool::Curve(q2q3Edge,U1,U2);
	gp_Pnt q2q3CentrePoint;
	q2q3LineCurve->D0(U2/2,q2q3CentrePoint);

	gp_Lin q2q3ParallelLine(leftBottomPointOfDividerWall,q2q3Vec);
	TopoDS_Edge q2q3ParallelEdge=BRepBuilderAPI_MakeEdge(q2q3ParallelLine);
	gp_Pnt q2q3ParallelPointOnq1q2=getMinimumDistancePoint(q1q2Edge,q2q3ParallelEdge);
	gp_Pnt q5=q2q3ParallelPointOnq1q2;
	gp_Pnt q6=q1q2minimumDistPnt;

	TopoDS_Edge q1q5Edge=BRepBuilderAPI_MakeEdge(q1,q5);
	TopoDS_Edge q1r1Edge=BRepBuilderAPI_MakeEdge(q1,r1);
	TopoDS_Edge r1q4Edge=BRepBuilderAPI_MakeEdge(r1,q4);
	TopoDS_Edge q4q5Edge=BRepBuilderAPI_MakeEdge(q4,q5);

	BRepBuilderAPI_MakeWire makeLeftConstantWire;
	makeLeftConstantWire.Add(q1q5Edge);
	makeLeftConstantWire.Add(q1r1Edge);
	makeLeftConstantWire.Add(r1q4Edge);
	makeLeftConstantWire.Add(q4q5Edge);

	TopoDS_Wire leftConstantWire=makeLeftConstantWire;
	double leftConstantArea=getSurfaceArea(leftConstantWire);
	double q4q3AngleToq5q4= q2q3Vec.Angle(leftDividerWallScaleVec);
	gp_Vec q2q3VecRevesed=q2q3Vec.Reversed();
	double q5q2AngleToq5q4= q2q3VecRevesed.Angle(q1q2Vec);
	gp_Vec q4q6Vec(q4,q6);
	double q4q5AngleToQ4q6=q2q3VecRevesed.Angle(q4q6Vec);



	CString str;
	str.Format(_T("angle %g"),q4q5AngleToQ4q6);
	AfxMessageBox(str);






	displayPoint(rightTopPointOfDividerWall);
	displayPoint(leftTopPointOfDividerWall);

	displayPoint(Quantity_NOC_GREEN,leftBottomPointOfDividerWall);
	displayPoint(Quantity_NOC_GREEN,rightBottomPointOfDividerWall);

	

	displayPoint(Quantity_NOC_YELLOW,parallelEdgesIntersectionPoint);

	TopoDS_Edge p0q0Edge=BRepBuilderAPI_MakeEdge(p0,q0);
	gp_Pnt tran=q1.Translated(q1q2Perpendiculer);




	TopoDS_Edge xAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(10,0,0));
	TopoDS_Edge yAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(0,10,0));
	TopoDS_Edge zAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(0,0,10));

	//m_pcoloredshapeList->Add(Quantity_NOC_RED,p1p2ParallelEdge);
	//m_pcoloredshapeList->Add(Quantity_NOC_RED,q1q2ParallelEdge);
	


	m_pcoloredshapeList->Add(Quantity_NOC_RED,rightDividerWallEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,leftDividerWallEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,p1p2Edge);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,q1q2Edge);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,p0p1Edge);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,q0q1Edge);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,xAxisEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,yAxisEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_GREEN,zAxisEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_IVORY,leftConstantWire);

	





	m_pcoloredshapeList->Display(myAISContext);
	TopoDS_Wire wire;
	CImportExportDoc::Fit();
	return wire;

}



void CZXDialog::displayPoint(gp_Pnt pnt)
{
	TopoDS_Vertex vert =BRepBuilderAPI_MakeVertex(pnt);


	m_pcoloredshapeList->Add(Quantity_NOC_MATRABLUE,vert);




}

void CZXDialog::displayPoint( Quantity_NameOfColor color,gp_Pnt pnt)
{
	TopoDS_Vertex vert =BRepBuilderAPI_MakeVertex(pnt);

	m_pcoloredshapeList->Add(color,vert);


}

gp_Pnt CZXDialog::getMinimumDistancePoint(TopoDS_Edge edge1,TopoDS_Edge edge2)
	{

		TopoDS_Wire wire1=BRepBuilderAPI_MakeWire(edge1);
		TopoDS_Wire wire2=BRepBuilderAPI_MakeWire(edge2);


		BRepExtrema_DistShapeShape minimumDist(wire1,wire2,Extrema_ExtFlag_MIN,Extrema_ExtAlgo_Grad);
		gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

		return divideIntersectionPointOnCurve;
	}

gp_Pnt CZXDialog::getMinimumDistancePoint(TopoDS_Edge edge,TopoDS_Vertex vertex)
	{

		BRepExtrema_DistShapeShape minimumDist(edge,vertex,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
		gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

		return divideIntersectionPointOnCurve;
	}

double CZXDialog::getSurfaceArea(TopoDS_Wire wireShape){

		BRepBuilderAPI_MakeFace wireFace(wireShape);
		TopoDS_Face face =wireFace;
		BRepGProp gprop; 
		GProp_GProps surfaceProps;
		gprop.SurfaceProperties(wireFace,surfaceProps);
		double surfaceArea=surfaceProps.Mass();
		return surfaceArea;

	}


void CZXDialog::OnBnClickedButtonVolute()
{


	makeTwoLobeVolute();
	// TODO: Add your control notification handler code here
}
