// DualVoluteDialog.cpp : implementation file
//

#include "stdafx.h"
#include "DualVoluteDialog.h"
#include "afxdialogex.h"


#include <Standard_Address.hxx>
#include <Standard_Stream.hxx>
#include <Standard_Transient.hxx>
#include <type_traits>
#include <Geom_TrimmedCurve.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeSegment.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgo_Fuse.hxx>
#include <Standard_Macro.hxx>
#include <Standard.hxx>
#include <Standard_Handle.hxx>
#include <Standard_Transient.hxx>
#include <Standard_OStream.hxx>
#include <typeinfo>
#include <Standard_Type.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <Geom2d_Ellipse.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_BuilderAlgo.hxx"
#include "BRepFeat_MakeCylindricalHole.hxx"
#include "windef.h"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#define PI 3.14159265

// CDualVoluteDialog dialog

IMPLEMENT_DYNAMIC(CDualVoluteDialog, CDialog)

	CDualVoluteDialog::CDualVoluteDialog(Handle_AIS_InteractiveContext myContext,CWnd* pParent /*=NULL*/)
	: CDialog(CDualVoluteDialog::IDD, pParent)
{
	m_pcoloredshapeList = new CColoredShapes();
	myAISContext = myContext;
	myAISContext->SetDisplayMode(AIS_Shaded);
}

CDualVoluteDialog::~CDualVoluteDialog()
{
}

void CDualVoluteDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CDualVoluteDialog, CDialog)
	ON_BN_CLICKED(IDOK, &CDualVoluteDialog::OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, &CDualVoluteDialog::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_BUTTON1, &CDualVoluteDialog::OnBnClickedButton1)
END_MESSAGE_MAP()


// CDualVoluteDialog message handlers


void CDualVoluteDialog::OnBnClickedOk()
{

	CString str="sdsdf";
	AfxMessageBox(str);
	// TODO: Add your control notification handler code here
	CDialog::OnOK();

}


void CDualVoluteDialog::OnBnClickedButton1()
{
	OnBearingVolute();
	CString str="sdsdf";
	AfxMessageBox(str);

	// TODO: Add your control notification handler code here
}


void CDualVoluteDialog::OnBearingVolute(){

	/*         q2                     p1
	\         |         /
	\		   |	    /
	\       |       /
	\      |      /
	\_____|_____/
	q1              p1

	q1q2- bearing side
	p1p2- exhaust side


	*/

	double width;
	double exhaustFlankHeight;
	double bearingFlankHeight;
	double bearingSideAngle;
	double exhaustSideAngle;
	double wholeVoluteArea;
	double wholeVoluteTrapziumHeight;
	double trapeziumWidth;
	double trapeziumAngle;
	double tipRadius;
	double dividerWallHeight;

	double dividerAngle;
	double exhaustThickness;
	double bearingThickness;
	double transitionPartLength;
	double exitPipeLength;
	double exitDividerAngle;
	double voluteRadius;
	double exitDividerWallWidth;
	double exitPipeRadius;
	double toungAreaPercentage;
	double bearingFlankHeightGap;

	width=10*2;
	bearingFlankHeight=13.526;
	exhaustFlankHeight=10;
	bearingFlankHeightGap=abs(bearingFlankHeight-exhaustFlankHeight);
	bearingSideAngle=41;
	exhaustSideAngle=30;
	wholeVoluteArea=2000;

	tipRadius=0.5;
	dividerWallHeight=15;

	dividerAngle=10;
	exhaustThickness=15;
	bearingThickness=15;
	toungAreaPercentage=0.01;

	transitionPartLength=-150;
	exitPipeLength=-200;
	exitDividerAngle=90;
	exitDividerWallWidth=15;
	exitPipeRadius=width*3;
	voluteRadius=width*5;


	// volute
	BRepOffsetAPI_ThruSections sections;
	gp_Trsf transfer;
	gp_Ax1 rotationAxis(gp_Pnt(-voluteRadius,0,0),gp_Dir(0,0,1));

	int numberOfCrossSections=36;
	double toungArea=wholeVoluteArea*0.01;
	TopoDS_Wire arrayOfCrossSections[36];
	TopoDS_Wire newCrossSection;
	double areaReductionRate=(wholeVoluteArea-toungArea)/numberOfCrossSections;
	double rotationAngle=(360/numberOfCrossSections)*PI/180;


	double expectedArea=wholeVoluteArea-(areaReductionRate*0);
	newCrossSection=getDualVoluteCrossSection(width,bearingFlankHeightGap, bearingFlankHeight, bearingSideAngle, exhaustSideAngle,expectedArea, tipRadius, dividerWallHeight, exhaustFlankHeight,dividerAngle, exhaustThickness,bearingThickness);
	transfer.SetRotation(rotationAxis,360*PI/180);
	BRepBuilderAPI_Transform rotated(newCrossSection,transfer);
	arrayOfCrossSections[0]=TopoDS::Wire(rotated.Shape());// 360 crossSection
	sections.AddWire(arrayOfCrossSections[0]);
	//	m_pcoloredshapeList->Add(Quantity_NOC_IVORY,newCrossSection);
	for(int i=1;i<=numberOfCrossSections-1;i++){

		double expectedArea=wholeVoluteArea-(areaReductionRate*i);
		newCrossSection=getDualVoluteCrossSection(width,bearingFlankHeightGap, bearingFlankHeight, bearingSideAngle, exhaustSideAngle,expectedArea, tipRadius, dividerWallHeight, exhaustFlankHeight,dividerAngle, exhaustThickness,bearingThickness);
		transfer.SetRotation(rotationAxis,rotationAngle*i);
		BRepBuilderAPI_Transform rotated(newCrossSection,transfer);
		arrayOfCrossSections[i]=TopoDS::Wire(rotated.Shape());
		sections.AddWire(arrayOfCrossSections[i]);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,arrayOfCrossSections[i]);

	}

	double expectedArea1=wholeVoluteArea-(areaReductionRate*numberOfCrossSections);
	newCrossSection=getDualVoluteCrossSection(width,bearingFlankHeightGap, bearingFlankHeight, bearingSideAngle, exhaustSideAngle,expectedArea1, tipRadius, dividerWallHeight, exhaustFlankHeight,dividerAngle, exhaustThickness,bearingThickness); 
	transfer.SetRotation(rotationAxis,360*PI/180);
	BRepBuilderAPI_Transform rotated1(newCrossSection,transfer);
	arrayOfCrossSections[numberOfCrossSections-1]=TopoDS::Wire(rotated1.Shape());

	sections.AddWire(arrayOfCrossSections[numberOfCrossSections-1]);//toung crossSection
	sections.Build();
	TopoDS_Shape VoluteShape(sections.Shape());


	TopoDS_Edge xAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(10,0,0));
	TopoDS_Edge yAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(0,10,0));
	TopoDS_Edge zAxisEdge=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(0,0,10));

	m_pcoloredshapeList->Add(Quantity_NOC_RED,xAxisEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,yAxisEdge);
	m_pcoloredshapeList->Add(Quantity_NOC_GREEN,zAxisEdge);


		m_pcoloredshapeList->Display(myAISContext);
	
	//exit pipe section

	TopoDS_Wire leftVoluteWire;
	TopoDS_Wire rightVoluteWire;
	TopoDS_Edge p2q2Edge;
	gp_Pnt leftTopPointofDividerWall;
	gp_Pnt rightTopPointOfDividerWall;
	gp_Vec leftBottomVector;
	gp_Vec divideLineVec;
	TopoDS_Wire leftTransition;
	TopoDS_Wire rightTransition;
	gp_Pnt centrePointOfCircle;

	//Transition
	createTrasitionVoluteCrossSection(leftVoluteWire,rightVoluteWire,p2q2Edge,leftTopPointofDividerWall,rightTopPointOfDividerWall,leftBottomVector,divideLineVec,width,bearingFlankHeightGap, bearingFlankHeight, bearingSideAngle, exhaustSideAngle,wholeVoluteArea, tipRadius, dividerWallHeight, exhaustFlankHeight,dividerAngle, exhaustThickness,bearingThickness);

	double maxWidthOfDividerWall=leftTopPointofDividerWall.Distance(rightTopPointOfDividerWall);	

	TopoDS_Shape leftTransitionPartShape;
	TopoDS_Shape rightTransitionPartShape;

	createTransitionExitPipePart(leftTransitionPartShape,rightTransitionPartShape,leftTransition,rightTransition,centrePointOfCircle,arrayOfCrossSections[27],leftVoluteWire,rightVoluteWire,p2q2Edge,leftBottomVector,divideLineVec,exitPipeRadius,exitDividerAngle,maxWidthOfDividerWall,transitionPartLength,voluteRadius);



	//ExitPipe

	TopoDS_Shape leftExitPipeEnding;
	TopoDS_Shape rightExitPipeEnding;

	createExitPipeEndigPart(leftExitPipeEnding,rightExitPipeEnding,leftTransition,rightTransition,centrePointOfCircle,leftBottomVector,divideLineVec,exitDividerWallWidth,exitPipeRadius,exitDividerAngle,transitionPartLength,exitPipeLength);

	TopoDS_Compound transitionCom;
	TopoDS_Builder trBuilder;
	trBuilder.MakeCompound(transitionCom);
	trBuilder.Add(transitionCom,leftTransitionPartShape);
	trBuilder.Add(transitionCom,rightTransitionPartShape);
	trBuilder.Add(transitionCom,leftExitPipeEnding);
	trBuilder.Add(transitionCom,rightExitPipeEnding);
	trBuilder.Add(transitionCom,VoluteShape);


	BRepTools::Write(VoluteShape,"D:/Breps/newVolute/Volute.brep");
	BRepTools::Write(leftTransitionPartShape,"D:/Breps/newVolute/leftTransitionPartShape.brep");
	BRepTools::Write(rightTransitionPartShape,"D:/Breps/newVolute/rightTransitionPartShape.brep");
	BRepTools::Write(leftExitPipeEnding,"D:/Breps/newVolute/leftExitPipeEnding.brep");
	BRepTools::Write(rightExitPipeEnding,"D:/Breps/newVolute/rightExitPipeEnding.brep");
	BRepTools::Write(transitionCom,"D:/Breps/newVolute/completeShape.brep");
	
	

}

double CDualVoluteDialog::getTrapezuimHeight(double area,double width, double ang1Tan, double ang2Tan){



	double c=area*2;
	double a=1/ang1Tan+1/ang2Tan;
	double b= 2*width;

	double determinaterPart =(b*b) + (4*a*c); 
	double diterminator= sqrt(determinaterPart);
	double root1=((-b)+(diterminator))/(2*a);

	return root1;


}


TopoDS_Wire CDualVoluteDialog::getDualVoluteCrossSection(double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle, double exhaustSideAngle, 
														 double expectedArea,double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle,
														 double exhaustThickness,double bearingThickness)
{


	//trapezium common to both lobes to create outer wires
	double trapeziumWidth;
	double trapeziumAngle;
	double wholeVoluteTrapziumHeight;

	trapeziumWidth=sqrt(width*width+bearingFlankHeightGap*bearingFlankHeightGap);
	trapeziumAngle=atan2(bearingFlankHeightGap,width)*180/PI;
	double bearingSideAngTan=tan((bearingSideAngle-trapeziumAngle)*PI/180);
	double exhaustSideAngTan=tan((exhaustSideAngle-trapeziumAngle)*PI/180);
	wholeVoluteTrapziumHeight=getTrapezuimHeight(expectedArea,trapeziumWidth,bearingSideAngTan,exhaustSideAngle);

	double h=wholeVoluteTrapziumHeight;
	double A1=exhaustSideAngle;
	double A2=bearingSideAngle;
	double alfa=trapeziumAngle;
	double initialLengthOfOuterLine=width/12;

	double cosA1=cos(A1*PI/180);
	double sinA1=sin(A1*PI/180);

	double cosA2=cos(A2*PI/180);
	double sinA2=cos(A2*PI/180);

	double tanAlfa=tan(alfa*PI/180);
	double cosAlfa=cos(alfa*PI/180);

	double p2x=(h*cosA1)/(sin((A1+alfa)*PI/180));
	double p2y=(h*sinA1)/(sin((A1+alfa)*PI/180));

	double q2x=(h*cosA2)/(sin((A2-alfa)*PI/180));
	double q2y=(h*sinA2)/(sin((A2-alfa)*PI/180));


	gp_Pnt p0(-exhaustFlankHeight,0,0);
	gp_Pnt q0(-exhaustFlankHeight,0,-width);
	gp_Vec leftHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,-1));
	gp_Vec rightHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,1));

	gp_Pnt p1(0,0,0);
	gp_Pnt p2(p2y,0,p2x);

	gp_Pnt q1(bearingFlankHeightGap,0,-width);
	gp_Pnt q2(bearingFlankHeightGap+q2y,0,-(width+q2x));
	gp_Pnt r1((exhaustThickness*cosA1),0,-(exhaustThickness*sinA1));

	gp_Ax1 q1Axis(q1,gp_Dir(0,1,0));
	gp_Vec q1q2Vec=leftHorizontalVec.Rotated(q1Axis,(360-A2)*PI/180);
	q1q2Vec.Normalize();
	q1q2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A2-alfa)*PI/180)));
	gp_Pnt q1End=q1.Translated(q1q2Vec);
	q2=q1End;

	gp_Ax1 p1Axis(p1,gp_Dir(0,1,0));
	gp_Vec p1p2Vec=rightHorizontalVec.Rotated(q1Axis,A1*PI/180);
	p1p2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A1)*PI/180)));
	gp_Pnt p1End=p1.Translated(p1p2Vec);
	TopoDS_Edge q1q2Strait=BRepBuilderAPI_MakeEdge(p1End,p1);

	p2=p1End;

	TopoDS_Edge p1p2Edge=BRepBuilderAPI_MakeEdge(p1,p2);
	TopoDS_Edge p1q1Edge=BRepBuilderAPI_MakeEdge(p1,q1);
	TopoDS_Edge q1q2Edge=BRepBuilderAPI_MakeEdge(q1,q2);
	TopoDS_Edge p2q2Edge=BRepBuilderAPI_MakeEdge(p2,q2);
	TopoDS_Edge p1r1Edge=BRepBuilderAPI_MakeEdge(p1,r1);
	TopoDS_Edge q1r1Edge=BRepBuilderAPI_MakeEdge(q1,r1);
	TopoDS_Edge p0p1Edge=BRepBuilderAPI_MakeEdge(p0,p1);
	TopoDS_Edge q0q1Edge=BRepBuilderAPI_MakeEdge(q0,q1);

	TopoDS_Edge hirizontalLine=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,-10),gp_Pnt(0,0,10));

	TopoDS_Vertex vp1=BRepBuilderAPI_MakeVertex(p1);
	TopoDS_Vertex vp2=BRepBuilderAPI_MakeVertex(p2);
	TopoDS_Vertex vq1=BRepBuilderAPI_MakeVertex(q1);
	TopoDS_Vertex vq2=BRepBuilderAPI_MakeVertex(q2);

	BRepBuilderAPI_MakeWire trape;
	trape.Add(p1p2Edge);
	trape.Add(p1q1Edge);
	trape.Add(q1q2Edge);
	trape.Add(p2q2Edge);

	TopoDS_Wire wire=trape;

	Standard_Real U1;
	Standard_Real U2;


	Handle_Geom_Curve p1q1LineCurve=BRep_Tool::Curve(p1q1Edge,U1,U2);

	//bottom points of splitter
	Handle_Geom_Curve q1q2LineCurve=BRep_Tool::Curve(q1q2Edge,U1,U2);
	gp_Pnt q1q2pointZero;
	gp_Vec q1q2Vec1;
	gp_Vec q1q2Vec2;
	gp_Vec oZ(gp_Dir(0,1,0));
	q1q2LineCurve->D0(U1,q1q2pointZero);
	TopoDS_Vertex q1q2vertZero=BRepBuilderAPI_MakeVertex(q1q2pointZero);
	q1q2LineCurve->D2(U1,q1q2pointZero,q1q2Vec1,q1q2Vec2);
	gp_Vec q1q2crossed=q1q2Vec1.Crossed(oZ);
	q1q2crossed.Normalize();
	q1q2crossed.Multiply(bearingThickness);
	q1q2pointZero.Translate(q1q2crossed);
	q1q2vertZero=BRepBuilderAPI_MakeVertex(q1q2pointZero);
	TopoDS_Edge bearingSideBottomEdge=BRepBuilderAPI_MakeEdge(q1,q1q2pointZero);
	gp_Pnt q1q2EndPoint;
	q1q2LineCurve->D0(U2,q1q2EndPoint);
	double q1q2Length=U2;

	Handle_Geom_Curve p1p2LineCurve=BRep_Tool::Curve(p1p2Edge,U1,U2);
	gp_Pnt p1p2pointZero;
	gp_Vec p1p2Vec1;
	gp_Vec p1p2Vec2;
	p1p2LineCurve->D0(U1,p1p2pointZero);
	TopoDS_Vertex p1p2vertZero=BRepBuilderAPI_MakeVertex(p1p2pointZero);
	p1p2LineCurve->D2(U1,p1p2pointZero,p1p2Vec1,p1p2Vec2);
	gp_Vec p1p2Vec1Reversed=p1p2Vec1.Reversed();
	gp_Vec p1p2crossed=p1p2Vec1Reversed.Crossed(oZ);
	p1p2crossed.Normalize();
	p1p2crossed.Multiply(exhaustThickness);
	p1p2pointZero.Translate(p1p2crossed);
	p1p2vertZero=BRepBuilderAPI_MakeVertex(p1p2pointZero);
	TopoDS_Edge exhaustSideBottomEdge=BRepBuilderAPI_MakeEdge(p1,p1p2pointZero);
	gp_Pnt p1p2EndPoint;
	p1p2LineCurve->D0(U2,p1p2EndPoint);
	double p1p2Length=U2;

	gp_Pnt p1q1pointZero;
	gp_Vec p1q1Vec1;
	gp_Vec p1q1Vec2;
	p1q1LineCurve->D2(U1,p1q1pointZero,p1q1Vec1,p1q1Vec2);

	//offset p1p2 and q1q2 
	//To get bottom point of the divider wall
	gp_Dir p1p2Direction(p1p2Vec1);
	gp_Lin p1p2ParallelLine(p1p2pointZero,p1p2Direction);
	TopoDS_Edge p1p2ParallelEdge=BRepBuilderAPI_MakeEdge(p1p2ParallelLine);

	gp_Dir q1q2Direction(q1q2Vec1);
	gp_Lin q1q2ParallelLine(q1q2pointZero,q1q2Direction);
	TopoDS_Edge q1q2ParallelEdge=BRepBuilderAPI_MakeEdge(q1q2ParallelLine);

	gp_Pnt bottomParallelEdgesIntersectionPnt=getMinimumDistancePoint(p1p2ParallelEdge,q1q2ParallelEdge);

	gp_Vec leftBottomVector=p1q1Vec1.Normalized();
	leftBottomVector.Multiply(tipRadius);
	gp_Pnt leftBottomPiontOfDividerWall=bottomParallelEdgesIntersectionPnt.Translated(leftBottomVector);

	gp_Vec rightBottomVector=leftBottomVector.Reversed();
	gp_Pnt rightBottomPointOfDividerWall=bottomParallelEdgesIntersectionPnt.Translated(rightBottomVector); 

	TopoDS_Vertex intersectionVert=BRepBuilderAPI_MakeVertex(bottomParallelEdgesIntersectionPnt);
	gp_Pnt dividerPointOnBase=getMinimumDistancePoint(p1q1Edge,intersectionVert);

	TopoDS_Vertex LeftbottomVert=BRepBuilderAPI_MakeVertex(leftBottomPiontOfDividerWall);
	TopoDS_Vertex RightbottpmVert=BRepBuilderAPI_MakeVertex(rightBottomPointOfDividerWall);
	TopoDS_Vertex dividerOnBaseVert=BRepBuilderAPI_MakeVertex(dividerPointOnBase);

	gp_Pnt minimumDistancePointOnp1p2=getMinimumDistancePoint(p1p2Edge,RightbottpmVert);
	gp_Pnt minimumDistancePointOnq1q2=getMinimumDistancePoint(q1q2Edge,LeftbottomVert);

	TopoDS_Vertex mdVert=BRepBuilderAPI_MakeVertex(minimumDistancePointOnq1q2);

	TopoDS_Edge mdline=BRepBuilderAPI_MakeEdge(minimumDistancePointOnq1q2,leftBottomPiontOfDividerWall);
	TopoDS_Edge mdline1=BRepBuilderAPI_MakeEdge(minimumDistancePointOnp1p2,rightBottomPointOfDividerWall);

	double bearingSideOuterWallScalar=minimumDistancePointOnq1q2.Distance(q2);//divider wall lengths
	double exhaustSideOuterWallScalar=minimumDistancePointOnp1p2.Distance(p2);

	gp_Ax1 rightAx(rightBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec divideLineVec(dividerPointOnBase,bottomParallelEdgesIntersectionPnt);
	divideLineVec.Normalize();
	gp_Vec rightRotatedVec=divideLineVec.Rotated(rightAx,((360-(dividerAngle/2)))*PI/180);
	rightRotatedVec.Multiply(exhaustSideOuterWallScalar);
	gp_Pnt rightTopPointOfDividerWall=rightBottomPointOfDividerWall.Translated(rightRotatedVec);
	TopoDS_Edge rightDividerWall=BRepBuilderAPI_MakeEdge(rightTopPointOfDividerWall,rightBottomPointOfDividerWall);

	gp_Ax1 leftAx(rightBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec leftRotateVec=divideLineVec.Rotated(leftAx,(dividerAngle/2)*PI/180); 
	leftRotateVec.Multiply(bearingSideOuterWallScalar);
	gp_Pnt leftTopPointofDividerWall=leftBottomPiontOfDividerWall.Translated(leftRotateVec);
	TopoDS_Edge leftDividerWall=BRepBuilderAPI_MakeEdge(leftTopPointofDividerWall,leftBottomPiontOfDividerWall);

	TopoDS_Edge topLeftBaseEdge=BRepBuilderAPI_MakeEdge(q1q2EndPoint,leftTopPointofDividerWall);
	Handle_Geom_Curve topLeftBaseLineCurve=BRep_Tool::Curve(topLeftBaseEdge,U1,U2);
	gp_Pnt leftCircleCentrePoint;
	topLeftBaseLineCurve->D0(U2/2,leftCircleCentrePoint);
	double topLeftCircleRadius=U2/2;
	gp_Ax2 leftCircleAxis(leftCircleCentrePoint,gp_Dir(0,1,0));
	gp_Circ TopLeftCircle(leftCircleAxis,topLeftCircleRadius);
	GC_MakeArcOfCircle topLeftHalfCircle(TopLeftCircle,leftTopPointofDividerWall,q1q2EndPoint,false);
	Handle(Geom_TrimmedCurve) leftCircleArc=topLeftHalfCircle;
	TopoDS_Edge leftCircle=BRepBuilderAPI_MakeEdge(leftCircleArc);
	TopoDS_Edge leftCircleArcBazierCurveEdge=convertTrimmToBezier(leftCircleArc,q1q2Vec1,leftRotateVec,bearingSideOuterWallScalar);
	leftCircle=leftCircleArcBazierCurveEdge;

	TopoDS_Edge topRightBaseEdge=BRepBuilderAPI_MakeEdge(p1p2EndPoint,rightTopPointOfDividerWall);
	Handle_Geom_Curve topRightBaseLineCurve=BRep_Tool::Curve(topRightBaseEdge,U1,U2);
	gp_Pnt rightCircleCentrePoint;
	topRightBaseLineCurve->D0(U2/2,rightCircleCentrePoint);
	double topRightCircleRaduis=U2/2;
	gp_Ax2 rightCircleAxis(rightCircleCentrePoint,gp_Dir(0,1,0));
	gp_Circ topRightCircle(rightCircleAxis,topRightCircleRaduis);
	GC_MakeArcOfCircle topRightHalfOfCircle(topRightCircle,p1p2EndPoint,rightTopPointOfDividerWall,false);
	Handle(Geom_TrimmedCurve) rightCircleArc=topRightHalfOfCircle;
	TopoDS_Edge rightCircle=BRepBuilderAPI_MakeEdge(rightCircleArc);
	TopoDS_Edge rightCircleArcBezierCurveEdge =convertTrimmToBezier(rightCircleArc,rightRotatedVec,p1p2Vec1,exhaustSideOuterWallScalar);
	rightCircle = rightCircleArcBezierCurveEdge;
	

	gp_Ax2 bottomCircleAxis(bottomParallelEdgesIntersectionPnt,gp_Dir(0,1,0));
	gp_Circ tipcircle(bottomCircleAxis,tipRadius);
	GC_MakeArcOfCircle bottomHalfCircle(tipcircle,leftBottomPiontOfDividerWall,rightBottomPointOfDividerWall,false);
	Handle(Geom_TrimmedCurve) tipCircleArc=bottomHalfCircle;
	TopoDS_Edge bottomTipCircle=BRepBuilderAPI_MakeEdge(tipCircleArc);

	TopoDS_Compound compound;
	TopoDS_Builder abuilder;
	abuilder.MakeCompound(compound);
	abuilder.Add(compound,p1p2Edge);
	abuilder.Add(compound,q1q2Edge);
	abuilder.Add(compound,p1q1Edge);
	abuilder.Add(compound,leftCircle);
	abuilder.Add(compound,rightCircle);
	abuilder.Add(compound,bottomTipCircle);
	abuilder.Add(compound,leftDividerWall);
	abuilder.Add(compound,rightDividerWall);
	abuilder.Add(compound,mdline);
	abuilder.Add(compound,mdline1);

	BRepBuilderAPI_MakeWire completeCrossSection;

	completeCrossSection.Add(q0q1Edge);
	completeCrossSection.Add(q1q2Edge);
	completeCrossSection.Add(leftCircle);
	completeCrossSection.Add(leftDividerWall);
	completeCrossSection.Add(bottomTipCircle);
	completeCrossSection.Add(rightDividerWall);
	completeCrossSection.Add(rightCircle);
	completeCrossSection.Add(p1p2Edge);
	completeCrossSection.Add(p0p1Edge);


	TopoDS_Wire completeCrossSectionWire=completeCrossSection;
	// 	m_pcoloredshapeList->Add(Quantity_NOC_GREEN, completeCrossSection);
	return completeCrossSectionWire;
}


void CDualVoluteDialog::createTrasitionVoluteCrossSection(TopoDS_Wire& w1,TopoDS_Wire& w2,TopoDS_Edge& p2q2e1,gp_Pnt& leftTopPoint,gp_Pnt& rightTopPoint,
														  gp_Vec& leftBottomVec,gp_Vec& dividerLineVec,double width,double bearingFlankHeightGap, double bearingFlankHeight, double bearingSideAngle,
														  double exhaustSideAngle, double wholeVoluteArea, double tipRadius,double dividerWallHeight, double exhaustFlankHeight, double dividerAngle, 
														  double exhaustThickness, double bearingThickness)
{

	double trapeziumWidth;
	double trapeziumAngle;
	double wholeVoluteTrapziumHeight;

	trapeziumWidth=sqrt(width*width+bearingFlankHeightGap*bearingFlankHeightGap);
	trapeziumAngle=atan2(bearingFlankHeightGap,width)*180/PI;
	double bearingSideAngTan=tan((bearingSideAngle-trapeziumAngle)*PI/180);
	double exhaustSideAngTan=tan((exhaustSideAngle-trapeziumAngle)*PI/180);
	wholeVoluteTrapziumHeight=getTrapezuimHeight(wholeVoluteArea,trapeziumWidth,bearingSideAngTan,exhaustSideAngle);

	double h=wholeVoluteTrapziumHeight;
	double A1=exhaustSideAngle;
	double A2=bearingSideAngle;
	double alfa=trapeziumAngle;
	double initialLengthOfOuterLine=width/12;

	double cosA1=cos(A1*PI/180);
	double sinA1=sin(A1*PI/180);

	double cosA2=cos(A2*PI/180);
	double sinA2=cos(A2*PI/180);

	double tanAlfa=tan(alfa*PI/180);
	double cosAlfa=cos(alfa*PI/180);

	double p2x=(h*cosA1)/(sin((A1+alfa)*PI/180));
	double p2y=(h*sinA1)/(sin((A1+alfa)*PI/180));

	double q2x=(h*cosA2)/(sin((A2-alfa)*PI/180));
	double q2y=(h*sinA2)/(sin((A2-alfa)*PI/180));

	gp_Vec leftHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,-1));
	gp_Vec rightHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(0,0,1));

	gp_Pnt p1(0,0,0);
	gp_Pnt p2(p2y,0,p2x);

	gp_Pnt q1(bearingFlankHeightGap,0,-width);
	gp_Pnt q2(bearingFlankHeightGap+q2y,0,-(width+q2x));
	gp_Pnt r1((exhaustThickness*cosA1),0,-(exhaustThickness*sinA1));

	gp_Ax1 q1Axis(q1,gp_Dir(0,1,0));
	gp_Vec q1q2Vec=leftHorizontalVec.Rotated(q1Axis,(360-A2)*PI/180);
	q1q2Vec.Normalize();
	q1q2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A2-alfa)*PI/180)));
	gp_Pnt q1End=q1.Translated(q1q2Vec);
	q2=q1End;

	gp_Ax1 p1Axis(p1,gp_Dir(0,1,0));
	gp_Vec p1p2Vec=rightHorizontalVec.Rotated(q1Axis,A1*PI/180);
	p1p2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A1)*PI/180)));
	gp_Pnt p1End=p1.Translated(p1p2Vec);
	TopoDS_Edge q1q2Strait=BRepBuilderAPI_MakeEdge(p1End,p1);

	p2=p1End;


	TopoDS_Edge p1p2Edge=BRepBuilderAPI_MakeEdge(p1,p2);
	TopoDS_Edge p1q1Edge=BRepBuilderAPI_MakeEdge(p1,q1);
	TopoDS_Edge q1q2Edge=BRepBuilderAPI_MakeEdge(q1,q2);
	TopoDS_Edge p2q2Edge=BRepBuilderAPI_MakeEdge(p2,q2);
	TopoDS_Edge p1r1Edge=BRepBuilderAPI_MakeEdge(p1,r1);
	TopoDS_Edge q1r1Edge=BRepBuilderAPI_MakeEdge(q1,r1);

	TopoDS_Edge hirizontalLine=BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,-10),gp_Pnt(0,0,10));

	TopoDS_Vertex vp1=BRepBuilderAPI_MakeVertex(p1);
	TopoDS_Vertex vp2=BRepBuilderAPI_MakeVertex(p2);
	TopoDS_Vertex vq1=BRepBuilderAPI_MakeVertex(q1);
	TopoDS_Vertex vq2=BRepBuilderAPI_MakeVertex(q2);

	BRepBuilderAPI_MakeWire trape;
	trape.Add(p1p2Edge);
	trape.Add(p1q1Edge);
	trape.Add(q1q2Edge);
	trape.Add(p2q2Edge);

	TopoDS_Wire wire=trape;

	Standard_Real U1;
	Standard_Real U2;
	gp_Vec p1q1Vec;
	gp_Pnt p1q1Pnt;


	Handle_Geom_Curve p1q1LineCurve=BRep_Tool::Curve(p1q1Edge,U1,U2);


	Handle_Geom_Curve q1q2LineCurve=BRep_Tool::Curve(q1q2Edge,U1,U2);
	gp_Pnt q1q2pointZero;
	gp_Vec q1q2Vec1;
	gp_Vec q1q2Vec2;
	gp_Vec oZ(gp_Dir(0,1,0));
	q1q2LineCurve->D0(U1,q1q2pointZero);
	TopoDS_Vertex q1q2vertZero=BRepBuilderAPI_MakeVertex(q1q2pointZero);
	q1q2LineCurve->D2(U1,q1q2pointZero,q1q2Vec1,q1q2Vec2);
	gp_Vec q1q2crossed=q1q2Vec1.Crossed(oZ);
	q1q2crossed.Normalize();
	q1q2crossed.Multiply(bearingThickness);
	q1q2pointZero.Translate(q1q2crossed);
	q1q2vertZero=BRepBuilderAPI_MakeVertex(q1q2pointZero);
	TopoDS_Edge bearingSideBottomEdge=BRepBuilderAPI_MakeEdge(q1,q1q2pointZero);
	gp_Pnt q1q2EndPoint;
	q1q2LineCurve->D0(U2,q1q2EndPoint);
	double q1q2Length=U2;

	Handle_Geom_Curve p1p2LineCurve=BRep_Tool::Curve(p1p2Edge,U1,U2);
	gp_Pnt p1p2pointZero;
	gp_Vec p1p2Vec1;
	gp_Vec p1p2Vec2;
	p1p2LineCurve->D0(U1,p1p2pointZero);
	TopoDS_Vertex p1p2vertZero=BRepBuilderAPI_MakeVertex(p1p2pointZero);
	p1p2LineCurve->D2(U1,p1p2pointZero,p1p2Vec1,p1p2Vec2);
	gp_Vec p1p2Vec1Reversed=p1p2Vec1.Reversed();
	gp_Vec p1p2crossed=p1p2Vec1Reversed.Crossed(oZ);
	p1p2crossed.Normalize();
	p1p2crossed.Multiply(exhaustThickness);
	p1p2pointZero.Translate(p1p2crossed);
	p1p2vertZero=BRepBuilderAPI_MakeVertex(p1p2pointZero);
	TopoDS_Edge exhaustSideBottomEdge=BRepBuilderAPI_MakeEdge(p1,p1p2pointZero);
	gp_Pnt p1p2EndPoint;
	p1p2LineCurve->D0(U2,p1p2EndPoint);
	double p1p2Length=U2;

	gp_Pnt p1q1pointZero;
	gp_Vec p1q1Vec1;//parallel Vector to base 
	gp_Vec p1q1Vec2;

	p1q1LineCurve->D2(U1,p1q1pointZero,p1q1Vec1,p1q1Vec2);

	p1q1pointZero.Translate(p1q1Vec1);

	TopoDS_Vertex vertp1q1 =BRepBuilderAPI_MakeVertex(p1q1pointZero);

	//offset

	gp_Dir p1p2Direction(p1p2Vec1);
	gp_Lin p1p2ParallelLine(p1p2pointZero,p1p2Direction);
	TopoDS_Edge p1p2ParallelEdge=BRepBuilderAPI_MakeEdge(p1p2ParallelLine);

	gp_Dir q1q2Direction(q1q2Vec1);
	gp_Lin q1q2ParallelLine(q1q2pointZero,q1q2Direction);
	TopoDS_Edge q1q2ParallelEdge=BRepBuilderAPI_MakeEdge(q1q2ParallelLine);

	gp_Pnt bottomParallelEdgesIntersectionPnt=getMinimumDistancePoint(p1p2ParallelEdge,q1q2ParallelEdge);

	gp_Vec leftBottomVector=p1q1Vec1.Normalized();
	leftBottomVector.Multiply(tipRadius);
	gp_Pnt leftBottomPiontOfDividerWall=bottomParallelEdgesIntersectionPnt.Translated(leftBottomVector);

	gp_Vec rightBottomVector=leftBottomVector.Reversed();
	gp_Pnt rightBottomPointOfDividerWall=bottomParallelEdgesIntersectionPnt.Translated(rightBottomVector); 

	TopoDS_Vertex intersectionVert=BRepBuilderAPI_MakeVertex(bottomParallelEdgesIntersectionPnt);
	gp_Pnt dividerPointOnBase=getMinimumDistancePoint(p1q1Edge,intersectionVert);

	TopoDS_Edge leftVoluteFillEdge=BRepBuilderAPI_MakeEdge(leftBottomPiontOfDividerWall,dividerPointOnBase);
	TopoDS_Edge rightVoluteFillEdge=BRepBuilderAPI_MakeEdge(rightBottomPointOfDividerWall,dividerPointOnBase);
	TopoDS_Edge leftVoluteBaseEdge =BRepBuilderAPI_MakeEdge(dividerPointOnBase,q1);
	TopoDS_Edge rightVoluteBaseEdge=BRepBuilderAPI_MakeEdge(dividerPointOnBase,p1);

	TopoDS_Vertex LeftbottomVert=BRepBuilderAPI_MakeVertex(leftBottomPiontOfDividerWall);
	TopoDS_Vertex RightbottpmVert=BRepBuilderAPI_MakeVertex(rightBottomPointOfDividerWall);
	TopoDS_Vertex dividerOnBaseVert=BRepBuilderAPI_MakeVertex(dividerPointOnBase);

	gp_Pnt minimumDistancePointOnp1p2=getMinimumDistancePoint(p1p2Edge,RightbottpmVert);
	gp_Pnt minimumDistancePointOnq1q2=getMinimumDistancePoint(q1q2Edge,LeftbottomVert);

	TopoDS_Vertex mdVert=BRepBuilderAPI_MakeVertex(minimumDistancePointOnq1q2);

	TopoDS_Edge mdline=BRepBuilderAPI_MakeEdge(minimumDistancePointOnq1q2,leftBottomPiontOfDividerWall);
	TopoDS_Edge mdline1=BRepBuilderAPI_MakeEdge(minimumDistancePointOnp1p2,rightBottomPointOfDividerWall);


	double bearingSideOuterWallScalar=minimumDistancePointOnq1q2.Distance(q2);
	double exhaustSideOuterWallScalar=minimumDistancePointOnp1p2.Distance(p2);

	gp_Ax1 rightAx(rightBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec divideLineVec(dividerPointOnBase,bottomParallelEdgesIntersectionPnt);
	divideLineVec.Normalize();
	gp_Vec rightRotatedVec=divideLineVec.Rotated(rightAx,((360-(dividerAngle/2)))*PI/180);
	rightRotatedVec.Multiply(exhaustSideOuterWallScalar);
	gp_Pnt rightTopPointOfDividerWall=rightBottomPointOfDividerWall.Translated(rightRotatedVec);
	TopoDS_Edge rightDividerWall=BRepBuilderAPI_MakeEdge(rightTopPointOfDividerWall,rightBottomPointOfDividerWall);

	gp_Ax1 leftAx(rightBottomPointOfDividerWall,gp_Dir(0,1,0));
	gp_Vec leftRotateVec=divideLineVec.Rotated(leftAx,(dividerAngle/2)*PI/180); 
	leftRotateVec.Multiply(bearingSideOuterWallScalar);
	gp_Pnt leftTopPointofDividerWall=leftBottomPiontOfDividerWall.Translated(leftRotateVec);
	TopoDS_Edge leftDividerWall=BRepBuilderAPI_MakeEdge(leftTopPointofDividerWall,leftBottomPiontOfDividerWall);


	TopoDS_Edge topLeftBaseEdge=BRepBuilderAPI_MakeEdge(q1q2EndPoint,leftTopPointofDividerWall);
	Handle_Geom_Curve topLeftBaseLineCurve=BRep_Tool::Curve(topLeftBaseEdge,U1,U2);
	gp_Pnt leftCircleCentrePoint;
	topLeftBaseLineCurve->D0(U2/2,leftCircleCentrePoint);
	double topLeftCircleRadius=U2/2;
	gp_Ax2 leftCircleAxis(leftCircleCentrePoint,gp_Dir(0,1,0));
	gp_Circ TopLeftCircle(leftCircleAxis,topLeftCircleRadius);
	GC_MakeArcOfCircle topLeftHalfCircle(TopLeftCircle,leftTopPointofDividerWall,q1q2EndPoint,false);
	Handle(Geom_TrimmedCurve) leftCircleArc=topLeftHalfCircle;
	TopoDS_Edge leftCircle=BRepBuilderAPI_MakeEdge(leftCircleArc);
	TopoDS_Edge leftCircleArcBazierCurveEdge=convertTrimmToBezier(leftCircleArc,q1q2Vec1,leftRotateVec,bearingSideOuterWallScalar);
	leftCircle=leftCircleArcBazierCurveEdge;
	BRepTools::Write(leftCircle,"D:/Breps/newVolute/curve.brep");

	TopoDS_Edge topRightBaseEdge=BRepBuilderAPI_MakeEdge(p1p2EndPoint,rightTopPointOfDividerWall);
	Handle_Geom_Curve topRightBaseLineCurve=BRep_Tool::Curve(topRightBaseEdge,U1,U2);
	gp_Pnt rightCircleCentrePoint;
	topRightBaseLineCurve->D0(U2/2,rightCircleCentrePoint);
	double topRightCircleRaduis=U2/2;
	gp_Ax2 rightCircleAxis(rightCircleCentrePoint,gp_Dir(0,1,0));
	gp_Circ topRightCircle(rightCircleAxis,topRightCircleRaduis);
	GC_MakeArcOfCircle topRightHalfOfCircle(topRightCircle,p1p2EndPoint,rightTopPointOfDividerWall,false);
	Handle(Geom_TrimmedCurve) rightCircleArc=topRightHalfOfCircle;
	TopoDS_Edge rightCircle=BRepBuilderAPI_MakeEdge(rightCircleArc);
	TopoDS_Edge rightCircleArcBezierCurveEdge =convertTrimmToBezier(rightCircleArc,rightRotatedVec,p1p2Vec1,exhaustSideOuterWallScalar);
	rightCircle = rightCircleArcBezierCurveEdge;
	BRepTools::Write(rightCircle,"D:/Breps/newVolute/curve1.brep");

	gp_Ax2 bottomCircleAxis(bottomParallelEdgesIntersectionPnt,gp_Dir(0,1,0));
	gp_Circ tipcircle(bottomCircleAxis,tipRadius);
	GC_MakeArcOfCircle bottomHalfCircle(tipcircle,leftBottomPiontOfDividerWall,rightBottomPointOfDividerWall,false);
	Handle(Geom_TrimmedCurve) tipCircleArc=bottomHalfCircle;
	TopoDS_Edge bottomTipCircle=BRepBuilderAPI_MakeEdge(tipCircleArc);


	BRepBuilderAPI_MakeWire leftVoluteWire;
	leftVoluteWire.Add(leftVoluteBaseEdge);
	leftVoluteWire.Add(leftVoluteFillEdge);
	leftVoluteWire.Add(leftDividerWall);
	leftVoluteWire.Add(leftCircle);
	leftVoluteWire.Add(q1q2Edge);

	BRepBuilderAPI_MakeWire rightVoluteWire;
	rightVoluteWire.Add(rightVoluteBaseEdge);
	rightVoluteWire.Add(rightVoluteFillEdge);
	rightVoluteWire.Add(rightDividerWall);
	rightVoluteWire.Add(rightCircle);
	rightVoluteWire.Add(p1p2Edge);

	BRepTools::Write(leftVoluteWire,"D:/Breps/newVolute/leftVoluteWire.brep");

	w1=leftVoluteWire;
	w2=rightVoluteWire;
	p2q2e1=p2q2Edge;
	leftTopPoint=leftTopPointofDividerWall;
	rightTopPoint=rightTopPointOfDividerWall;
	leftBottomVec=leftBottomVector;
	dividerLineVec=divideLineVec;

}

gp_Pnt CDualVoluteDialog::getMinimumDistancePoint(TopoDS_Edge edge,TopoDS_Vertex vertex)
{

	BRepExtrema_DistShapeShape minimumDist(edge,vertex,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
	gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

	return divideIntersectionPointOnCurve;
}

gp_Pnt CDualVoluteDialog::getMinimumDistancePoint(TopoDS_Edge edge1,TopoDS_Edge edge2)
{

	TopoDS_Wire wire1=BRepBuilderAPI_MakeWire(edge1);
	TopoDS_Wire wire2=BRepBuilderAPI_MakeWire(edge2);


	BRepExtrema_DistShapeShape minimumDist(wire1,wire2,Extrema_ExtFlag_MIN,Extrema_ExtAlgo_Grad);
	gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

	return divideIntersectionPointOnCurve;
}

double  CDualVoluteDialog::getMaximumHeight(TopoDS_Wire section){

	TopoDS_Edge edges[10];
	TopExp_Explorer edgeExplorer(section,TopAbs_EDGE);
	int i=0;
	while(edgeExplorer.More()){
		TopoDS_Edge aEdge =TopoDS::Edge(edgeExplorer.Current());
		edges[i]=aEdge;
		edgeExplorer.Next();
		i++;

	}

	TopoDS_Vertex vertexArrayOfEdge0[2];
	TopExp_Explorer vertexExplorer(edges[0],TopAbs_VERTEX);
	int j=0;
	while(vertexExplorer.More()){
		TopoDS_Vertex aVertex=TopoDS::Vertex(vertexExplorer.Current());
		vertexArrayOfEdge0[j]=aVertex;
		vertexExplorer.Next();
		j++;
	}

	TopoDS_Vertex vertexArrayOfEdge8[2];
	TopExp_Explorer vertexExplorer8(edges[8],TopAbs_VERTEX);
	j=0;
	while(vertexExplorer8.More()){
		TopoDS_Vertex aVertex=TopoDS::Vertex(vertexExplorer8.Current());
		vertexArrayOfEdge8[j]=aVertex;
		vertexExplorer8.Next();
		j++;
	}


	BRepExtrema_ExtPC maxDistancePoint(vertexArrayOfEdge0[0],edges[2]);
	gp_Pnt maxPoint=maxDistancePoint.Point(1);
	TopoDS_Vertex maxVertex=BRepBuilderAPI_MakeVertex(maxPoint);
	BRepExtrema_DistShapeShape maxDistExtrema(edges[2],maxVertex);
	double value=maxDistExtrema.Value();
	gp_Pnt PointOnEdge2=BRep_Tool::Pnt(vertexArrayOfEdge0[0]);

	double distanceValue=PointOnEdge2.Distance(maxPoint);




	/*m_pcoloredshapeList->Add(Quantity_NOC_GREEN,section);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,maxVertex);
	m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,vertexArray[0]);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,edges[8]);
	m_pcoloredshapeList->Add(Quantity_NOC_RED,edges[0]);
	m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vertexArrayOfEdge0[0]);
	m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vertexArrayOfEdge8[0]);*/


	return distanceValue;
} 

void CDualVoluteDialog::getExitPipeTransitionCrossSection(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,  
														  double exitDividerWallWidth,double exitPipeRadius,gp_Vec divideLineVec)
{

	gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,1,0));
	gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
	TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


	gp_Vec straitExitDividerVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
	straitExitDividerVec.Rotate(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),90*PI/180);
	gp_Pnt zeroPnt(0,0,0);
	zeroPnt.Translate(straitExitDividerVec);


	gp_Dir exitDividerDir(divideLineVec);
	//gp_Dir exitDividerDir(straitExitDividerVec);

	gp_Vec leftWidthVectorOfExit=leftBottomVector.Normalized();
	leftWidthVectorOfExit.Multiply(exitDividerWallWidth/2);
	gp_Pnt leftWidthPointOfExitPipe=centrePointOfCircle.Translated(leftWidthVectorOfExit);

	gp_Vec rightWidthVectorOfExit=leftWidthVectorOfExit.Reversed();
	gp_Pnt rightWidthPointOfExitPipe=centrePointOfCircle.Translated(rightWidthVectorOfExit);

	GeomConvert convert;
	gp_Lin leftExitDividerLine(leftWidthPointOfExitPipe,exitDividerDir);
	TopoDS_Edge leftExitDividerEdge = BRepBuilderAPI_MakeEdge(leftExitDividerLine);
	TopoDS_Wire circleWire = BRepBuilderAPI_MakeWire(exitCircleEdge);
	TopoDS_Wire leftdividerWire = BRepBuilderAPI_MakeWire(leftExitDividerEdge);
	BRepExtrema_DistShapeShape LeftminimumDist(circleWire,leftdividerWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
	gp_Pnt topLeftExitPnt = LeftminimumDist.PointOnShape1(1);
	gp_Pnt bottomLeftExitPnt = LeftminimumDist.PointOnShape1(2);
	leftExitDividerEdge = BRepBuilderAPI_MakeEdge(topLeftExitPnt,bottomLeftExitPnt);
	GC_MakeArcOfCircle leftExitCircleArc(exitPipeCircle,topLeftExitPnt,bottomLeftExitPnt,false);
	Handle_Geom_TrimmedCurve leftExitCircleArcTrimmedCurve = leftExitCircleArc;
	TopoDS_Edge leftExitArcEdge=BRepBuilderAPI_MakeEdge(leftExitCircleArcTrimmedCurve);

	TopoDS_Vertex topVert=BRepBuilderAPI_MakeVertex(topLeftExitPnt);
	TopoDS_Vertex bottomVert=BRepBuilderAPI_MakeVertex(bottomLeftExitPnt);

	double filletRadius=8;

	gp_Pln filletPlane(centrePointOfCircle,gp_Dir(0,1,0));
	ChFi2d_FilletAPI leftTopFillet(leftExitDividerEdge,leftExitArcEdge,filletPlane);
	leftTopFillet.Perform(filletRadius);
	TopoDS_Edge leftTopFilletEdge=leftTopFillet.Result(topLeftExitPnt,leftExitDividerEdge,leftExitArcEdge,-1);

	ChFi2d_FilletAPI leftBottomFillet(leftExitArcEdge,leftExitDividerEdge,filletPlane);
	leftBottomFillet.Perform(filletRadius);
	TopoDS_Edge leftBottomFilletEdge=leftBottomFillet.Result(bottomLeftExitPnt,leftExitArcEdge,leftExitDividerEdge,-1);

	gp_Lin rightExitDividerLine(rightWidthPointOfExitPipe,exitDividerDir);
	TopoDS_Edge rightExitDividerEdge=BRepBuilderAPI_MakeEdge(rightExitDividerLine);
	TopoDS_Wire rightDividerWire=BRepBuilderAPI_MakeWire(rightExitDividerEdge);
	BRepExtrema_DistShapeShape rightMinimumDist(circleWire,rightDividerWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
	gp_Pnt topRightExitPnt= rightMinimumDist.PointOnShape1(1);
	gp_Pnt bottomRightExitPnt= rightMinimumDist.PointOnShape1(2);
	rightExitDividerEdge=BRepBuilderAPI_MakeEdge(topRightExitPnt,bottomRightExitPnt);
	GC_MakeArcOfCircle rightExitCircleArc(exitPipeCircle,bottomRightExitPnt,topRightExitPnt,false);
	Handle_Geom_TrimmedCurve rightExitCircleArcTrimmedCurve = rightExitCircleArc;
	TopoDS_Edge rightExitArcEdge = BRepBuilderAPI_MakeEdge(rightExitCircleArcTrimmedCurve);

	ChFi2d_FilletAPI rightTopFillet(rightExitArcEdge,rightExitDividerEdge,filletPlane);
	rightTopFillet.Perform(filletRadius);
	TopoDS_Edge rightTopFilletEdge=rightTopFillet.Result(topRightExitPnt,rightExitArcEdge,rightExitDividerEdge,-1);

	TopoDS_Vertex ToprightFilletVertex[2];
	TopExp_Explorer vertexExplorer(rightTopFilletEdge,TopAbs_VERTEX);
	int j=0;
	while(vertexExplorer.More()){
		TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
		ToprightFilletVertex[j]=aVertex;
		vertexExplorer.Next();
		j++;

	}

	TopoDS_Vertex bottomRightExitVert=BRepBuilderAPI_MakeVertex(bottomRightExitPnt);
	TopoDS_Edge rightExitDividerAfterFillet=BRepBuilderAPI_MakeEdge(ToprightFilletVertex[1],bottomRightExitVert);

	ChFi2d_FilletAPI rightBottomFillet(rightExitDividerAfterFillet,rightExitArcEdge,filletPlane);
	rightBottomFillet.Perform(filletRadius);
	TopoDS_Edge rightBottomFilletEdge=rightBottomFillet.Result(bottomRightExitPnt,rightExitDividerAfterFillet,rightExitArcEdge,-1);


	TopoDS_Vertex BottomrightFilletVertex[2];
	TopExp_Explorer vertexExplorer1(rightBottomFilletEdge,TopAbs_VERTEX);
	j=0;
	while(vertexExplorer1.More()){
		TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer1.Current());
		BottomrightFilletVertex[j]=aVertex;
		vertexExplorer1.Next();
		j++;

	}
	rightExitArcEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle,BottomrightFilletVertex[1],ToprightFilletVertex[0]);

	BRepBuilderAPI_MakeWire leftPartOfExitPipe;
	leftPartOfExitPipe.Add(leftTopFilletEdge);
	leftPartOfExitPipe.Add(leftExitDividerEdge);
	leftPartOfExitPipe.Add(leftBottomFilletEdge);
	leftPartOfExitPipe.Add(leftExitArcEdge);

	w1=leftPartOfExitPipe;

	BRepBuilderAPI_MakeWire rightPartOfExitPipe;
	rightPartOfExitPipe.Add(rightTopFilletEdge);
	rightPartOfExitPipe.Add(rightExitDividerAfterFillet);
	rightPartOfExitPipe.Add(rightBottomFilletEdge);
	rightPartOfExitPipe.Add(rightExitArcEdge);
	w2=rightPartOfExitPipe;

		TopoDS_Compound compound;
	BRep_Builder builder;
	builder.MakeCompound(compound);
	builder.Add(compound,w1);
	builder.Add(compound, w2);

	BRepTools::Write(compound,"D:/Breps/newVolute/transitionCom.brep");


}


double CDualVoluteDialog::getDividerWallMaximumWidth(TopoDS_Wire section)
{

	TopoDS_Edge edges[10];
	TopExp_Explorer edgeExplorer(section,TopAbs_EDGE);
	int i=0;
	while(edgeExplorer.More()){
		TopoDS_Edge aEdge =TopoDS::Edge(edgeExplorer.Current());
		edges[i]=aEdge;
		edgeExplorer.Next();
		i++;


	}
	return 0;
}

void CDualVoluteDialog::getExitPipeEnding(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,
										  double exitDividerWallWidth,double exitPipeRadius,gp_Vec divideLineVec)
{

	gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,1,0));
	gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
	TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


	gp_Vec straitExitDividerVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
	straitExitDividerVec.Rotate(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),90*PI/180);
	gp_Pnt zeroPnt(0,0,0);
	zeroPnt.Translate(straitExitDividerVec);


	gp_Dir exitDividerDir(divideLineVec);
	//gp_Dir exitDividerDir(straitExitDividerVec);

	gp_Vec leftWidthVectorOfExit=leftBottomVector.Normalized();
	leftWidthVectorOfExit.Multiply(exitDividerWallWidth/2);
	gp_Pnt leftWidthPointOfExitPipe=centrePointOfCircle.Translated(leftWidthVectorOfExit);

	gp_Vec rightWidthVectorOfExit=leftWidthVectorOfExit.Reversed();
	gp_Pnt rightWidthPointOfExitPipe=centrePointOfCircle.Translated(rightWidthVectorOfExit);

	GeomConvert convert;
	gp_Lin leftExitDividerLine(leftWidthPointOfExitPipe,exitDividerDir);
	TopoDS_Edge leftExitDividerEdge = BRepBuilderAPI_MakeEdge(leftExitDividerLine);
	TopoDS_Wire circleWire = BRepBuilderAPI_MakeWire(exitCircleEdge);
	TopoDS_Wire leftdividerWire = BRepBuilderAPI_MakeWire(leftExitDividerEdge);
	BRepExtrema_DistShapeShape LeftminimumDist(circleWire,leftdividerWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
	gp_Pnt topLeftExitPnt = LeftminimumDist.PointOnShape1(1);
	gp_Pnt bottomLeftExitPnt = LeftminimumDist.PointOnShape1(2);
	leftExitDividerEdge = BRepBuilderAPI_MakeEdge(topLeftExitPnt,bottomLeftExitPnt);
	GC_MakeArcOfCircle leftExitCircleArc(exitPipeCircle,topLeftExitPnt,bottomLeftExitPnt,false);
	Handle_Geom_TrimmedCurve leftExitCircleArcTrimmedCurve = leftExitCircleArc;
	TopoDS_Edge leftExitArcEdge=BRepBuilderAPI_MakeEdge(leftExitCircleArcTrimmedCurve);

	TopoDS_Vertex topVert=BRepBuilderAPI_MakeVertex(topLeftExitPnt);
	TopoDS_Vertex bottomVert=BRepBuilderAPI_MakeVertex(bottomLeftExitPnt);

	double filletRadius=8;

	gp_Pln filletPlane(centrePointOfCircle,gp_Dir(0,1,0));
	ChFi2d_FilletAPI leftTopFillet(leftExitDividerEdge,leftExitArcEdge,filletPlane);
	leftTopFillet.Perform(filletRadius);
	TopoDS_Edge leftTopFilletEdge=leftTopFillet.Result(topLeftExitPnt,leftExitDividerEdge,leftExitArcEdge,-1);

	ChFi2d_FilletAPI leftBottomFillet(leftExitArcEdge,leftExitDividerEdge,filletPlane);
	leftBottomFillet.Perform(filletRadius);
	TopoDS_Edge leftBottomFilletEdge=leftBottomFillet.Result(bottomLeftExitPnt,leftExitArcEdge,leftExitDividerEdge,-1);

	gp_Lin rightExitDividerLine(rightWidthPointOfExitPipe,exitDividerDir);
	TopoDS_Edge rightExitDividerEdge=BRepBuilderAPI_MakeEdge(rightExitDividerLine);
	TopoDS_Wire rightDividerWire=BRepBuilderAPI_MakeWire(rightExitDividerEdge);
	BRepExtrema_DistShapeShape rightMinimumDist(circleWire,rightDividerWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
	gp_Pnt topRightExitPnt= rightMinimumDist.PointOnShape1(1);
	gp_Pnt bottomRightExitPnt= rightMinimumDist.PointOnShape1(2);
	rightExitDividerEdge=BRepBuilderAPI_MakeEdge(topRightExitPnt,bottomRightExitPnt);
	GC_MakeArcOfCircle rightExitCircleArc(exitPipeCircle,bottomRightExitPnt,topRightExitPnt,false);
	Handle_Geom_TrimmedCurve rightExitCircleArcTrimmedCurve = rightExitCircleArc;
	TopoDS_Edge rightExitArcEdge = BRepBuilderAPI_MakeEdge(rightExitCircleArcTrimmedCurve);

	ChFi2d_FilletAPI rightTopFillet(rightExitArcEdge,rightExitDividerEdge,filletPlane);
	rightTopFillet.Perform(filletRadius);
	TopoDS_Edge rightTopFilletEdge=rightTopFillet.Result(topRightExitPnt,rightExitArcEdge,rightExitDividerEdge,-1);

	TopoDS_Vertex ToprightFilletVertex[2];
	TopExp_Explorer vertexExplorer(rightTopFilletEdge,TopAbs_VERTEX);
	int j=0;
	while(vertexExplorer.More()){
		TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
		ToprightFilletVertex[j]=aVertex;
		vertexExplorer.Next();
		j++;

	}

	TopoDS_Vertex bottomRightExitVert=BRepBuilderAPI_MakeVertex(bottomRightExitPnt);
	TopoDS_Edge rightExitDividerAfterFillet=BRepBuilderAPI_MakeEdge(ToprightFilletVertex[1],bottomRightExitVert);

	ChFi2d_FilletAPI rightBottomFillet(rightExitDividerAfterFillet,rightExitArcEdge,filletPlane);
	rightBottomFillet.Perform(filletRadius);
	TopoDS_Edge rightBottomFilletEdge=rightBottomFillet.Result(bottomRightExitPnt,rightExitDividerAfterFillet,rightExitArcEdge,-1);


	TopoDS_Vertex BottomrightFilletVertex[2];
	TopExp_Explorer vertexExplorer1(rightBottomFilletEdge,TopAbs_VERTEX);
	j=0;
	while(vertexExplorer1.More()){
		TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer1.Current());
		BottomrightFilletVertex[j]=aVertex;
		vertexExplorer1.Next();
		j++;

	}
	rightExitArcEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle,BottomrightFilletVertex[1],ToprightFilletVertex[0]);

	BRepBuilderAPI_MakeWire leftPartOfExitPipe;
	leftPartOfExitPipe.Add(leftTopFilletEdge);
	leftPartOfExitPipe.Add(leftExitDividerEdge);
	leftPartOfExitPipe.Add(leftBottomFilletEdge);
	leftPartOfExitPipe.Add(leftExitArcEdge);

	w1=leftPartOfExitPipe;

	BRepBuilderAPI_MakeWire rightPartOfExitPipe;
	rightPartOfExitPipe.Add(rightTopFilletEdge);
	rightPartOfExitPipe.Add(rightExitDividerAfterFillet);
	rightPartOfExitPipe.Add(rightBottomFilletEdge);
	rightPartOfExitPipe.Add(rightExitArcEdge);
	w2=rightPartOfExitPipe;

}


TopoDS_Edge CDualVoluteDialog::convertTrimmToBezier(Handle_Geom_Curve curve,gp_Vec leftVec,gp_Vec rightVec,double scaler)
{
	GeomConvert convert;
	Handle_Geom_BSplineCurve ArcBSplineCurve=convert.CurveToBSplineCurve(curve);
	GeomConvert_BSplineCurveToBezierCurve bsplineToBezierCurve(ArcBSplineCurve);
	Handle_Geom_BezierCurve Arc1BezierCurve=bsplineToBezierCurve.Arc(1);
	Handle_Geom_BezierCurve Arc2BezierCurve=bsplineToBezierCurve.Arc(2);

	gp_Pnt startPointOfBSplineCurve=ArcBSplineCurve->StartPoint();
	gp_Pnt endPointOfBSplineCurve = ArcBSplineCurve->EndPoint();

	double lengthOfCurve =startPointOfBSplineCurve.Distance(endPointOfBSplineCurve);


	leftVec.Normalize();
	rightVec.Normalize();

	leftVec.Multiply(scaler/10);
	rightVec.Multiply(scaler/10);

	TopoDS_Edge curveEdge=BRepBuilderAPI_MakeEdge( ArcBSplineCurve);
	TopoDS_Edge e1=BRepBuilderAPI_MakeEdge(Arc1BezierCurve);


	gp_Pnt pointArray1[4];
	TopoDS_Vertex vertArray[4];
	for(int i=0;i<Arc1BezierCurve->NbPoles();i++){

		pointArray1[i]=Arc1BezierCurve->Pole(i+1);
		vertArray[i]=BRepBuilderAPI_MakeVertex(pointArray1[i]);

	}	
	gp_Pnt startPoint=pointArray1[0];
	gp_Pnt leftMiddlePoint=pointArray1[1];
	gp_Pnt leftTanPoint=pointArray1[0].Translated(leftVec);

	Arc1BezierCurve->InsertPoleAfter(1,leftTanPoint);
	e1=BRepBuilderAPI_MakeEdge(Arc1BezierCurve);

	gp_Pnt pointArray2[4];
	for(int i=0;i<Arc2BezierCurve->NbPoles();i++){

		pointArray2[i]=Arc2BezierCurve->Pole(i+1);
		vertArray[i]=BRepBuilderAPI_MakeVertex(pointArray2[i]);

	}	
	gp_Pnt middlePoint= pointArray2[0];
	gp_Pnt endPoint =pointArray2[2];
	gp_Pnt rightMiddlePoint=pointArray2[1];

	gp_Pnt rightTanPoint=pointArray2[2].Translated(rightVec);

	Arc2BezierCurve->InsertPoleBefore(3,rightTanPoint);
	TopoDS_Edge e2=BRepBuilderAPI_MakeEdge(Arc2BezierCurve);
	TopoDS_Wire curveWire=BRepBuilderAPI_MakeWire(e1,e2);  

	TopoDS_Edge mergedEdge;
	mergedEdge=MergeEdges(e1,e2);

	return mergedEdge;
}

void CDualVoluteDialog::createTransitionExitPipePart(TopoDS_Shape& s1,TopoDS_Shape& s2,TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt& centre,TopoDS_Wire crossSection27,TopoDS_Wire leftVoluteWire, TopoDS_Wire rightVoluteWire ,TopoDS_Edge p2q2Edge,gp_Vec leftBottomVector ,gp_Vec divideLineVec , double exitPipeRadius,double exitDividerAngle,double maxWidthOfDividerWall,double transitionPartLength, double voluteRadius ){

	Standard_Real U1;
	Standard_Real U2;

	Handle_Geom_Curve p2q2LineCurve=BRep_Tool::Curve(p2q2Edge,U1,U2);
	gp_Pnt centrePointOfCircle;
	p2q2LineCurve->D0(U2/2,centrePointOfCircle);

	TopoDS_Vertex p2q2Vert=BRepBuilderAPI_MakeVertex(centrePointOfCircle); 

	gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,1,0));
	gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
	TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


	TopoDS_Wire leftPartOfExitPipeTransition;
	TopoDS_Wire rightPartOfExitPipeTransition;

	getExitPipeTransitionCrossSection(leftPartOfExitPipeTransition,rightPartOfExitPipeTransition,centrePointOfCircle,exitDividerAngle,leftBottomVector,maxWidthOfDividerWall,exitPipeRadius,divideLineVec);

	//Transition
	double tRLength=getMaximumHeight(crossSection27);
	transitionPartLength=-voluteRadius-tRLength;
	gp_Vec exitPipeTranslateVec(gp_Dir(0,1,0));
	exitPipeTranslateVec.Multiply(transitionPartLength);
	gp_Trsf translate;
	translate.SetTranslation(exitPipeTranslateVec);

	BRepBuilderAPI_Transform leftTransform( leftPartOfExitPipeTransition,translate,true);
	TopoDS_Wire leftTransition= TopoDS::Wire(leftTransform.ModifiedShape( leftPartOfExitPipeTransition));

	BRepBuilderAPI_Transform rightTransform(rightPartOfExitPipeTransition,translate,true);
	TopoDS_Wire rightTransition=TopoDS::Wire(rightTransform.ModifiedShape(rightPartOfExitPipeTransition));

	BRepTools::Write(leftVoluteWire,"D:/Breps/newVolute/leftVoluteWire.brep");
	BRepTools::Write(rightVoluteWire,"D:/Breps/newVolute/rightVoluteWire.brep");
	BRepTools::Write(leftTransition,"D:/Breps/newVolute/lefttrans.brep");
	BRepTools::Write(rightTransition,"D:/Breps/newVolute/righttrans.brep");

	m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftTransition);


		m_pcoloredshapeList->Display(myAISContext);

	BRepOffsetAPI_ThruSections leftTransitionThruSection;
	leftTransitionThruSection.AddWire(leftVoluteWire);
	leftTransitionThruSection.AddWire(leftTransition);
	leftTransitionThruSection.Build();


	BRepOffsetAPI_ThruSections rightTransitionThruSection;
	rightTransitionThruSection.AddWire(rightVoluteWire);
	rightTransitionThruSection.AddWire(rightTransition);
	rightTransitionThruSection.Build();

	TopoDS_Shape leftTransitionPartShape=leftTransitionThruSection.Shape();
	TopoDS_Shape rightTransitionPartShape=rightTransitionThruSection.Shape();

	s1=leftTransitionPartShape;
	s2=rightTransitionPartShape;
	w1=leftTransition;
	w2=rightTransition;
	centre=centrePointOfCircle;

}

void CDualVoluteDialog::createExitPipeEndigPart(TopoDS_Shape& s1,TopoDS_Shape& s2,TopoDS_Wire leftTransition,TopoDS_Wire rightTransition,gp_Pnt centrePointOfCircle,gp_Vec leftBottomVector,gp_Vec divideLineVec ,double exitDividerWallWidth,double exitPipeRadius,double exitDividerAngle,double transitionPartLength,double exitPipeLength)
{

	gp_Vec exitPipeTranslateVec1(gp_Dir(0,1,0));
	gp_Trsf translate;
	exitPipeTranslateVec1.Multiply(transitionPartLength+exitPipeLength);
	translate.SetTranslation(exitPipeTranslateVec1);

	TopoDS_Wire leftPartOfExitPipeEnding;
	TopoDS_Wire rightPartOfExitPipeEnding;

	getExitPipeEnding(leftPartOfExitPipeEnding,rightPartOfExitPipeEnding,centrePointOfCircle,exitDividerAngle,leftBottomVector,exitDividerWallWidth,exitPipeRadius,divideLineVec);

	BRepBuilderAPI_Transform leftExitPipeTrsf(leftPartOfExitPipeEnding,translate,true);
	TopoDS_Wire endOfLeftExitPipe=TopoDS::Wire(leftExitPipeTrsf.ModifiedShape(leftPartOfExitPipeEnding));

	BRepBuilderAPI_Transform rightExitPipeTrsf(rightPartOfExitPipeEnding,translate,true);
	TopoDS_Wire endOfRightExitPipe=TopoDS::Wire(rightExitPipeTrsf.ModifiedShape(rightPartOfExitPipeEnding));

	BRepOffsetAPI_ThruSections leftExitPipeSections;
	leftExitPipeSections.AddWire(leftTransition);
	leftExitPipeSections.AddWire(endOfLeftExitPipe);
	leftExitPipeSections.Build();

	BRepOffsetAPI_ThruSections rightExitPipeSection;
	rightExitPipeSection.AddWire(rightTransition);
	rightExitPipeSection.AddWire(endOfRightExitPipe);
	rightExitPipeSection.Build();

	TopoDS_Compound compound;
	BRep_Builder builder;
	builder.MakeCompound(compound);
	builder.Add(compound,leftPartOfExitPipeEnding);
	builder.Add(compound, rightPartOfExitPipeEnding);

	BRepTools::Write(compound,"D:/Breps/newVolute/endingCom.brep");

	TopoDS_Shape leftExitPipeEnding=leftExitPipeSections.Shape();
	TopoDS_Shape rightExitPipeEnding=rightExitPipeSection.Shape();

	s1=leftExitPipeEnding;
	s2=rightExitPipeEnding;

}


	TopoDS_Edge CDualVoluteDialog::MergeEdges(TopoDS_Edge edge1,TopoDS_Edge edge2)
{
       Handle_Geom_Curve pGeoCurve;
       BRep_Builder BrepBuilder;
       TopoDS_Compound topoNewCurves;
       Handle_AIS_Shape pNewAISCurves;
       double first, last;
       double tolerance=0.001;
      
              BrepBuilder.MakeCompound(topoNewCurves);
            bool suceed = true;
              bool firstedge = true;
              TopoDS_Edge edge =edge1;
              pGeoCurve = BRep_Tool::Curve(edge, first, last);
              Handle(Geom_TrimmedCurve) curve = new Geom_TrimmedCurve(pGeoCurve, pGeoCurve->FirstParameter(), pGeoCurve->LastParameter());
              GeomConvert_CompCurveToBSplineCurve final_spline(curve);
            
                                  edge = edge2;
                                  pGeoCurve = BRep_Tool::Curve(edge, first, last);
                                  curve = new Geom_TrimmedCurve(pGeoCurve, pGeoCurve->FirstParameter(), pGeoCurve->LastParameter());
								  if (!final_spline.Add(curve, tolerance, Standard_False))
                                  {
                                         suceed = false;
                                                                    
              }

				TopoDS_Edge merged;
              if (suceed)
              {
                     BRepBuilderAPI_MakeEdge theEdgeBuilder_final(final_spline.BSplineCurve());
                     BrepBuilder.Add(topoNewCurves, theEdgeBuilder_final.Edge());
					 BRepTools::Write(theEdgeBuilder_final,"D:/Breps/newVolute/merged.brep");
					 merged=theEdgeBuilder_final;
              }
       
			  return merged;

}







void CDualVoluteDialog::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	CDialog::OnCancel();
}


