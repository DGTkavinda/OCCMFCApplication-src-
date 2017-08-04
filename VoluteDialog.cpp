// VoluteDialog.cpp : implementation file
//

#include "stdafx.h"
#include "VoluteDialog.h"
#include "afxdialogex.h"
#include "AIS_InteractiveObject.hxx"
#include "AIS_Shape.hxx"

#define PI 3.14159265
// CVoluteDialog dialog

IMPLEMENT_DYNAMIC(CVoluteDialog, CDialog)

	CVoluteDialog::CVoluteDialog(Handle_AIS_InteractiveContext myContext,CWnd* pParent /*=NULL*/)
	: CDialog(CVoluteDialog::IDD, pParent)
{
	myAISContext = myContext;
	myAISContext->SetDisplayMode(AIS_Shaded);
}

CVoluteDialog::~CVoluteDialog()
{
}

void CVoluteDialog::DoDataExchange(CDataExchange* pDX)
{

	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX,IDC_EDIT_H,m_edtHeight);
	DDX_Control(pDX,IDC_EDIT_W,m_edtWidth);
	DDX_Control(pDX,IDC_EDIT_A1,m_edtAngle1);
	DDX_Control(pDX,IDC_EDIT_A2,m_edtAngle2);
	DDX_Control(pDX,IDC_EDIT_P,m_edtPercentage);
	DDX_Control(pDX,IDC_EDIT_Fr,m_edtFraction);
	DDX_Control(pDX,IDC_EDIT_L,m_edtLength);
	DDX_Control(pDX,IDC_EDIT_Inlet,m_edtInlet);


}


BEGIN_MESSAGE_MAP(CVoluteDialog, CDialog)
	ON_BN_CLICKED(1, &CVoluteDialog::OnBnClicked1)
END_MESSAGE_MAP()


// CVoluteDialog message handlers


void CVoluteDialog::OnBnClicked1()
{

	double angle1=45.00;
	double angle2=45.00;
	double height=25+5;
	double initialHeight=height;
	double width=10;
	double startPointx=0;
	double ang1Tan= tan(angle1*PI/180.0);
	double ang2Tan= tan(angle2*PI/180.0);
	double area=150;
	double degree_0_percentage=5;
	double exitAreaFraction;
	double exitPipeLength;
	double inletRaduis;
	double newWireWholeArea;


	angle1=45.00;
	angle2=45.00;
	height=10;
	initialHeight=height;
	width=5;
	exitAreaFraction=2;
	exitPipeLength=20;
	inletRaduis=15;

	startPointx=0;
	ang1Tan=tan(angle1*PI/180.0);
	ang2Tan=tan(angle2*PI/180.0);
	degree_0_percentage=20;
	

	CString str;

	m_edtHeight.GetWindowText(str);
	height=_wtof(str);

	m_edtWidth.GetWindowText(str);
	width=_wtof(str);

	m_edtAngle1.GetWindowText(str);
	angle1=_wtof(str);

	m_edtAngle2.GetWindowText(str);
	angle2=_wtof(str);

	m_edtPercentage.GetWindowText(str);
	degree_0_percentage=_wtof(str);

	m_edtFraction.GetWindowText(str);
	exitAreaFraction=_wtof(str);

	m_edtLength.GetWindowText(str);
	exitPipeLength=_wtof(str);

	m_edtInlet.GetWindowText(str);
	inletRaduis=_wtof(str);


	if(height<=0||width<=0||angle1<=0||angle2<=0||degree_0_percentage<=0||exitAreaFraction<=0||exitPipeLength<=0){

		AfxMessageBox(_T("Values cannot be minus or null!"));

	} else if(angle1>80||angle2>80){
		AfxMessageBox(_T("Angles cannot Exceed 80 Degrees!"));

	}else{

		ang1Tan=tan(angle1*PI/180.0);
		ang2Tan=tan(angle2*PI/180.0);


		gp_Pnt r1(width/2.,height*2.,0);
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);
		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);
		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);
		gp_Pnt p3(width+((height+height/10)/ang1Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang1Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang1Tan),(height+height/3),0);

		BRepBuilderAPI_MakeEdge line1(p1,p2);
		BRepBuilderAPI_MakeEdge line2(q1,q2);
		BRepBuilderAPI_MakeEdge line3(q1,p1);
		BRepBuilderAPI_MakeEdge line4(q2,p2);

		BRepBuilderAPI_MakeWire voluteBaseWire1(line1,line3,line2);
		BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2,line4);

		//----------------------------------

		TColgp_Array1OfPnt qCurvePoles(1,5); 
		qCurvePoles(1) = q2; 
		qCurvePoles(2) = q4; 
		qCurvePoles(3) = r1;
		qCurvePoles(4) = p4; 
		qCurvePoles(5) = p2; 
		Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
		TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 
		BRepBuilderAPI_MakeWire initialCurveWire(curveEdge,line4);
		BRepBuilderAPI_MakeWire initialwholeWire (voluteBaseWire1,curveEdge);

		BRepTools::Write( initialwholeWire,"D:/Breps/tRapAxis.brep");
		//BRepBuilderAPI_MakeWire wire2(wire1,mainE);
		//BRepBuilderAPI_MakeWire voluteCross(wire2,pE)
		TopoDS_Shape sh=trapeziumWire;
		//Fit();

		//for 50% piece
		double initialTrapaziumArea = getSurfaceArea(trapeziumWire);
		double initailWholeArea= getSurfaceArea(initialwholeWire);
		double expectedWholeArea = initailWholeArea*(degree_0_percentage/100);//Expected whole area

		double intialCurveArea=getSurfaceArea(initialCurveWire);

		
		double ExpectedTrapeziumArea = initialTrapaziumArea*(degree_0_percentage/100);//user percentage of the 0 point area
		double newTrapeziumHeight= getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);//calculating the height
		double ExpectedCurveArea = intialCurveArea*(degree_0_percentage/100);
		TopoDS_Wire newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);

		newWireWholeArea=getSurfaceArea(newWholeWire);

		double degree_0_area=newWireWholeArea;

		TopoDS_Shape rotatedShape;
		gp_Trsf transfer;
		BRepOffsetAPI_ThruSections sections;

		gp_Ax1 axis(gp_Pnt(0,-inletRaduis,0),gp_Dir(1,0,0));
		//initial wire rotation

		double ang=1;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated(newWholeWire,transfer);
		rotatedShape=rotated.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);

	


		double areaIncreasingFactor=(100-degree_0_percentage)/18;

		CString areaIncreacingFactor;

		for(int i=2;i<=18;i++){

			double ang=20*(i-1);
			transfer.SetRotation(axis,ang*PI/180);

			double i1=i;
			double expectedWholeArea = initailWholeArea*((areaIncreasingFactor/100)*i);
			double ExpectedTrapeziumArea = initialTrapaziumArea*((areaIncreasingFactor/100)*i);
			double newTrapeziumHeight=getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);
			newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);


			BRepBuilderAPI_Transform rotated(newWholeWire,transfer);

			rotatedShape=rotated.Shape();
			newWholeWire=TopoDS::Wire(rotatedShape);
			sections.AddWire(newWholeWire);

			double newWholArea=getSurfaceArea(newWholeWire);
		CString initialWholeAreaString;
		CString expectedWholeAreaString;
		CString newWholeAreaString;
		CString msg;
		initialWholeAreaString.Format(_T("initialWholeAreaString %g \n" ),initailWholeArea);
		expectedWholeAreaString.Format(_T("expectedWholeAreaString %g \n"),expectedWholeArea);
		newWholeAreaString.Format(_T("newWholeArea %g \n"),i1);
		
		msg=initialWholeAreaString+expectedWholeAreaString+newWholeAreaString;
		//AfxMessageBox(msg);


		}


		ang=360;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated1(initialwholeWire,transfer);
		rotatedShape=rotated1.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);
		sections.Build();

		TopoDS_Shape completeVolute=sections.Shape();

		//circle
		gp_Dir dir(0,0,1); 
		double exitPipeRadius=sqrt((initailWholeArea*exitAreaFraction)/PI);
		gp_Pnt point(width/2,height*0.7,exitPipeLength);
		gp_Circ circle1(gp_Ax2( point, dir), exitPipeRadius);
		BRepBuilderAPI_MakeEdge circle(circle1);
		BRepBuilderAPI_MakeWire exitPipeCircleWire(circle);
		BRepOffsetAPI_ThruSections exitSection;

		exitSection.AddWire(initialwholeWire);
		exitSection.AddWire(exitPipeCircleWire);
		exitSection.Build();
		TopoDS_Shape exitPipe=exitSection.Shape();

		TopoDS_Compound completeShape;
		BRep_Builder aBuilder;
		aBuilder.MakeCompound (completeShape);
		aBuilder.Add(completeShape,completeVolute);
		aBuilder.Add(completeShape,exitPipe);

		TopoDS_Shape voluteAndExit=completeShape;

		BRepTools::Write(voluteAndExit,"D:/Breps/whoaleShapewithExit.brep");	
		
		Handle(AIS_InteractiveObject) obj1 = (new AIS_Shape(voluteAndExit));

		myAISContext->RemoveAll();
		myAISContext->Display(obj1);	
		OCC_3dBaseDoc::Fit();

		double newWholArea=getSurfaceArea(newWholeWire);
		CString initialWholeAreaString;
		CString expectedWholeAreaString;
		CString newWholeAreaString;
		CString msg;
		initialWholeAreaString.Format(_T("initialWholeAreaString %g \n" ),initailWholeArea);
		expectedWholeAreaString.Format(_T("expectedWholeAreaString %g \n"),expectedWholeArea);
		newWholeAreaString.Format(_T("newWholeArea %g \n"),newWireWholeArea);
		
		msg=initialWholeAreaString+expectedWholeAreaString+newWholeAreaString;
		AfxMessageBox(msg);

	}

}

TopoDS_Wire CVoluteDialog:: createNewShapeWithRightArea(double height,double width,double ang1Tan,double ang2Tan,double expectedArea,double initialHieght){


	double variationTolerance=expectedArea*0.001;
	double heightVariation=height*0.5;
	double r1Height;

	r1Height=height*2;
	gp_Pnt r1(width/2.,height*2.,0);
	gp_Pnt q1(0,0,0); 
	gp_Pnt q2(-(height/ang2Tan),height,0);
	gp_Pnt p1(width,0,0);
	gp_Pnt p2(((height/ang1Tan)+width),height,0);
	gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
	gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
	gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);
	gp_Pnt p3(width+((height+height/10)/ang1Tan),(height+height/10),0);
	gp_Pnt p4(width+((height+height/4)/ang1Tan),(height+height/4),0);
	gp_Pnt p5((width+height/ang1Tan),(height+height/3),0);

	BRepBuilderAPI_MakeEdge line1(p1,p2);
	BRepBuilderAPI_MakeEdge line2(q1,q2);
	BRepBuilderAPI_MakeEdge line3(q1,p1);
	BRepBuilderAPI_MakeEdge line4(q2,p2);


	TColgp_Array1OfPnt qCurvePoles(1,5);
	qCurvePoles(1) = q2; 
	qCurvePoles(2) = q4; 
	qCurvePoles(3) = r1;
	qCurvePoles(4) = p4; 
	qCurvePoles(5) = p2; 


	Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
	TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 
	BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2);
	BRepBuilderAPI_MakeWire wholeWire ( trapeziumWire,curveEdge);
	BRepTools::Write(wholeWire,"C:/Users/Dell/Desktop/Shapes/xAxis.brep");

	double up=height*4;
	double bottom=height;
	double middle=(up-bottom)/2.0;
	double surfaceArea =getSurfaceArea(wholeWire);
	double left= expectedArea-(expectedArea*variationTolerance);
	double right=expectedArea+(expectedArea*variationTolerance);
	double Tolerance =  0.01*expectedArea; 
	double topMiddle=up;
	double bottomMiddle=bottom; 

	if(fabs(expectedArea-surfaceArea)<Tolerance){
		return wholeWire;

	}else{

		while(fabs(expectedArea-surfaceArea)>Tolerance){

			if(surfaceArea<expectedArea){

				up=up;
				bottom=middle;
				middle=(up+middle)/2;

			}else{

				bottom=bottom;
				up=middle;
				middle=(middle+bottom)/2;

			}
			wholeWire=createNewShapeAccordingToR1height(height,width,ang1Tan,ang2Tan,middle);
			surfaceArea=getSurfaceArea(wholeWire);

		}

	}

	return wholeWire;

}

TopoDS_Wire CVoluteDialog::createNewShapeAccordingToR1height(double height,double width,double ang1Tan,double ang2Tan,double r1Height){



	gp_Pnt r1(width/2.,r1Height,0);
	gp_Pnt q1(0,0,0); 
	gp_Pnt q2(-(height/ang2Tan),height,0);
	gp_Pnt p1(width,0,0);
	gp_Pnt p2(((height/ang1Tan)+width),height,0);
	gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
	gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
	gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);
	gp_Pnt p3(width+((height+height/10)/ang1Tan),(height+height/10),0);
	gp_Pnt p4(width+((height+height/4)/ang1Tan),(height+height/4),0);
	gp_Pnt p5((width+height/ang1Tan),(height+height/3),0);

	BRepBuilderAPI_MakeEdge line1(p1,p2);
	BRepBuilderAPI_MakeEdge line2(q1,q2);
	BRepBuilderAPI_MakeEdge line3(q1,p1);
	BRepBuilderAPI_MakeEdge line4(q2,p2);


	TColgp_Array1OfPnt qCurvePoles(1,5); 
	qCurvePoles(1) = q2; 
	qCurvePoles(2) = q4; 
	qCurvePoles(3) = r1;
	qCurvePoles(4) = p4; 
	qCurvePoles(5) = p2; 

	Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
	TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 

	TopoDS_Wire wholeWire ;
	wholeWire=BRepBuilderAPI_MakeWire (line1,line3,line2,curveEdge);

	return wholeWire;


}


double CVoluteDialog::getSurfaceArea(TopoDS_Wire wireShape){

	BRepBuilderAPI_MakeFace wireFace(wireShape);
	TopoDS_Face face =wireFace;
	BRepGProp gprop; 
	GProp_GProps surfaceProps;

	gprop.SurfaceProperties(wireFace,surfaceProps);

	double surfaceArea=surfaceProps.Mass();
	CString fivs;
	fivs.Format(_T("surface %g\n",surfaceArea));

	return surfaceArea;

}

double  CVoluteDialog::getTrapezuimHeight(double area,double width, double ang1Tan, double ang2Tan){


	double c=area*2;
	double a=1/ang1Tan+1/ang2Tan;
	double b= 2*width;
	double determinaterPart =(b*b) + (4*a*c); 
	double diterminator= sqrt(determinaterPart);
	double root1=((-b)+(diterminator))/(2*a);

	return root1;


}



