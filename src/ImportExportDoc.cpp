// ImportExportDoc.cpp : implementation of the CImportExportDoc class
//


#include "stdafx.h"
#include "ImportExportApp.h"

#include "ImportExportDoc.h"

#include <ImportExport/ImportExport.h>

#include <AISDialogs.h>
#include "res/resource.h"

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
//#include "BOPAlgo_Algo.hxx"

//#include <Standard_Transient.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <BRepLib.hxx>
#define PI 3.14159265


#include <BRepBuilderAPI_Sewing.hxx>



#ifdef _DEBUG
//#define new DEBUG_NEW  // by cascade
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CImportExportDoc

IMPLEMENT_DYNCREATE(CImportExportDoc, OCC_3dDoc)

	BEGIN_MESSAGE_MAP(CImportExportDoc, OCC_3dDoc)
		//{{AFX_MSG_MAP(CImportExportDoc)
		ON_COMMAND(ID_FILE_IMPORT_BREP, OnFileImportBrep)
		ON_COMMAND(ID_FILE_IMPORT_IGES, OnFileImportIges)
		ON_COMMAND(ID_FILE_EXPORT_IGES, OnFileExportIges)
		ON_COMMAND(ID_FILE_IMPORT_STEP, OnFileImportStep)
		ON_COMMAND(ID_FILE_EXPORT_STEP, OnFileExportStep)
		ON_COMMAND(ID_FILE_EXPORT_VRML, OnFileExportVrml)
		ON_COMMAND(ID_FILE_EXPORT_STL, OnFileExportStl)
		ON_COMMAND(ID_BOX, OnBox)
		ON_COMMAND(ID_Cylinder, OnCylinder)
		ON_COMMAND(ID_OBJECT_REMOVE, OnObjectRemove)
		ON_COMMAND(ID_OBJECT_ERASE, OnObjectErase)
		ON_COMMAND(ID_OBJECT_DISPLAYALL, OnObjectDisplayall)
		ON_COMMAND(ID_BUTTONBox1,OnBox1)
		ON_COMMAND(ID_BUTTONbottle,OnBottle)
		ON_COMMAND(ID_BUTTONTest,OnTest)
		ON_COMMAND(ID_BUTTONImportIGESNew,OnIGESFile)
		ON_COMMAND(ID_BUTTONImportBREPNew,OnBREPFile)
		ON_COMMAND(ID_BUTTONImportSTEPNew,OnSTEPFile)
		ON_COMMAND(ID_BUTTONFillet,OnFillet)
		ON_COMMAND(ID_BUTTONFilletDialog,OnFilletDialog)
		ON_COMMAND(ID_BUTTONCut,OnCut)
		ON_COMMAND(ID_BUTTON_BoxMake,OnMakeBoxDrill)
		ON_COMMAND(ID_BUTTON_Volute,OnVolute)
		ON_COMMAND(ID_BUTTON_Thickness,OnThickness)
		ON_COMMAND(ID_BUTTON_FaceOffSet,OnOffSet)
		ON_COMMAND(ID_BUTTON_ImportVolute,OnImportVolute)
		ON_COMMAND(ID_BUTTON_Dual_Volute,OnDualVolute)
		ON_COMMAND(ID_BUTTON_Bearing_Volute,OnBearingVolute)

		//ON_COMMAND(ID_BUTTONImportBREPNew,OnBREPFile)

		//}}AFX_MSG_MAP

	END_MESSAGE_MAP()

	/////////////////////////////////////////////////////////////////////////////
	// CImportExportDoc construction/destruction

	CImportExportDoc::CImportExportDoc()
		: OCC_3dDoc (false)
	{
		/*
		// TRIHEDRON
		Handle(AIS_Trihedron) aTrihedron;
		Handle(Geom_Axis2Placement) aTrihedronAxis=new Geom_Axis2Placement(gp::XOY());
		aTrihedron=new AIS_Trihedron(aTrihedronAxis);
		myAISContext->Display(aTrihedron);
		*/

		m_pcoloredshapeList = new CColoredShapes();



		//boxPointer= &box;
	}
	BRepPrimAPI_MakeBox box(200.,150.,100.);
	BRepPrimAPI_MakeBox* filletBox;
	//filletBox = &box;



	CImportExportDoc::~CImportExportDoc()
	{
		if( m_pcoloredshapeList ) delete m_pcoloredshapeList;
	}


	/////////////////////////////////////////////////////////////////////////////
	// CSerializeDoc serialization

	void CImportExportDoc::Serialize(CArchive& ar)
	{
		if (ar.IsStoring())
		{
			// Put the curent CColoredShape in the archive
			ar << m_pcoloredshapeList;
		}
		else
		{
			// Read from the archive the current CColoredShape
			ar >> m_pcoloredshapeList;

			// Display the new object
			m_pcoloredshapeList->Display(myAISContext);
		}
	}


	/*
	void CImportExportDoc::OnWindowNew3d() 
	{
	((CImportExportApp*)AfxGetApp())->CreateView3D(this);	
	}
	*/

	//  nCmdShow could be :    ( default is SW_RESTORE ) 
	// SW_HIDE   SW_SHOWNORMAL   SW_NORMAL   
	// SW_SHOWMINIMIZED     SW_SHOWMAXIMIZED    
	// SW_MAXIMIZE          SW_SHOWNOACTIVATE   
	// SW_SHOW              SW_MINIMIZE         
	// SW_SHOWMINNOACTIVE   SW_SHOWNA           
	// SW_RESTORE           SW_SHOWDEFAULT      
	// SW_MAX    

	// use pViewClass = RUNTIME_CLASS( CImportExportView3D ) for 3D Views

	void CImportExportDoc::ActivateFrame(CRuntimeClass* pViewClass,int nCmdShow)
	{
		POSITION position = GetFirstViewPosition();
		while (position != (POSITION)NULL)
		{
			CView* pCurrentView = (CView*)GetNextView(position);
			if(pCurrentView->IsKindOf(pViewClass) )
			{
				ASSERT_VALID(pCurrentView);
				CFrameWnd* pParentFrm = pCurrentView->GetParentFrame();
				ASSERT(pParentFrm != (CFrameWnd *)NULL);
				// simply make the frame window visible
				pParentFrm->ActivateFrame(nCmdShow);
			}
		}

	}

	/////////////////////////////////////////////////////////////////////////////
	// CImportExportDoc diagnostics

#ifdef _DEBUG
	void CImportExportDoc::AssertValid() const
	{
		CDocument::AssertValid();
	}

	void CImportExportDoc::Dump(CDumpContext& dc) const
	{
		CDocument::Dump(dc);
	}
#endif //_DEBUG

	/////////////////////////////////////////////////////////////////////////////
	// CImportExportDoc commands







	void CImportExportDoc::OnFileImportBrep()
	{
		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadBREP();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{

			m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, aSeqOfShape->Value(i));
			m_pcoloredshapeList->Display(myAISContext);
		}
		Fit();
	}


	TopoDS_Wire CImportExportDoc::getDualVoluteCrossSection(double area)
	{

		

		double width;
		double bearingFlankHeight;
		double bearingFlankHeightGap;
		double exhaustFlankHeight;
		double bearingSideAngle;
		double exhaustSideAngle;
		double wholeVoluteArea;
		double wholeVoluteTrapziumHeight;
		double trapeziumWidth;
		double trapeziumAngle;
		double tipRadius;
		double dividerWallHeight;
		double areaDivideRatio;
		double dividerAngle;
		double exhaustThickness;
		double bearingThickness;


		width=10*2;
		bearingFlankHeightGap=1.763*2;
		bearingFlankHeight=8;
		bearingSideAngle=41;
		exhaustSideAngle=30;
		wholeVoluteArea=area;
		areaDivideRatio=0.5;
		tipRadius=0.5;
		dividerWallHeight=15;
		exhaustFlankHeight=4;
		dividerAngle=10;
		exhaustThickness=15;
		bearingThickness=15;

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


		gp_Pnt p0(0,-exhaustFlankHeight,0);
		gp_Pnt q0(-width,-exhaustFlankHeight,0);
		gp_Vec leftHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(-1,0,0));
		gp_Vec rightHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));

		gp_Pnt p1(0,0,0);
		gp_Pnt p2(p2x,p2y,0);

		gp_Pnt q1(-width,bearingFlankHeightGap,0);
		gp_Pnt q2(-(width+q2x),bearingFlankHeightGap+q2y,0);
		gp_Pnt r1(-(exhaustThickness*sinA1),(exhaustThickness*cosA1),0);

		gp_Ax1 q1Axis(q1,gp_Dir(0,0,1));
		gp_Vec q1q2Vec=leftHorizontalVec.Rotated(q1Axis,(360-A2)*PI/180);
		q1q2Vec.Normalize();
		q1q2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A2-alfa)*PI/180)));
		gp_Pnt q1End=q1.Translated(q1q2Vec);
		q2=q1End;

		gp_Ax1 p1Axis(p1,gp_Dir(0,0,1));
		gp_Vec p1p2Vec=rightHorizontalVec.Rotated(q1Axis,A1*PI/180);
		p1p2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A1)*PI/180)));
		gp_Pnt p1End=p1.Translated(p1p2Vec);
		TopoDS_Edge q1q2Strait=BRepBuilderAPI_MakeEdge(p1End,p1);
		
		p2=p1End;

		double extraPntX=-(width+10);
		double extraPntY=(width+10)*tan(alfa*PI/180);
		gp_Pnt extraPnt(extraPntX,extraPntY,0);
		TopoDS_Edge extraLine=BRepBuilderAPI_MakeEdge(extraPnt,q1);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,extraLine);

		double angleCheck=atan2(1,1)*180/PI;/////////////////////////////////////
		CString angle;
		angle.Format(_T("angle %g \n"),alfa);

		TopoDS_Edge p1p2Edge=BRepBuilderAPI_MakeEdge(p1,p2);
		TopoDS_Edge p1q1Edge=BRepBuilderAPI_MakeEdge(p1,q1);
		TopoDS_Edge q1q2Edge=BRepBuilderAPI_MakeEdge(q1,q2);
		TopoDS_Edge p2q2Edge=BRepBuilderAPI_MakeEdge(p2,q2);
		TopoDS_Edge p1r1Edge=BRepBuilderAPI_MakeEdge(p1,r1);
		TopoDS_Edge q1r1Edge=BRepBuilderAPI_MakeEdge(q1,r1);
		TopoDS_Edge p0p1Edge=BRepBuilderAPI_MakeEdge(p0,p1);
		TopoDS_Edge q0q1Edge=BRepBuilderAPI_MakeEdge(q0,q1);

		TopoDS_Edge hirizontalLine=BRepBuilderAPI_MakeEdge(gp_Pnt(-10,0,0),gp_Pnt(10,0,0));

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
		gp_Pnt dividePoint;

		Handle_Geom_Curve p1q1LineCurve=BRep_Tool::Curve(p1q1Edge,U1,U2);
		p1q1LineCurve->D1(U2*areaDivideRatio,dividePoint,p1q1Vec);
		p1q1Vec.Multiply(10);
		double magnitude=p1q1Vec.Magnitude();

		CString p1q1;
		p1q1.Format(_T("p1q1 %g \n"),U2);
		double p1q1Length=U2;

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;
		gp_Pnt surfacePnt;

		BRepFill_Filling filledFace;
		filledFace.Add(p1p2Edge,GeomAbs_C0);
		filledFace.Add(p1q1Edge,GeomAbs_C0);
		filledFace.Add(q1q2Edge,GeomAbs_C0);
		filledFace.Add(p2q2Edge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,surfacePnt,uVectorSurface,vVectorSurface);


		//bottom points of splitter
		Handle_Geom_Curve q1q2LineCurve=BRep_Tool::Curve(q1q2Edge,U1,U2);
		gp_Pnt q1q2pointZero;
		gp_Vec q1q2Vec1;
		gp_Vec q1q2Vec2;
		gp_Vec oZ(gp_Dir(0,0,1));
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


		CString q1q2;
		q1q2.Format(_T("q1q2 %g \n"),U2);


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

		CString p1p2;
		p1p2.Format(_T("p1p2 %g \n"),U2);


		gp_Pnt p1q1pointZero;
		gp_Vec p1q1Vec1;
		gp_Vec p1q1Vec2;

		p1q1LineCurve->D2(U1,p1q1pointZero,p1q1Vec1,p1q1Vec2);

		p1q1pointZero.Translate(p1q1Vec1);

		TopoDS_Vertex vertp1q1 =BRepBuilderAPI_MakeVertex(p1q1pointZero);


		//gp_Pnt intersectionPointOfBottomEdgeLines=getMinimumDistancePoint(bearingSideBottomEdge,exhaustSideBottomEdge);

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

		TopoDS_Vertex LeftbottomVert=BRepBuilderAPI_MakeVertex(leftBottomPiontOfDividerWall);
		TopoDS_Vertex RightbottpmVert=BRepBuilderAPI_MakeVertex(rightBottomPointOfDividerWall);
		TopoDS_Vertex dividerOnBaseVert=BRepBuilderAPI_MakeVertex(dividerPointOnBase);

		gp_Pnt minimumDistancePointOnp1p2=getMinimumDistancePoint(p1p2Edge,RightbottpmVert);
		gp_Pnt minimumDistancePointOnq1q2=getMinimumDistancePoint(q1q2Edge,LeftbottomVert);

		TopoDS_Vertex mdVert=BRepBuilderAPI_MakeVertex(minimumDistancePointOnq1q2);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,mdVert);
		TopoDS_Edge mdline=BRepBuilderAPI_MakeEdge(minimumDistancePointOnq1q2,leftBottomPiontOfDividerWall);
		TopoDS_Edge mdline1=BRepBuilderAPI_MakeEdge(minimumDistancePointOnp1p2,rightBottomPointOfDividerWall);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY, mdline1);

		double bearingSideOuterWallScalar=minimumDistancePointOnq1q2.Distance(q2);
		double exhaustSideOuterWallScalar=minimumDistancePointOnp1p2.Distance(p2);


		gp_Ax1 rightAx(rightBottomPointOfDividerWall,gp_Dir(0,0,1));
		gp_Vec divideLineVec(dividerPointOnBase,bottomParallelEdgesIntersectionPnt);
		divideLineVec.Normalize();
		gp_Vec rightRotatedVec=divideLineVec.Rotated(rightAx,((360-(dividerAngle/2)))*PI/180);
		rightRotatedVec.Multiply(exhaustSideOuterWallScalar);
		gp_Pnt rightTopPointOfDividerWall=rightBottomPointOfDividerWall.Translated(rightRotatedVec);
		TopoDS_Edge rightDividerWall=BRepBuilderAPI_MakeEdge(rightTopPointOfDividerWall,rightBottomPointOfDividerWall);

		gp_Ax1 leftAx(rightBottomPointOfDividerWall,gp_Dir(0,0,1));
		gp_Vec leftRotateVec=divideLineVec.Rotated(leftAx,(dividerAngle/2)*PI/180); 
		leftRotateVec.Multiply(bearingSideOuterWallScalar);
		gp_Pnt leftTopPointofDividerWall=leftBottomPiontOfDividerWall.Translated(leftRotateVec);
		TopoDS_Edge leftDividerWall=BRepBuilderAPI_MakeEdge(leftTopPointofDividerWall,leftBottomPiontOfDividerWall);

		TopoDS_Edge topLeftBaseEdge=BRepBuilderAPI_MakeEdge(q1q2EndPoint,leftTopPointofDividerWall);
		Handle_Geom_Curve topLeftBaseLineCurve=BRep_Tool::Curve(topLeftBaseEdge,U1,U2);
		gp_Pnt leftCircleCentrePoint;
		topLeftBaseLineCurve->D0(U2/2,leftCircleCentrePoint);
		double topLeftCircleRadius=U2/2;
		gp_Ax2 leftCircleAxis(leftCircleCentrePoint,gp_Dir(0,0,1));
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
		gp_Ax2 rightCircleAxis(rightCircleCentrePoint,gp_Dir(0,0,1));
		gp_Circ topRightCircle(rightCircleAxis,topRightCircleRaduis);
		GC_MakeArcOfCircle topRightHalfOfCircle(topRightCircle,p1p2EndPoint,rightTopPointOfDividerWall,false);
		Handle(Geom_TrimmedCurve) rightCircleArc=topRightHalfOfCircle;
		TopoDS_Edge rightCircle=BRepBuilderAPI_MakeEdge(rightCircleArc);
		TopoDS_Edge rightCircleArcBezierCurveEdge =convertTrimmToBezier(rightCircleArc,rightRotatedVec,p1p2Vec1,exhaustSideOuterWallScalar);
		rightCircle = rightCircleArcBezierCurveEdge;


		gp_Ax2 bottomCircleAxis(bottomParallelEdgesIntersectionPnt,gp_Dir(0,0,1));
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
		return completeCrossSectionWire;
	}

	void CImportExportDoc::displayPoint(gp_Pnt pnt)
	{
		TopoDS_Vertex vert =BRepBuilderAPI_MakeVertex(pnt);
		m_pcoloredshapeList->Add(Quantity_NOC_MATRABLUE,vert);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();


	}

	void CImportExportDoc::displayPoint( Quantity_NameOfColor color,gp_Pnt pnt)
	{
		TopoDS_Vertex vert =BRepBuilderAPI_MakeVertex(pnt);
		m_pcoloredshapeList->Add(color,vert);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();


	}


	TopoDS_Edge CImportExportDoc::convertTrimmToBezier(Handle_Geom_Curve curve,gp_Vec leftVec,gp_Vec rightVec,double scaler)
	{
		GeomConvert convert;
		Handle_Geom_BSplineCurve ArcBSplineCurve=convert.CurveToBSplineCurve(curve);
		GeomConvert_BSplineCurveToBezierCurve bsplineToBezierCurve(ArcBSplineCurve);
		Handle_Geom_BezierCurve Arc1BezierCurve=bsplineToBezierCurve.Arc(1);
		Handle_Geom_BezierCurve Arc2BezierCurve=bsplineToBezierCurve.Arc(2);

		leftVec.Normalize();
		rightVec.Normalize();

		leftVec.Multiply(scaler/10);
		rightVec.Multiply(scaler/10);
		
		TopoDS_Edge curveEdge=BRepBuilderAPI_MakeEdge( ArcBSplineCurve);
		TopoDS_Edge e1=BRepBuilderAPI_MakeEdge(Arc1BezierCurve);
		double nbPoles=Arc1BezierCurve->NbPoles();

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

		BRepLib_FuseEdges fuse(BRepAlgoAPI_Fuse(e1, e2).Shape());
		TopoDS_Shape e3=fuse.Shape();
		
		TopoDS_Edge newEdge;
		TopExp_Explorer EdgeExplorer(e3,TopAbs_EDGE);
		int i=0;
		while(EdgeExplorer.More()){
			TopoDS_Edge aEdge=TopoDS::Edge(EdgeExplorer.Current());
			newEdge=aEdge;
			EdgeExplorer.Next();

		}

		double lastParaOfArc1=Arc1BezierCurve->LastParameter();
		double lastParaOfArc2=Arc2BezierCurve->LastParameter();
		double Arc1paraRatio=lastParaOfArc1/10;
		double Arc2paraRatio=lastParaOfArc2/10;

		gp_Pnt arc1PointArray[10];
		gp_Pnt arc2PointArray[10];
		TColgp_Array1OfPnt bsplinePoles(1,21);
		bsplinePoles(1)=Arc1BezierCurve->StartPoint();
		for(int i=0;i<10;i++){
			Arc1BezierCurve->D0(Arc1paraRatio*(i+1),arc1PointArray[i]);
			bsplinePoles(i+2)=arc1PointArray[i];
		}

		for(int i=0;i<10;i++){
			Arc2BezierCurve->D0(Arc2paraRatio*(i+1),arc2PointArray[i]);
			bsplinePoles(i+12)=arc2PointArray[i];
		}

		GeomAPI_PointsToBSpline bspline(bsplinePoles); 
		TopoDS_Edge bsplineEdge = BRepBuilderAPI_MakeEdge(bspline.Curve());
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,bsplineEdge);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();
		return bsplineEdge;
	}

	void CImportExportDoc::OnBearingVolute(){

		//		 q2            p1
		//		   \          /
		//			\        /
		//			 \______/
		//			q1     p1


		double width;
		double bearingFlankHeightGap;
		double exhaustFlankHeight;
		double bearingSideAngle;
		double exhaustSideAngle;
		double wholeVoluteArea;
		double wholeVoluteTrapziumHeight;
		double trapeziumWidth;
		double trapeziumAngle;
		double tipRadius;
		double dividerWallHeight;
		double areaDivideRatio;
		double dividerAngle;
		double exhaustThickness;
		double bearingThickness;
		double transitionPartLength;
		double exitPipeLength;
		double exitDividerAngle;
		double voluteRadius;
		double exitDividerWallWidth=6;
		double exitPipeRadius=width*3;

		width=10*2;
		bearingFlankHeightGap=1.763*2;
		bearingSideAngle=41;
		exhaustSideAngle=30;
		wholeVoluteArea=2000;
		areaDivideRatio=0.5;
		tipRadius=0.5;
		dividerWallHeight=15;
		exhaustFlankHeight=10;
		dividerAngle=10;
		exhaustThickness=15;
		bearingThickness=15;
		transitionPartLength=-150;
		exitPipeLength=-60;
		exitDividerAngle=90;
		exitDividerWallWidth=6;
		exitPipeRadius=width*3;

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

		gp_Vec leftHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(-1,0,0));
		gp_Vec rightHorizontalVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));

		gp_Pnt p1(0,0,0);
		gp_Pnt p2(p2x,p2y,0);

		gp_Pnt q1(-width,bearingFlankHeightGap,0);
		gp_Pnt q2(-(width+q2x),bearingFlankHeightGap+q2y,0);
		gp_Pnt r1(-(exhaustThickness*sinA1),(exhaustThickness*cosA1),0);

		gp_Ax1 q1Axis(q1,gp_Dir(0,0,1));
		gp_Vec q1q2Vec=leftHorizontalVec.Rotated(q1Axis,(360-A2)*PI/180);
		q1q2Vec.Normalize();
		q1q2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A2-alfa)*PI/180)));
		gp_Pnt q1End=q1.Translated(q1q2Vec);
		q2=q1End;

		gp_Ax1 p1Axis(p1,gp_Dir(0,0,1));
		gp_Vec p1p2Vec=rightHorizontalVec.Rotated(q1Axis,A1*PI/180);
		p1p2Vec.Multiply(initialLengthOfOuterLine+(h/sin((A1)*PI/180)));
		gp_Pnt p1End=p1.Translated(p1p2Vec);
		TopoDS_Edge q1q2Strait=BRepBuilderAPI_MakeEdge(p1End,p1);
		
		p2=p1End;
		
		double extraPntX=-(width+10);
		double extraPntY=(width+10)*tan(alfa*PI/180);
		gp_Pnt extraPnt(extraPntX,extraPntY,0);
		TopoDS_Edge extraLine=BRepBuilderAPI_MakeEdge(extraPnt,q1);
		
		double angleCheck=atan2(1,1)*180/PI;
		CString angle;
		angle.Format(_T("angle %g \n"),alfa);

		TopoDS_Edge p1p2Edge=BRepBuilderAPI_MakeEdge(p1,p2);
		TopoDS_Edge p1q1Edge=BRepBuilderAPI_MakeEdge(p1,q1);
		TopoDS_Edge q1q2Edge=BRepBuilderAPI_MakeEdge(q1,q2);
		TopoDS_Edge p2q2Edge=BRepBuilderAPI_MakeEdge(p2,q2);
		TopoDS_Edge p1r1Edge=BRepBuilderAPI_MakeEdge(p1,r1);
		TopoDS_Edge q1r1Edge=BRepBuilderAPI_MakeEdge(q1,r1);

		TopoDS_Edge hirizontalLine=BRepBuilderAPI_MakeEdge(gp_Pnt(-10,0,0),gp_Pnt(10,0,0));

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
		gp_Pnt dividePoint;

		Handle_Geom_Curve p1q1LineCurve=BRep_Tool::Curve(p1q1Edge,U1,U2);
		p1q1LineCurve->D1(U2*areaDivideRatio,dividePoint,p1q1Vec);
		p1q1Vec.Multiply(10);
		double magnitude=p1q1Vec.Magnitude();

		CString p1q1;
		p1q1.Format(_T("p1q1 %g \n"),U2);
		double p1q1Length=U2;

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;
		gp_Pnt surfacePnt;

		BRepFill_Filling filledFace;
		filledFace.Add(p1p2Edge,GeomAbs_C0);
		filledFace.Add(p1q1Edge,GeomAbs_C0);
		filledFace.Add(q1q2Edge,GeomAbs_C0);
		filledFace.Add(p2q2Edge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,surfacePnt,uVectorSurface,vVectorSurface);

		double uParaOfDividePoint=U2*areaDivideRatio;
		double bottomScaler=(dividerWallHeight-exhaustFlankHeight)/cosAlfa-uParaOfDividePoint*tanAlfa+tipRadius;
		
		//gp_Vec dividerWallBottomVec=vVectorSurface.Multiplied(bottomScaler);
		//gp_Pnt dividerBottomPnt=dividePoint.Translated(dividerWallBottomVec);
		//gp_Ax1 rightAx(dividerBottomPnt,gp_Dir(0,0,1));
		//gp_Vec leftRotatedVector=dividerWallBottomVec.Rotated(rightAx,dividerAngle*PI/180);
		//gp_Pnt leftDividerWallBottomPnt=dividePoint.Translated(leftRotatedVector);
		//gp_Vec rightRotatedVector=dividerWallBottomVec.Rotated(rightAx,(360-dividerAngle)*PI/180);
		//gp_Pnt rightDividerWallBottomPnt=dividePoint.Translated(rightRotatedVector);
		//TopoDS_Vertex vert=BRepBuilderAPI_MakeVertex(leftDividerWallBottomPnt);
		//TopoDS_Vertex vert1=BRepBuilderAPI_MakeVertex(rightDividerWallBottomPnt);
		//TopoDS_Edge divideEdge=BRepBuilderAPI_MakeEdge(dividePoint,dividerBottomPnt);
		//TopoDS_Vertex middle=BRepBuilderAPI_MakeVertex(r1);
		//bottom points of splitter

		Handle_Geom_Curve q1q2LineCurve=BRep_Tool::Curve(q1q2Edge,U1,U2);
		gp_Pnt q1q2pointZero;
		gp_Vec q1q2Vec1;
		gp_Vec q1q2Vec2;
		gp_Vec oZ(gp_Dir(0,0,1));
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


		CString q1q2;
		q1q2.Format(_T("q1q2 %g \n"),U2);


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

		CString p1p2;
		p1p2.Format(_T("p1p2 %g \n"),U2);


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

		gp_Ax1 rightAx(rightBottomPointOfDividerWall,gp_Dir(0,0,1));
		gp_Vec divideLineVec(dividerPointOnBase,bottomParallelEdgesIntersectionPnt);
		divideLineVec.Normalize();
		gp_Vec rightRotatedVec=divideLineVec.Rotated(rightAx,((360-(dividerAngle/2)))*PI/180);
		rightRotatedVec.Multiply(exhaustSideOuterWallScalar);
		gp_Pnt rightTopPointOfDividerWall=rightBottomPointOfDividerWall.Translated(rightRotatedVec);
		TopoDS_Edge rightDividerWall=BRepBuilderAPI_MakeEdge(rightTopPointOfDividerWall,rightBottomPointOfDividerWall);

		gp_Ax1 leftAx(rightBottomPointOfDividerWall,gp_Dir(0,0,1));
		gp_Vec leftRotateVec=divideLineVec.Rotated(leftAx,(dividerAngle/2)*PI/180); 
		leftRotateVec.Multiply(bearingSideOuterWallScalar);
		gp_Pnt leftTopPointofDividerWall=leftBottomPiontOfDividerWall.Translated(leftRotateVec);
		TopoDS_Edge leftDividerWall=BRepBuilderAPI_MakeEdge(leftTopPointofDividerWall,leftBottomPiontOfDividerWall);
		

		TopoDS_Edge topLeftBaseEdge=BRepBuilderAPI_MakeEdge(q1q2EndPoint,leftTopPointofDividerWall);
		Handle_Geom_Curve topLeftBaseLineCurve=BRep_Tool::Curve(topLeftBaseEdge,U1,U2);
		gp_Pnt leftCircleCentrePoint;
		topLeftBaseLineCurve->D0(U2/2,leftCircleCentrePoint);
		double topLeftCircleRadius=U2/2;
		gp_Ax2 leftCircleAxis(leftCircleCentrePoint,gp_Dir(0,0,1));
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
		gp_Ax2 rightCircleAxis(rightCircleCentrePoint,gp_Dir(0,0,1));
		gp_Circ topRightCircle(rightCircleAxis,topRightCircleRaduis);
		GC_MakeArcOfCircle topRightHalfOfCircle(topRightCircle,p1p2EndPoint,rightTopPointOfDividerWall,false);
		Handle(Geom_TrimmedCurve) rightCircleArc=topRightHalfOfCircle;
		TopoDS_Edge rightCircle=BRepBuilderAPI_MakeEdge(rightCircleArc);
		TopoDS_Edge rightCircleArcBezierCurveEdge =convertTrimmToBezier(rightCircleArc,rightRotatedVec,p1p2Vec1,exhaustSideOuterWallScalar);
		rightCircle = rightCircleArcBezierCurveEdge;

		gp_Ax2 bottomCircleAxis(bottomParallelEdgesIntersectionPnt,gp_Dir(0,0,1));
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
		TopoDS_Compound compound;
		TopoDS_Builder abuilder;

		/*abuilder.MakeCompound(compound);
		abuilder.Add(compound,p1p2Edge);
		abuilder.Add(compound,q1q2Edge);
		abuilder.Add(compound,p1q1Edge);
		abuilder.Add(compound,leftCircle);
		abuilder.Add(compound,rightCircle);
		abuilder.Add(compound,bottomTipCircle);
		abuilder.Add(compound,leftDividerWall);
		abuilder.Add(compound,rightDividerWall);
		abuilder.Add(compound,mdline);
		abuilder.Add(compound,mdline1);*/


		BRepBuilderAPI_MakeWire completeCrossSection;
		completeCrossSection.Add(p1p2Edge);
		completeCrossSection.Add(p1q1Edge);
		completeCrossSection.Add(q1q2Edge);
		completeCrossSection.Add(leftCircle);
		completeCrossSection.Add(leftDividerWall);
		completeCrossSection.Add(bottomTipCircle);
		completeCrossSection.Add(rightDividerWall);
		completeCrossSection.Add(rightCircle);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,completeCrossSection);

		voluteRadius=width*5;
		BRepOffsetAPI_ThruSections sections;
		gp_Trsf transfer;
		gp_Ax1 rotationAxis(gp_Pnt(0,-voluteRadius,0),gp_Dir(1,0,0));

		int numberOfCrossSections=36;
		double toungArea=wholeVoluteArea*0.01;
		TopoDS_Wire arrayOfCrossSections[36];
		TopoDS_Wire newCrossSection;
		double areaReductionRate=(wholeVoluteArea-toungArea)/numberOfCrossSections;
		double rotationAngle=(360/numberOfCrossSections)*PI/180;


		double expectedArea=wholeVoluteArea-(areaReductionRate*0);
		newCrossSection=getDualVoluteCrossSection(expectedArea);
		transfer.SetRotation(rotationAxis,360*PI/180);
		BRepBuilderAPI_Transform rotated(newCrossSection,transfer);
		arrayOfCrossSections[0]=TopoDS::Wire(rotated.Shape());//
		
		sections.AddWire(arrayOfCrossSections[0]);

		for(int i=1;i<=numberOfCrossSections-1;i++){

			double expectedArea=wholeVoluteArea-(areaReductionRate*i);
			newCrossSection=getDualVoluteCrossSection(expectedArea);
			transfer.SetRotation(rotationAxis,rotationAngle*i);
			BRepBuilderAPI_Transform rotated(newCrossSection,transfer);
			arrayOfCrossSections[i]=TopoDS::Wire(rotated.Shape());
			sections.AddWire(arrayOfCrossSections[i]);
	
		}

		double expectedArea1=wholeVoluteArea-(areaReductionRate*numberOfCrossSections);
		newCrossSection=getDualVoluteCrossSection(expectedArea1); 
		transfer.SetRotation(rotationAxis,360*PI/180);
		BRepBuilderAPI_Transform rotated1(newCrossSection,transfer);
		arrayOfCrossSections[numberOfCrossSections-1]=TopoDS::Wire(rotated1.Shape());
		
		sections.AddWire(arrayOfCrossSections[numberOfCrossSections-1]);
		sections.Build();
		TopoDS_Shape VoluteShape(sections.Shape());
		BRepTools::Write(VoluteShape,"D:/Breps/DualVolute/Volute.brep");
		TopoDS_Wire completeCrossSectionWire=completeCrossSection;

		Handle_Geom_Curve p2q2LineCurve=BRep_Tool::Curve(p2q2Edge,U1,U2);
		gp_Pnt centrePointOfCircle;
		p2q2LineCurve->D0(U2/2,centrePointOfCircle);

		TopoDS_Vertex p2q2Vert=BRepBuilderAPI_MakeVertex(centrePointOfCircle); 

		gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,0,1));
		gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
		TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


		double maxWidthOfDividerWall=leftTopPointofDividerWall.Distance(rightTopPointOfDividerWall);

		TopoDS_Wire leftPartOfExitPipeTransition;
		TopoDS_Wire rightPartOfExitPipeTransition;

		getExitPipeTrainsition(leftPartOfExitPipeTransition,rightPartOfExitPipeTransition,centrePointOfCircle,exitDividerAngle,leftBottomVector,maxWidthOfDividerWall,exitPipeRadius,divideLineVec);

		//Transition
		double tRLength=getMaximumHeight(arrayOfCrossSections[27]);
		transitionPartLength=-voluteRadius-tRLength;
		gp_Vec exitPipeTranslateVec(gp_Dir(0,0,1));
		exitPipeTranslateVec.Multiply(transitionPartLength);
		gp_Trsf translate;
		translate.SetTranslation(exitPipeTranslateVec);

		BRepBuilderAPI_Transform leftTransform( leftPartOfExitPipeTransition,translate,true);
		TopoDS_Wire leftTransition= TopoDS::Wire(leftTransform.ModifiedShape( leftPartOfExitPipeTransition));


		BRepBuilderAPI_Transform rightTransform(rightPartOfExitPipeTransition,translate,true);
		TopoDS_Wire rightTransition=TopoDS::Wire(rightTransform.ModifiedShape(rightPartOfExitPipeTransition));

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

		//ExitPipe
		gp_Vec exitPipeTranslateVec1(gp_Dir(0,0,1));

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

		TopoDS_Shape leftExitPipeEnding=leftExitPipeSections.Shape();
		TopoDS_Shape rightExitPipeEnding=rightExitPipeSection.Shape();

		TopoDS_Compound transitionCom;
		TopoDS_Builder trBuilder;
		trBuilder.MakeCompound(transitionCom);
		trBuilder.Add(transitionCom,leftTransitionPartShape);
		trBuilder.Add(transitionCom,rightTransitionPartShape);
		trBuilder.Add(transitionCom,leftExitPipeEnding);
		trBuilder.Add(transitionCom,rightExitPipeEnding);
		trBuilder.Add(transitionCom,VoluteShape);

		TopoDS_Wire * exitPipetransitionSection;

		exitPipetransitionSection= getExitPipeSection();

		/*m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftTransitionPartShape);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,rightTransitionPartShape);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftVoluteWire);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,rightTopFilletEdge);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,rightBottomFilletEdge);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,rightExitArcEdge);
		m_pcoloredshapeList->Add(Quantity_NOC_RED, rightExitDividerAfterFillet);
		m_pcoloredshapeList->Add(Quantity_NOC_RED, leftTransition);
		m_pcoloredshapeList->Add(Quantity_NOC_RED, rightTransition);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,bottomRightExitVert);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,BottomrightFilletVertex[1]);*/
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,transitionCom);
		//gp_Ax2 rightCircleAxis(rightCircleCentrePoint,gp_Dir(0,0,1));
		//gp_Circ topRightCircle(rightCircleAxis,topRightCircleRaduis);
		//GC_MakeArcOfCircle topRightHalfOfCircle(topRightCircle,p1p2EndPoint,rightTopPointOfDividerWall,false);
		//abuilder.Add(compound,divideEdge);

		BRepTools::Write(leftTransitionPartShape,"D:/Breps/DualVolute/leftTransitionPartShape.brep");
		BRepTools::Write(rightTransitionPartShape,"D:/Breps/DualVolute/rightTransitionPartShape.brep");

		BRepTools::Write(leftExitPipeEnding,"D:/Breps/DualVolute/leftExitPipeEnding.brep");
		BRepTools::Write(rightExitPipeEnding,"D:/Breps/DualVolute/rightExitPipeEnding.brep");

		BRepTools::Write(transitionCom,"D:/Breps/DualVolute/completeShape.brep");
		CString str;



		//str.Format(_T("upara %g \n"),endUpara);

		//AfxMessageBox(str);
		/*m_pcoloredshapeList->Add(Quantity_NOC_GREEN,p2q2Vert);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,leftPartOfExitPipe);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,rightPartOfExitPipe);*/


		//m_pcoloredshapeList->Add(Quantity_NOC_RED,topPartOfLeftDividerEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,exitCircleEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,leftExitDividerEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,rightExitDividerEdge);
		//	m_pcoloredshapeList->Add(Quantity_NOC_IVORY,topPartOfLeftArcEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,rightExitCircleArcEdge);

		/*m_pcoloredshapeList->Add(Quantity_NOC_IVORY,leftTopFilletEdge);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,leftBottomFilletEdge);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,edge1);*/
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,finalShape);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,compound);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,p1p2Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,p1q1Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,q1q2Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,p2q2Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,vp1);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,vp2);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vq1);  
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vq2);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vert);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,vert1);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,divideEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,hirizontalLine);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,q1q2vertZero);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,p1p2vertZero);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,q1r1Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,p1r1Edge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,bearingSideBottomEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,exhaustSideBottomEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vertp1q1);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,bottomVertLeft);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,intersectionVert);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,bottpmVertRight);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,dividerOnBaseVert);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,p1p2ParallelEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,q1q2ParallelEdge);
		/*	m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,rightDividerWall);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,leftDividerWall);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftCircle);
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,rightCircle);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,bottomTipCircle);*/
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,completeCrossSectionWire);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();

	}

	double CImportExportDoc::getMaximumHeight(TopoDS_Wire section){

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


		//gp_Pnt maxDistPntOnEdge2=maxDistExtremaOnEdge2.PointOnE1(1);
		//TopoDS_Vertex vert2=BRepBuilderAPI_MakeVertex(maxDistPntOnEdge2); 


		BRepExtrema_ExtPC maxDistancePoint(vertexArrayOfEdge0[0],edges[2]);
		gp_Pnt maxPoint=maxDistancePoint.Point(1);
		TopoDS_Vertex maxVertex=BRepBuilderAPI_MakeVertex(maxPoint);
		BRepExtrema_DistShapeShape maxDistExtrema(edges[2],maxVertex);
		double value=maxDistExtrema.Value();
		gp_Pnt PointOnEdge2=BRep_Tool::Pnt(vertexArrayOfEdge0[0]);

		double distanceValue=PointOnEdge2.Distance(maxPoint);


		CString str;
		str.Format(_T("unit %g"), distanceValue);
		AfxMessageBox(str);

		/*m_pcoloredshapeList->Add(Quantity_NOC_GREEN,section);
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,maxVertex);
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,vertexArray[0]);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,edges[8]);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,edges[0]);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vertexArrayOfEdge0[0]);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,vertexArrayOfEdge8[0]);*/
		m_pcoloredshapeList->Display(myAISContext);
		Fit();

		return distanceValue;
	} 


	TopoDS_Wire* CImportExportDoc::getExitPipeSection(){

		TopoDS_Wire crossSection[2];

		return crossSection;


	}

	double CImportExportDoc::getDividerWallMaximumWidth(TopoDS_Wire section)
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

		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,edges[3]);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,edges[5]);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();
		return 0;

	}

	void CImportExportDoc::getExitPipeTrainsition(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,double exitDividerWallWidth,double exitPipeRadius,gp_Vec divideLineVec)
	{

		gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,0,1));
		gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
		TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


		gp_Vec straitExitDividerVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
		straitExitDividerVec.Rotate(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),90*PI/180);
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
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,topVert);

		/*Handle_Geom_Curve leftExitArcCurve=convert.CurveToBSplineCurve(leftExitCircleArcTrimmedCurve);
		TopoDS_Edge leftExitArcEdgeBS=BRepBuilderAPI_MakeEdge(leftExitArcCurve);
		leftExitArcCurve=BRep_Tool::Curve(leftExitDividerEdge,U1,U2);
		Handle_Geom_BSplineCurve topPartOfLeftExitArc=Handle_Geom_BSplineCurve::DownCast(leftExitArcCurve->Copy()); 
		Handle_Geom_BSplineCurve bottomPartOfLeftExitArc=Handle_Geom_BSplineCurve::DownCast(leftExitArcCurve->Copy()); 
		topPartOfLeftExitArc->Segment(U2/2,U2);
		TopoDS_Edge topPartOfLeftArcEdgeBS=BRepBuilderAPI_MakeEdge(topPartOfLeftExitArc);*/

		double filletRadius=8;

		gp_Pln filletPlane(centrePointOfCircle,gp_Dir(0,0,1));
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




	void CImportExportDoc::getExitPipeEnding(TopoDS_Wire& w1,TopoDS_Wire& w2,gp_Pnt centrePointOfCircle,double exitDividerAngle,gp_Vec leftBottomVector,double exitDividerWallWidth,double exitPipeRadius,gp_Vec divideLineVec)
	{

		gp_Ax2 exitePipeCircleAxis(centrePointOfCircle,gp_Dir(0,0,1));
		gp_Circ exitPipeCircle(exitePipeCircleAxis,exitPipeRadius);
		TopoDS_Edge exitCircleEdge=BRepBuilderAPI_MakeEdge(exitPipeCircle);


		gp_Vec straitExitDividerVec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
		straitExitDividerVec.Rotate(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),exitDividerAngle*PI/180);
		gp_Pnt zeroPnt(0,0,0);
		zeroPnt.Translate(straitExitDividerVec);

		gp_Dir exitDividerDir(straitExitDividerVec);

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

		gp_Pln filletPlane(centrePointOfCircle,gp_Dir(0,0,1));
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



	gp_Pnt CImportExportDoc::getMaxDistancePointBetweenVertexAndEdge(TopoDS_Vertex vertex,TopoDS_Edge edge){
		BRepExtrema_ExtPC maxDistancePoint(vertex,edge);
		gp_Pnt maxPoint=maxDistancePoint.Point(1);
		return maxPoint;

	}

	gp_Pnt CImportExportDoc::getMinimumDistancePoint(TopoDS_Edge edge,TopoDS_Vertex vertex)
	{

		BRepExtrema_DistShapeShape minimumDist(edge,vertex,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);
		gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

		return divideIntersectionPointOnCurve;
	}

	gp_Pnt CImportExportDoc::getMinimumDistancePoint(TopoDS_Edge edge1,TopoDS_Edge edge2)
	{

		TopoDS_Wire wire1=BRepBuilderAPI_MakeWire(edge1);
		TopoDS_Wire wire2=BRepBuilderAPI_MakeWire(edge2);


		BRepExtrema_DistShapeShape minimumDist(wire1,wire2,Extrema_ExtFlag_MIN,Extrema_ExtAlgo_Grad);
		gp_Pnt divideIntersectionPointOnCurve=minimumDist.PointOnShape1(1);

		return divideIntersectionPointOnCurve;
	}


	TopoDS_Wire CImportExportDoc::createSplitterFromThreeEdges(double dividePointPara,double topWidth,double bottomWidth,double heightRatio,TopoDS_Wire importedWire )
	{

		double leftStartParaOfSplitterOnBase;
		double rightStartParaOfSplitterOnBase;

		double bottomWidthLeftParaOnBase;
		double bottomWidthRightParaOnBase;

		double tipRadius=15;
		double topRadius=3;
		leftStartParaOfSplitterOnBase=dividePointPara-(topWidth/20);
		rightStartParaOfSplitterOnBase=dividePointPara+(topWidth/20);

		bottomWidthLeftParaOnBase=dividePointPara-(bottomWidth/20);
		bottomWidthRightParaOnBase=dividePointPara+(bottomWidth/20);

		TopoDS_Edge edges[4];
		TopExp_Explorer anEdgeExplorer(importedWire, TopAbs_EDGE);
		int i=0;
		while(anEdgeExplorer.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			edges[i]=anEdge;
			anEdgeExplorer.Next();
			i++;

		}

		TopoDS_Edge importedEdge=edges[0];
		Standard_Real U1=0;
		Standard_Real U2=1;
		Standard_Real uParaOfIntersectionPoint;

		gp_Pnt pnt;
		gp_Vec V1;
		gp_Vec V2;
		gp_Vec V3;
		Handle_Geom_Curve curve= BRep_Tool::Curve(importedEdge,U1, U2);


		Standard_Real leftTurningPntPara=U2*(44/U2);
		Standard_Real rightTurningPntPara=U2*(44/U2)*3;

		Handle(Geom_BSplineCurve) straitLeft = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) middleCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) straitRight = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());

		straitLeft->Segment(U1,leftTurningPntPara);
		middleCurve->Segment(leftTurningPntPara,rightTurningPntPara);
		straitRight->Segment(rightTurningPntPara,U2);

		TopoDS_Edge straitLeftBaseEdge=BRepBuilderAPI_MakeEdge(straitLeft);
		TopoDS_Edge straitRightBaseEdge=BRepBuilderAPI_MakeEdge(straitRight);
		TopoDS_Edge middleCurveEdge=BRepBuilderAPI_MakeEdge(middleCurve);

		TopoDS_Vertex vertex[2];
		TopExp_Explorer vertexExplorer(importedWire,TopAbs_VERTEX);
		int j=0;
		while(vertexExplorer.More()){
			TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
			vertex[j]=aVertex;
			vertexExplorer.Next();
			j++;

		}

		Standard_Real basePara1;
		Standard_Real basePara2;
		gp_Pnt dividePnt;

		gp_Pnt leftStartOfSplitterOnBasePnt;
		gp_Pnt rightStartOfSplitterOnBasePnt;
		gp_Pnt bottomWidthLeftPntOnBase;
		gp_Pnt bottomWidthRightPntOnBase;


		gp_Vec devidePointVec;

		TopoDS_Edge horizontalBaseEdge=BRepBuilderAPI_MakeEdge(vertex[0],vertex[1]);
		Handle_Geom_Curve horizontalBaseLine=BRep_Tool::Curve(horizontalBaseEdge,basePara1,basePara2);

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;

		BRepFill_Filling filledFace;
		filledFace.Add(edges[0],GeomAbs_C0);
		filledFace.Add(horizontalBaseEdge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,dividePnt,uVectorSurface,vVectorSurface);

		horizontalBaseLine->D1(dividePointPara,dividePnt,devidePointVec);
		horizontalBaseLine->D0(leftStartParaOfSplitterOnBase,leftStartOfSplitterOnBasePnt);
		horizontalBaseLine->D0(rightStartParaOfSplitterOnBase,rightStartOfSplitterOnBasePnt);
		horizontalBaseLine->D0(bottomWidthLeftParaOnBase,bottomWidthLeftPntOnBase);
		horizontalBaseLine->D0(bottomWidthRightParaOnBase,bottomWidthRightPntOnBase);


		BRepBuilderAPI_MakeVertex divideVertex(dividePnt);



		gp_Dir dir(uVectorSurface); 
		gp_Lin areaDivideLine(dividePnt,dir);
		gp_Lin leftSplitterStartLine(leftStartOfSplitterOnBasePnt,dir);
		gp_Lin rightSplitterStartLine(rightStartOfSplitterOnBasePnt,dir);
		gp_Lin bottomWidthLeftLine(bottomWidthLeftPntOnBase,dir);
		gp_Lin bottomWidthRightLine(bottomWidthRightPntOnBase,dir);



		TopoDS_Edge areaDivideEdge= BRepBuilderAPI_MakeEdge(areaDivideLine);
		TopoDS_Edge leftSplitterStartEdge=BRepBuilderAPI_MakeEdge(leftSplitterStartLine);
		TopoDS_Edge rightSplitterStartEdge=BRepBuilderAPI_MakeEdge(rightSplitterStartLine);
		TopoDS_Edge bottomWidthLeftEdge=BRepBuilderAPI_MakeEdge(bottomWidthLeftLine);
		TopoDS_Edge bottomWidthRightEdge=BRepBuilderAPI_MakeEdge(bottomWidthRightLine);

		gp_Pnt leftStartOfSplitterPnt = getMinimumDistancePoint(middleCurveEdge,leftSplitterStartEdge);
		gp_Pnt rightStartOfSplitterPnt = getMinimumDistancePoint(middleCurveEdge,rightSplitterStartEdge);
		gp_Pnt dividePntOnCurve = getMinimumDistancePoint(middleCurveEdge,areaDivideEdge);
		gp_Pnt bottomWidthLeftPntOnCurve = getMinimumDistancePoint(middleCurveEdge,bottomWidthLeftEdge);
		gp_Pnt bottomWidthRightPntOnCurve = getMinimumDistancePoint(middleCurveEdge,bottomWidthRightEdge);

		GeomLib_Tool::Parameter(middleCurve,leftStartOfSplitterPnt,1,uParaOfIntersectionPoint);
		Handle(Geom_BSplineCurve) leftPartOfMiddleCurve = Handle(Geom_BSplineCurve)::DownCast(middleCurve->Copy());
		leftPartOfMiddleCurve->Segment(leftTurningPntPara,uParaOfIntersectionPoint);

		GeomLib_Tool::Parameter(middleCurve,rightStartOfSplitterPnt,1,uParaOfIntersectionPoint);
		Handle(Geom_BSplineCurve) rightPartOfMiddleCurve = Handle(Geom_BSplineCurve)::DownCast(middleCurve->Copy());
		rightPartOfMiddleCurve->Segment(uParaOfIntersectionPoint,rightTurningPntPara);

		TopoDS_Edge leftPartOfMiddleCurveEdge=BRepBuilderAPI_MakeEdge(leftPartOfMiddleCurve);
		TopoDS_Edge rightPartOfMiddleCurveEdge=BRepBuilderAPI_MakeEdge(rightPartOfMiddleCurve);

		TopoDS_Edge lineToGetBottomLeftPointOfSplitter=BRepBuilderAPI_MakeEdge(bottomWidthLeftPntOnBase,bottomWidthLeftPntOnCurve);
		TopoDS_Edge lineToGetBottomRightPointOfSplitter=BRepBuilderAPI_MakeEdge(bottomWidthRightPntOnBase,bottomWidthRightPntOnCurve);

		double voluteHeight=dividePnt.Distance(dividePntOnCurve);
		double splitterHeight=voluteHeight-(voluteHeight*heightRatio);

		Handle_Geom_Curve curveToGetBottomLeftPointOfSplitter=BRep_Tool::Curve(lineToGetBottomLeftPointOfSplitter,basePara1,basePara2);
		Handle_Geom_Curve curveToGetBottomRightPointOfSplitter=BRep_Tool::Curve(lineToGetBottomRightPointOfSplitter,basePara1,basePara2);
		Handle_Geom_Curve curveToGetBottomCentrePointOfSplitter=BRep_Tool::Curve(areaDivideEdge,basePara1,basePara2);

		gp_Pnt leftBottomPointOfSplitter;
		gp_Pnt rightBottomPointOfSplitter;
		gp_Pnt centreBottomPointOfSplitter;

		curveToGetBottomLeftPointOfSplitter->D0(splitterHeight,leftBottomPointOfSplitter);
		curveToGetBottomRightPointOfSplitter->D0(splitterHeight,rightBottomPointOfSplitter);
		curveToGetBottomCentrePointOfSplitter->D0(splitterHeight,centreBottomPointOfSplitter);

		TopoDS_Edge leftSplitterEdge=BRepBuilderAPI_MakeEdge(leftStartOfSplitterPnt,leftBottomPointOfSplitter);
		TopoDS_Edge rightSplitterEdge=BRepBuilderAPI_MakeEdge(rightStartOfSplitterPnt,rightBottomPointOfSplitter);

		gp_Vec normalVec=uVectorSurface.Crossed(vVectorSurface);
		gp_Dir dir2(normalVec);
		gp_Ax2 axisForTipEllipse(centreBottomPointOfSplitter,dir2);

		gp_Elips tipEllipse(axisForTipEllipse,tipRadius/10,bottomWidth/20);//error when raduis<width/2

		Handle(Geom_TrimmedCurve) tipellipseArc =GC_MakeArcOfEllipse(tipEllipse,leftBottomPointOfSplitter,rightBottomPointOfSplitter,true);


		TopoDS_Vertex bottomVert1=BRepBuilderAPI_MakeVertex(leftBottomPointOfSplitter);
		TopoDS_Vertex bottomVert2=BRepBuilderAPI_MakeVertex(rightBottomPointOfSplitter);
		TopoDS_Vertex bottomVert3=BRepBuilderAPI_MakeVertex(centreBottomPointOfSplitter);

		TopoDS_Vertex leftStartOfSplitterVert=BRepBuilderAPI_MakeVertex(leftStartOfSplitterPnt);
		TopoDS_Vertex rightStartOfSplitterVert=BRepBuilderAPI_MakeVertex(rightStartOfSplitterPnt);

		TopoDS_Edge elipsEdge=BRepBuilderAPI_MakeEdge(tipellipseArc);

		TopoDS_Wire baseWire=BRepBuilderAPI_MakeWire(straitLeftBaseEdge,horizontalBaseEdge,straitRightBaseEdge);
		TopoDS_Edge divideEdge=BRepBuilderAPI_MakeEdge(dividePnt,dividePntOnCurve);


		TopoDS_Edge leftPartOfMiddleCurveFilletedEdge;
		TopoDS_Edge leftSplitterFilletedEdge;

		gp_Pln filletPlane(dividePnt,dir2);
		ChFi2d_FilletAPI leftFillet(leftPartOfMiddleCurveEdge,leftSplitterEdge,filletPlane);
		leftFillet.Perform(topRadius);
		TopoDS_Edge leftFilletEdge=leftFillet.Result(leftStartOfSplitterPnt,leftPartOfMiddleCurveFilletedEdge,leftSplitterFilletedEdge,-1);


		TopoDS_Edge rightPartOfMiddleCurveFilletedEdge;
		TopoDS_Edge rightSplitterFilletedEdge;


		ChFi2d_FilletAPI rightFillet(rightPartOfMiddleCurveEdge,rightSplitterEdge,filletPlane);
		rightFillet.Perform(topRadius);
		TopoDS_Edge rightFilletEdge=rightFillet.Result(rightStartOfSplitterPnt,rightPartOfMiddleCurveFilletedEdge,rightSplitterFilletedEdge,-1);


		BRepBuilderAPI_MakeWire crossSectionWire;

		crossSectionWire.Add(leftPartOfMiddleCurveFilletedEdge);
		crossSectionWire.Add(leftFilletEdge);
		crossSectionWire.Add(leftSplitterFilletedEdge);
		crossSectionWire.Add(elipsEdge);
		crossSectionWire.Add(rightSplitterFilletedEdge);
		crossSectionWire.Add(rightFilletEdge);
		crossSectionWire.Add(rightPartOfMiddleCurveFilletedEdge);
		crossSectionWire.Add(baseWire);

		TopoDS_Wire completeCrossSection=crossSectionWire;



		TopoDS_Compound compoundWire;
		BRep_Builder abuilder;
		abuilder.MakeCompound(compoundWire);
		abuilder.Add(compoundWire,divideEdge);
		abuilder.Add(compoundWire,completeCrossSection);

		BRepTools::Write( completeCrossSection,"D:/Breps/2ndwireSet.brep");


		//TopoDS_Edge leftFilletEdge=fillet.AddFillet(leftStartOfSplitterVert,0.03);
		/*CString str;
		CString str1;
		str1.Format(_T("volute height %g \n"), voluteHeight);
		str.Format(_T("volute height distance method %g \n"), voluteheight2);
		AfxMessageBox(str1+str);*/



		TopoDS_Wire testWire=BRepBuilderAPI_MakeWire(areaDivideEdge);

		gp_Pnt divideIntersectionPointOnCurve=getMinimumDistancePoint(importedEdge,areaDivideEdge);
		TopoDS_Vertex intersectionVertexOnCurve=BRepBuilderAPI_MakeVertex(divideIntersectionPointOnCurve);




		//m_pcoloredshapeList->Add(Quantity_NOC_RED,rightSplitterStartEdge);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,straitLeftBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,straitRightBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftPartOfMiddleCurveEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,rightPartOfMiddleCurveEdge);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,bottomVert3);
		//m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,leftSplitterEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,rightSplitterEdge);


		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,elipsEdge);

		////m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,testWire);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,horizontalBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,intersectionVertexOnCurve);

		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,completeCrossSection);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,baseWire);



		//m_pcoloredshapeList->Add(Quantity_NOC_RED,rightSplitterFilletedEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftSplitterFilletedEdge);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();
		return completeCrossSection;

	}




	TopoDS_Wire CImportExportDoc::createSplitter(double dividePointPara,double topWidth,double bottomWidth,double heightRatio,TopoDS_Wire importedWire )
	{

		double leftStartParaOfSplitterOnBase;
		double rightStartParaOfSplitterOnBase;

		double bottomWidthLeftParaOnBase;
		double bottomWidthRightParaOnBase;

		double tipRadius=15;
		double topRadius=3;
		leftStartParaOfSplitterOnBase=dividePointPara-(topWidth/20);
		rightStartParaOfSplitterOnBase=dividePointPara+(topWidth/20);

		bottomWidthLeftParaOnBase=dividePointPara-(bottomWidth/20);
		bottomWidthRightParaOnBase=dividePointPara+(bottomWidth/20);

		TopoDS_Edge edges[4];
		TopExp_Explorer anEdgeExplorer(importedWire, TopAbs_EDGE);
		int i=0;
		while(anEdgeExplorer.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			edges[i]=anEdge;
			anEdgeExplorer.Next();
			i++;

		}

		TopoDS_Edge importedEdge=edges[0];
		Standard_Real U1=0;
		Standard_Real U2=1;
		Standard_Real uParaOfIntersectionPoint;

		gp_Pnt pnt;
		gp_Vec V1;
		gp_Vec V2;
		gp_Vec V3;
		Handle_Geom_Curve curve= BRep_Tool::Curve(importedEdge,U1, U2);

		TopoDS_Vertex vertex[2];
		TopExp_Explorer vertexExplorer(importedWire,TopAbs_VERTEX);
		int j=0;
		while(vertexExplorer.More()){
			TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
			vertex[j]=aVertex;
			vertexExplorer.Next();
			j++;

		}

		Standard_Real basePara1;
		Standard_Real basePara2;
		gp_Pnt dividePnt;

		gp_Pnt leftStartOfSplitterOnBasePnt;
		gp_Pnt rightStartOfSplitterOnBasePnt;
		gp_Pnt bottomWidthLeftPntOnBase;
		gp_Pnt bottomWidthRightPntOnBase;


		gp_Vec devidePointVec;

		TopoDS_Edge horizontalBaseEdge=BRepBuilderAPI_MakeEdge(vertex[0],vertex[1]);
		Handle_Geom_Curve horizontalBaseLine=BRep_Tool::Curve(horizontalBaseEdge,basePara1,basePara2);

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;


		BRepFill_Filling filledFace;
		filledFace.Add(edges[0],GeomAbs_C0);
		filledFace.Add(horizontalBaseEdge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,dividePnt,uVectorSurface,vVectorSurface);

		horizontalBaseLine->D1(dividePointPara,dividePnt,devidePointVec);
		horizontalBaseLine->D0(leftStartParaOfSplitterOnBase,leftStartOfSplitterOnBasePnt);
		horizontalBaseLine->D0(rightStartParaOfSplitterOnBase,rightStartOfSplitterOnBasePnt);	
		horizontalBaseLine->D0(bottomWidthLeftParaOnBase,bottomWidthLeftPntOnBase);
		horizontalBaseLine->D0(bottomWidthRightParaOnBase,bottomWidthRightPntOnBase);


		BRepBuilderAPI_MakeVertex divideVertex(dividePnt);

		Standard_Real leftTurningPntPara=U2*(44/U2);
		Standard_Real rightTurningPntPara=U2*(44/U2)*3;

		Handle(Geom_BSplineCurve) straitLeft = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) middleCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) straitRight = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());

		straitLeft->Segment(U1,leftTurningPntPara);
		middleCurve->Segment(leftTurningPntPara,rightTurningPntPara);
		straitRight->Segment(rightTurningPntPara,U2);

		gp_Dir dir(uVectorSurface); 
		gp_Lin areaDivideLine(dividePnt,dir);
		gp_Lin leftSplitterStartLine(leftStartOfSplitterOnBasePnt,dir);
		gp_Lin rightSplitterStartLine(rightStartOfSplitterOnBasePnt,dir);
		gp_Lin bottomWidthLeftLine(bottomWidthLeftPntOnBase,dir);
		gp_Lin bottomWidthRightLine(bottomWidthRightPntOnBase,dir);

		TopoDS_Edge straitLeftBaseEdge=BRepBuilderAPI_MakeEdge(straitLeft);
		TopoDS_Edge straitRightBaseEdge=BRepBuilderAPI_MakeEdge(straitRight);
		TopoDS_Edge middleCurveEdge=BRepBuilderAPI_MakeEdge(middleCurve);

		TopoDS_Edge areaDivideEdge= BRepBuilderAPI_MakeEdge(areaDivideLine);
		TopoDS_Edge leftSplitterStartEdge=BRepBuilderAPI_MakeEdge(leftSplitterStartLine);
		TopoDS_Edge rightSplitterStartEdge=BRepBuilderAPI_MakeEdge(rightSplitterStartLine);
		TopoDS_Edge bottomWidthLeftEdge=BRepBuilderAPI_MakeEdge(bottomWidthLeftLine);
		TopoDS_Edge bottomWidthRightEdge=BRepBuilderAPI_MakeEdge(bottomWidthRightLine);

		gp_Pnt leftStartOfSplitterPnt=getMinimumDistancePoint(importedEdge,leftSplitterStartEdge);
		gp_Pnt rightStartOfSplitterPnt=getMinimumDistancePoint(importedEdge,rightSplitterStartEdge);
		gp_Pnt dividePntOnCurve=getMinimumDistancePoint(importedEdge,areaDivideEdge);
		gp_Pnt bottomWidthLeftPntOnCurve=getMinimumDistancePoint(importedEdge,bottomWidthLeftEdge);
		gp_Pnt bottomWidthRightPntOnCurve=getMinimumDistancePoint(importedEdge,bottomWidthRightEdge);

		GeomLib_Tool::Parameter(curve,leftStartOfSplitterPnt,1,uParaOfIntersectionPoint);
		Handle(Geom_BSplineCurve) leftPartOfMiddleCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		leftPartOfMiddleCurve->Segment(leftTurningPntPara,uParaOfIntersectionPoint);

		GeomLib_Tool::Parameter(curve,rightStartOfSplitterPnt,1,uParaOfIntersectionPoint);
		Handle(Geom_BSplineCurve) rightPartOfMiddleCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		rightPartOfMiddleCurve->Segment(uParaOfIntersectionPoint,rightTurningPntPara);

		TopoDS_Edge leftPartOfMiddleCurveEdge=BRepBuilderAPI_MakeEdge(leftPartOfMiddleCurve);
		TopoDS_Edge rightPartOfMiddleCurveEdge=BRepBuilderAPI_MakeEdge(rightPartOfMiddleCurve);

		TopoDS_Edge lineToGetBottomLeftPointOfSplitter=BRepBuilderAPI_MakeEdge(bottomWidthLeftPntOnBase,bottomWidthLeftPntOnCurve);
		TopoDS_Edge lineToGetBottomRightPointOfSplitter=BRepBuilderAPI_MakeEdge(bottomWidthRightPntOnBase,bottomWidthRightPntOnCurve);

		double voluteHeight=dividePnt.Distance(dividePntOnCurve);
		double splitterHeight=voluteHeight-(voluteHeight*heightRatio);

		Handle_Geom_Curve curveToGetBottomLeftPointOfSplitter=BRep_Tool::Curve(lineToGetBottomLeftPointOfSplitter,basePara1,basePara2);
		Handle_Geom_Curve curveToGetBottomRightPointOfSplitter=BRep_Tool::Curve(lineToGetBottomRightPointOfSplitter,basePara1,basePara2);
		Handle_Geom_Curve curveToGetBottomCentrePointOfSplitter=BRep_Tool::Curve(areaDivideEdge,basePara1,basePara2);

		gp_Pnt leftBottomPointOfSplitter;
		gp_Pnt rightBottomPointOfSplitter;
		gp_Pnt centreBottomPointOfSplitter;

		curveToGetBottomLeftPointOfSplitter->D0(splitterHeight,leftBottomPointOfSplitter);
		curveToGetBottomRightPointOfSplitter->D0(splitterHeight,rightBottomPointOfSplitter);
		curveToGetBottomCentrePointOfSplitter->D0(splitterHeight,centreBottomPointOfSplitter);

		TopoDS_Edge leftSplitterEdge=BRepBuilderAPI_MakeEdge(leftStartOfSplitterPnt,leftBottomPointOfSplitter);
		TopoDS_Edge rightSplitterEdge=BRepBuilderAPI_MakeEdge(rightStartOfSplitterPnt,rightBottomPointOfSplitter);

		gp_Vec normalVec=uVectorSurface.Crossed(vVectorSurface);
		gp_Dir dir2(normalVec);
		gp_Ax2 axisForTipEllipse(centreBottomPointOfSplitter,dir2);

		gp_Elips tipEllipse(axisForTipEllipse,tipRadius/10,bottomWidth/20);//error when raduis<width/2

		Handle(Geom_TrimmedCurve) tipellipseArc =GC_MakeArcOfEllipse(tipEllipse,leftBottomPointOfSplitter,rightBottomPointOfSplitter,true);


		TopoDS_Vertex bottomVert1=BRepBuilderAPI_MakeVertex(leftBottomPointOfSplitter);
		TopoDS_Vertex bottomVert2=BRepBuilderAPI_MakeVertex(rightBottomPointOfSplitter);
		TopoDS_Vertex bottomVert3=BRepBuilderAPI_MakeVertex(centreBottomPointOfSplitter);

		TopoDS_Vertex leftStartOfSplitterVert=BRepBuilderAPI_MakeVertex(leftStartOfSplitterPnt);
		TopoDS_Vertex rightStartOfSplitterVert=BRepBuilderAPI_MakeVertex(rightStartOfSplitterPnt);

		TopoDS_Edge elipsEdge=BRepBuilderAPI_MakeEdge(tipellipseArc);

		TopoDS_Wire baseWire=BRepBuilderAPI_MakeWire(straitLeftBaseEdge,horizontalBaseEdge,straitRightBaseEdge);
		TopoDS_Edge divideEdge=BRepBuilderAPI_MakeEdge(dividePnt,dividePntOnCurve);


		TopoDS_Edge leftPartOfMiddleCurveFilletedEdge;
		TopoDS_Edge leftSplitterFilletedEdge;

		gp_Pln filletPlane(dividePnt,dir2);
		ChFi2d_FilletAPI leftFillet(leftPartOfMiddleCurveEdge,leftSplitterEdge,filletPlane);
		leftFillet.Perform(topRadius);
		TopoDS_Edge leftFilletEdge=leftFillet.Result(leftStartOfSplitterPnt,leftPartOfMiddleCurveFilletedEdge,leftSplitterFilletedEdge,-1);


		TopoDS_Edge rightPartOfMiddleCurveFilletedEdge;
		TopoDS_Edge rightSplitterFilletedEdge;


		ChFi2d_FilletAPI rightFillet(rightPartOfMiddleCurveEdge,rightSplitterEdge,filletPlane);
		rightFillet.Perform(topRadius);
		TopoDS_Edge rightFilletEdge=rightFillet.Result(rightStartOfSplitterPnt,rightPartOfMiddleCurveFilletedEdge,rightSplitterFilletedEdge,-1);


		BRepBuilderAPI_MakeWire crossSectionWire;

		crossSectionWire.Add(leftPartOfMiddleCurveFilletedEdge);
		crossSectionWire.Add(leftFilletEdge);
		crossSectionWire.Add(leftSplitterFilletedEdge);
		crossSectionWire.Add(elipsEdge);
		crossSectionWire.Add(rightSplitterFilletedEdge);
		crossSectionWire.Add(rightFilletEdge);
		crossSectionWire.Add(rightPartOfMiddleCurveFilletedEdge);
		crossSectionWire.Add(baseWire);

		TopoDS_Wire completeCrossSection=crossSectionWire;



		TopoDS_Compound compoundWire;
		BRep_Builder abuilder;
		abuilder.MakeCompound(compoundWire);
		abuilder.Add(compoundWire,divideEdge);
		abuilder.Add(compoundWire,completeCrossSection);

		BRepTools::Write( completeCrossSection,"D:/Breps/2ndwireSet.brep");


		//TopoDS_Edge leftFilletEdge=fillet.AddFillet(leftStartOfSplitterVert,0.03);
		/*CString str;
		CString str1;
		str1.Format(_T("volute height %g \n"), voluteHeight);
		str.Format(_T("volute height distance method %g \n"), voluteheight2);
		AfxMessageBox(str1+str);*/



		TopoDS_Wire testWire=BRepBuilderAPI_MakeWire(areaDivideEdge);

		gp_Pnt divideIntersectionPointOnCurve=getMinimumDistancePoint(importedEdge,areaDivideEdge);
		TopoDS_Vertex intersectionVertexOnCurve=BRepBuilderAPI_MakeVertex(divideIntersectionPointOnCurve);




		//m_pcoloredshapeList->Add(Quantity_NOC_RED,rightSplitterStartEdge);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,straitLeftBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,straitRightBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftPartOfMiddleCurveEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,rightPartOfMiddleCurveEdge);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,bottomVert3);
		//m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,leftSplitterEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,rightSplitterEdge);


		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,elipsEdge);

		////m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,testWire);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,horizontalBaseEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,intersectionVertexOnCurve);

		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,completeCrossSection);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,baseWire);



		//m_pcoloredshapeList->Add(Quantity_NOC_RED,rightSplitterFilletedEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,leftSplitterFilletedEdge);
		m_pcoloredshapeList->Display(myAISContext);
		Fit();
		return completeCrossSection;

	}

	TopoDS_Wire CImportExportDoc::createLeftPartOfDualVoluteWire(double ratio,double dividePointUParameter,TopoDS_Wire importedWire){


		TopoDS_Edge edges[4];
		TopExp_Explorer anEdgeExplorer(importedWire, TopAbs_EDGE);
		int i=0;
		while(anEdgeExplorer.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			edges[i]=anEdge;
			anEdgeExplorer.Next();
			i++;

		}

		TopoDS_Edge importedEdge=edges[0];

		TopoDS_Vertex vertex[2];
		TopExp_Explorer vertexExplorer(importedEdge,TopAbs_VERTEX);
		int j=0;
		while(vertexExplorer.More()){
			TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
			vertex[j]=aVertex;
			vertexExplorer.Next();
			j++;
		}

		Standard_Real basePara1;
		Standard_Real basePara2;
		gp_Pnt dividePnt;
		gp_Vec devidePointVec;

		TopoDS_Edge horizontalBaseEdge=BRepBuilderAPI_MakeEdge(vertex[0],vertex[1]);
		Handle_Geom_Curve horizontalBaseLine=BRep_Tool::Curve(horizontalBaseEdge,basePara1,basePara2);

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;

		BRepFill_Filling filledFace;
		filledFace.Add(importedEdge,GeomAbs_C0);
		filledFace.Add(horizontalBaseEdge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,dividePnt,uVectorSurface,vVectorSurface);

		horizontalBaseLine->D1(dividePointUParameter,dividePnt,devidePointVec);

		BRepBuilderAPI_MakeVertex divideVertex(dividePnt);
		gp_Dir dir(uVectorSurface); 
		gp_Lin testLine(dividePnt,dir);

		TopoDS_Edge testEdge= BRepBuilderAPI_MakeEdge(testLine);

		TopoDS_Wire DivideWire=BRepBuilderAPI_MakeWire(testEdge);
		TopoDS_Shape testShape=DivideWire;

		BRepExtrema_DistShapeShape minimumDist(importedWire,DivideWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);

		gp_Pnt intersectionPointOnCurve=minimumDist.PointOnShape1(1);

		TopoDS_Vertex intersectionVertexOnCurve=BRepBuilderAPI_MakeVertex(intersectionPointOnCurve);

		TopoDS_Edge middleEdge=BRepBuilderAPI_MakeEdge(dividePnt,intersectionPointOnCurve);

		return importedWire;
	}

	void CImportExportDoc::OnDualVolute()
	{

		double leftAreaPercentage;
		double rightAreaPercentage;
		double ratio;
		double wholeWireSurfaceArea;
		double leftWireExpectedArea;
		double rightWireExpectedArea;

		leftAreaPercentage=20;
		rightAreaPercentage=80;

		ratio=leftAreaPercentage/100;

		TopoDS_Shape importedShape;

		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadBREP();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			importedShape=aSeqOfShape->Value(i);

		}

		TopoDS_Wire importedWire=TopoDS::Wire(importedShape);

		TopoDS_Edge edges[4];
		TopExp_Explorer anEdgeExplorer(importedWire, TopAbs_EDGE);
		int i=0;
		while(anEdgeExplorer.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			edges[i]=anEdge;
			anEdgeExplorer.Next();
			i++;

		}

		TopoDS_Vertex vertex[2];
		TopExp_Explorer vertexExplorer(importedWire,TopAbs_VERTEX);
		int j=0;
		while(vertexExplorer.More()){
			TopoDS_Vertex aVertex =TopoDS::Vertex(vertexExplorer.Current());
			vertex[j]=aVertex;
			vertexExplorer.Next();
			j++;

		}

		Standard_Real basePara1;
		Standard_Real basePara2;
		gp_Pnt dividePnt;
		gp_Vec devidePointVec;

		TopoDS_Edge horizontalBaseEdge=BRepBuilderAPI_MakeEdge(vertex[0],vertex[1]);
		Handle_Geom_Curve horizontalBaseLine=BRep_Tool::Curve(horizontalBaseEdge,basePara1,basePara2);

		Standard_Real uParaSurface;
		Standard_Real vParaSurface;
		gp_Vec uVectorSurface;
		gp_Vec vVectorSurface;

		BRepFill_Filling filledFace;
		filledFace.Add(edges[0],GeomAbs_C0);
		filledFace.Add(horizontalBaseEdge,GeomAbs_C0);
		filledFace.Build();
		TopoDS_Face fillFace=filledFace.Face();

		BRepAdaptor_Surface aface(fillFace);
		aface.D1(uParaSurface,vParaSurface,dividePnt,uVectorSurface,vVectorSurface);

		horizontalBaseLine->D1(basePara2*ratio,dividePnt,devidePointVec);
		BRepBuilderAPI_MakeVertex divideVertex(dividePnt);
		gp_Dir dir(1,0,0); 
		gp_Lin testLine(dividePnt,dir);
		TopoDS_Edge testEdge= BRepBuilderAPI_MakeEdge(testLine);

		TopoDS_Wire testWire=BRepBuilderAPI_MakeWire(testEdge);
		TopoDS_Shape testShape=testWire;

		gp_Dir xDir(0,1,0);
		gp_Ax1 xAxis(dividePnt, xDir);
		//BRepPrimAPI_MakeRevol revol(testShape,xAxis);  
		//revol.Build();

		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,revol.Shape());



		BRepExtrema_DistShapeShape minimumDist(importedWire,testWire,Extrema_ExtFlag_MINMAX,Extrema_ExtAlgo_Grad);

		gp_Pnt intersectionPointOnCurve=minimumDist.PointOnShape1(1);

		TopoDS_Vertex intersectionVertexOnCurve=BRepBuilderAPI_MakeVertex(intersectionPointOnCurve);


		TopoDS_Wire splitteWire=createSplitter(basePara2*ratio,60,20,0.3,importedWire);

		//Standard_Real distValue=minimumDist.DistValue();
		//TopoDS_Edge importedEdge=TopoDS::Edge(importedWire);

		//double U1=0;
		//double U2=1;
		Standard_Real U1=0;
		Standard_Real U2=1;
		Standard_Real uParaOfIntersectionPoint;

		gp_Pnt pnt;
		gp_Vec V1;
		gp_Vec V2;
		gp_Vec V3;
		Handle_Geom_Curve curve= BRep_Tool::Curve(edges[0],U1, U2);


		GeomAPI_ProjectPointOnCurve projection(intersectionPointOnCurve,curve);

		Quantity_Parameter uParaOfIntersection=projection.Parameter(1);

		GeomLib_Tool::Parameter(curve,intersectionPointOnCurve,1,uParaOfIntersectionPoint);

		//	Handle_Geom_BSplineCurve splineCurve=BRep_Tool::Curve(edges[0],U1,U2);

		Handle(Geom_BSplineCurve) splineCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) splineBase1 = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());

		Handle(Geom_BSplineCurve) splineBase2 = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) leftSegmentOfCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) rightSegmentOfCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());



		//Handle_Geom_BSplineCurve splineCurve= Handle_Geom_BSplineCurve::DownCast(curve->Copy());

		splineCurve->Segment(47,131);

		splineBase1->Segment(0,47);

		splineBase2->Segment(131,U2);

		leftSegmentOfCurve->Segment(U1,uParaOfIntersectionPoint);
		rightSegmentOfCurve->Segment(uParaOfIntersectionPoint,U2);

		TopoDS_Edge splineEdge= BRepBuilderAPI_MakeEdge(splineCurve);
		TopoDS_Edge baseEdge1= BRepBuilderAPI_MakeEdge(splineBase1);
		TopoDS_Edge baseEdge2=BRepBuilderAPI_MakeEdge(splineBase2);



		Standard_Real leftTurningPntPara=U2*(44/U2);
		Standard_Real rightTurningPntPara=U2*(44/U2)*3;



		Handle(Geom_BSplineCurve) straitLeft = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) middleCurve = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
		Handle(Geom_BSplineCurve) straitRight = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());

		straitLeft->Segment(U1,leftTurningPntPara);
		middleCurve->Segment(leftTurningPntPara,rightTurningPntPara);
		straitRight->Segment(rightTurningPntPara,U2);

		TopoDS_Edge straitLeftEdge=BRepBuilderAPI_MakeEdge(straitLeft);
		TopoDS_Edge straitRightEdge=BRepBuilderAPI_MakeEdge(straitRight);
		TopoDS_Edge middleCurveEdge=BRepBuilderAPI_MakeEdge(middleCurve);

		//TopoDS_Edge leftSegmentOfCurveEdge=BRepBuilderAPI_MakeEdge(leftSegmentOfCurve);
		//TopoDS_Edge rightSegmentOfCurveEdge=BRepBuilderAPI_MakeEdge(rightSegmentOfCurve);

		curve->D0(uParaOfIntersectionPoint,intersectionPointOnCurve);

		CString str;
		//str.Format(_T("dist %g \n"),distValue);
		CString str1;
		str1.Format(_T("upara %g \n"), uParaOfIntersectionPoint);

		AfxMessageBox(str+str1);


		ofstream out;
		out.open("D:/text.txt");

		for(int i=0; i<U2;i++){

			Standard_Real U=i;
			curve->D1(U,pnt,V1);
			Standard_Real x=V1.X();

			out<<"\n  U   : "<<U<<endl;
			out<<" v1 : "<<x<<endl;
		}


		//TopoDS_Edge leFaceTopBottomEdge = TopoDS::Edge(Expl.Current());
		//Handle_Geom_Curve cleFaceBottom = BRep_Tool::Curve(leFaceTopBottomEdge,U1, U2);
		//SurfaceLine leFaceTopBottomSurfaceLine = faceBuilder.makeSurfaceLineFromEdge(leFaceTopBottomEdge);
		//Node nodeTemp;
		//if(data.machineType != 1)
		//	/*double zVal*/nodeTemp = leFaceTopBottomSurfaceLine.minZPoint();//minZ();
		//else
		//	nodeTemp = leFaceTopBottomSurfaceLine.maxZPoint();;
		////Node nodeTemp = leFaceTopBottomSurfaceLine.atZ(zVal);
		//TopoDS_Vertex vminZ;
		//gp_Pnt projectedMinZPoint;
		//gp_Pnt minZPoint(nodeTemp.lp().cp().x(),nodeTemp.lp().cp().y(),nodeTemp.lp().cp().z());
		////gp_Pnt minZPoint = cleFaceBottom->Value((U1+U2)/2);
		//GeomAPI_ProjectPointOnCurve pointProject2(minZPoint,cleFaceBottom);
		//Standard_Real vProjPoint1,vProjPoint2;
		//if(pointProject2.NbPoints()>0)
		//{


		Standard_Real turningPnt=U2*(44/U2);
		Standard_Real turningPnt2=U2*(44/U2)*3;

		curve->D1(turningPnt,pnt,V1);
		BRepBuilderAPI_MakeVertex vert(pnt);

		curve->D1(turningPnt2,pnt,V1);

		BRepBuilderAPI_MakeVertex vert2(pnt);

		//m_pcoloredshapeList->Add(Quantity_NOC_BLACK,edges[0]);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,splineEdge);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,baseEdge1);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,vert);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,vert2);
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,horizontalBaseEdge);

		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,importedWire);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,leftSegmentOfCurveEdge);


		m_pcoloredshapeList->Display(myAISContext);
		Fit();

	}



	void CImportExportDoc::OnImportVolute()
	{

		TopoDS_Compound compoundShape;
		//BRep_Builder aBuilder;
		TopoDS_Shape shapeArray [10];
		std::vector<TopoDS_Shape> shapes;
		TopoDS_Shape importedShape;

		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadBREP();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			importedShape=aSeqOfShape->Value(i);

		}

		BRepPrimAPI_MakeBox(15,10,21);

		TopoDS_Face faces [50];
		BRepBuilderAPI_Sewing connectedFaces(1.0e-06,true,true,true,true);
		BRepBuilderAPI_Sewing originalVolute(1.0e-06,true,true,true,true);



		int i=0;
		for(TopExp_Explorer aFaceExploree(importedShape,TopAbs_FACE);aFaceExploree.More();aFaceExploree.Next()){

			TopoDS_Face aFace= TopoDS::Face(aFaceExploree.Current());
			faces[i]=aFace;

			if(i<4){
				connectedFaces.Add(aFace);
			}
			if(i==4){
				connectedFaces.Add(aFace);
			}
			if(i==6){
				connectedFaces.Add(aFace);
			}
			if(i==9){
				connectedFaces.Add(aFace);
			}
			if(i==14){
				connectedFaces.Add(aFace);
			}
			if(i==15){
				connectedFaces.Add(aFace);
			}
			if(i==16){
				connectedFaces.Add(aFace);
			}
			if(i==18){
				connectedFaces.Add(aFace);
			}
			if(i==19){
				connectedFaces.Add(aFace);
			}
			if(i==17){
				connectedFaces.Add(aFace);
			}

			i++;
		}

		for(int i=0;i<20;i++){
			originalVolute.Add(faces[i]);

		}

		for(TopExp_Explorer aFaceExploree(importedShape,TopAbs_FACE);aFaceExploree.More();aFaceExploree.Next()){

			TopoDS_Face aFace= TopoDS::Face(aFaceExploree.Current());

		}

		connectedFaces.Perform();
		TopoDS_Shape connectedShapes=connectedFaces.SewedShape();
		GeomAbs_JoinType joinType = GeomAbs_Arc;
		//BRepOffsetAPI_MakeOffsetShape exitOffSet(faces[6],-0.01,1e-6,BRepOffset_Skin,false,false,joinType);

		double offSetThickness=-0.03;
		BRepOffsetAPI_MakeOffsetShape voluteOffSet(connectedShapes,offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet13(faces[13],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet12(faces[12],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet10(faces[10],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet11(faces[11],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet7(faces[7],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet8(faces[8],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet4(faces[4],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);
		BRepOffsetAPI_MakeOffsetShape curvedOffSet5(faces[5],offSetThickness,1e-6,BRepOffset_Skin,false,false,joinType);


		//m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet4);
		//		m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet8);
		m_pcoloredshapeList->Add(Quantity_NOC_IVORY,faces[8]);


		TopExp_Explorer anEdgeExplorer12(curvedOffSet12, TopAbs_EDGE);
		TopoDS_Edge edges12[4];
		int iE12=0;
		while(anEdgeExplorer12.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer12.Current());
			edges12[iE12]=anEdge;
			anEdgeExplorer12.Next();
			iE12++;

		}
		TopoDS_Edge edges11[4];
		TopExp_Explorer anEdgeExplorer11(curvedOffSet11, TopAbs_EDGE);
		int iE11=0;
		while(anEdgeExplorer11.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer11.Current());
			edges11[iE11]=anEdge;
			anEdgeExplorer11.Next();
			iE11++;
		}
		TopoDS_Edge edges13[4];
		TopExp_Explorer anEdgeExplorer13(curvedOffSet13, TopAbs_EDGE);
		int iE13=0;
		while(anEdgeExplorer13.More()){
			TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer13.Current());
			edges13[iE13]=anEdge;
			anEdgeExplorer13.Next();
			iE13++;

		}
		TopoDS_Vertex vertex12 [2];
		int iV12=0;
		for(TopExp_Explorer expVert(edges12[3], TopAbs_VERTEX); expVert.More(); expVert.Next())
		{
			TopoDS_Vertex vert = TopoDS::Vertex(expVert.Current());
			vertex12[iV12]=vert;
			iV12++;
		}
		TopoDS_Vertex vertex11_12 [2];
		int iV11=0;
		for(TopExp_Explorer expVert(edges11[2], TopAbs_VERTEX); expVert.More(); expVert.Next())
		{
			TopoDS_Vertex vert = TopoDS::Vertex(expVert.Current());
			vertex11_12[iV11]=vert;
			iV11++;
		}


		TopoDS_Vertex vertex11_13 [2];
		int iV11_2=0;
		for(TopExp_Explorer expVert(edges11[0], TopAbs_VERTEX); expVert.More(); expVert.Next())
		{
			TopoDS_Vertex vert = TopoDS::Vertex(expVert.Current());
			vertex11_13[iV11_2]=vert;
			iV11_2++;
		}

		TopoDS_Vertex vertex13 [2];
		int iV13=0;
		for(TopExp_Explorer expVert(edges13[1],TopAbs_VERTEX);expVert.More();expVert.Next())
		{
			TopoDS_Vertex vert=TopoDS::Vertex(expVert.Current());
			vertex13[iV13]=vert;
			iV13++;
		}

		TopoDS_Edge edge13_11=BRepBuilderAPI_MakeEdge(vertex11_13[0],vertex13[0]);
		TopoDS_Edge edge12_11= BRepBuilderAPI_MakeEdge(vertex11_12[0],vertex12[0]);

		TopoDS_Wire triangle12=BRepBuilderAPI_MakeWire(edges11[2],edges12[3],edge12_11);
		TopoDS_Wire triangle13=BRepBuilderAPI_MakeWire(edges11[0],edges13[1],edge13_11);

		TopoDS_Face triangle12Face=BRepBuilderAPI_MakeFace(triangle12);
		TopoDS_Face triangle13Face=BRepBuilderAPI_MakeFace(triangle13);

		BRepBuilderAPI_Sewing manualyOffSetted(1.0e-06,true,true,true,true);
		manualyOffSetted.Add(curvedOffSet12);
		manualyOffSetted.Add(curvedOffSet11);
		manualyOffSetted.Add(curvedOffSet13);
		manualyOffSetted.Add(triangle12Face);
		manualyOffSetted.Add(triangle13Face);

		manualyOffSetted.Perform();
		originalVolute.Perform();

		m_pcoloredshapeList->Add(Quantity_NOC_RED,manualyOffSetted.SewedShape());

		TopoDS_Shape manualyOffsettedShape=manualyOffSetted.SewedShape();
		TopoDS_Shape voluteOffSet1=voluteOffSet.Shape();
		TopoDS_Compound shapeWithThickness;

		BRep_Builder aBuilder;
		aBuilder.MakeCompound (shapeWithThickness);
		aBuilder.Add(shapeWithThickness,voluteOffSet1);
		aBuilder.Add(shapeWithThickness,manualyOffsettedShape);
		aBuilder.Add(shapeWithThickness,originalVolute.SewedShape());

		BRepTools::Write( shapeWithThickness,"D:/Breps/thickness.brep");

		//m_pcoloredshapeList->Add(Quantity_NOC_BLACK,edges13[0]);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,triangleFace);

		//BRepTools_ReShape 
		//m_pcoloredshapeList->Add(Quantity_NOC_GREEN,importedShape);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,exitOffSet);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,voluteOffSet.Shape());
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,connectedShapes);

		//m_pcoloredshapeList->Add(Quantity_NOC_RED,faces[17]);
		/*m_pcoloredshapeList->Add(Quantity_NOC_RED,faces[14]);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,faces[15]);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,faces[16]);*/

		//m_pcoloredshapeList->Add(Quantity_NOC_GOLD,faces[4]);
		//m_pcoloredshapeList->Add(Quantity_NOC_GOLD,faces[5]);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet10);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet11);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet8);
		//m_pcoloredshapeList->Add(Quantity_NOC_RED,curvedOffSet7);

		//m_pcoloredshapeList->Add(Quantity_NOC_ORANGE,faces[7]);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,faces[13]);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,faces[13]);
		//m_pcoloredshapeList->Add(Quantity_NOC_IVORY,faces[10]);
		/*m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,curvedOffSet13);
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,curvedOffSet10);
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,curvedOffSet11);*/

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

	}




	double CImportExportDoc::getSurfaceArea(TopoDS_Wire wireShape){

		BRepBuilderAPI_MakeFace wireFace(wireShape);
		TopoDS_Face face =wireFace;
		BRepGProp gprop; 
		GProp_GProps surfaceProps;



		gprop.SurfaceProperties(wireFace,surfaceProps);

		double surfaceArea=surfaceProps.Mass();
		CString fivs;

		//AfxMessageBox(fivs);

		return surfaceArea;

	}


	double CImportExportDoc::getTrapezuimHeight(double area,double width, double ang1Tan, double ang2Tan){



		double c=area*2;
		double a=1/ang1Tan+1/ang2Tan;
		double b= 2*width;

		double determinaterPart =(b*b) + (4*a*c); 
		double diterminator= sqrt(determinaterPart);
		double root1=((-b)+(diterminator))/(2*a);




		return root1;


	}

	void CImportExportDoc::OnVolute(){


		voluteDlg = new CVoluteDialog(myAISContext);
		voluteDlg->Create(IDD_DIALOG_Volute);
		voluteDlg->ShowWindow(SW_SHOW);


	}


	void CImportExportDoc::OnOffSet()
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
		double thickness=1;


		gp_Pnt r1(width/2.,height*2.,0);
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);
		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);
		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);


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
		Fit();

		//for 50% piece

		double initialTrapaziumArea = getSurfaceArea(trapeziumWire);
		double initailWholeArea= getSurfaceArea(initialwholeWire);
		double expectedWholeArea = initailWholeArea*(degree_0_percentage/100);//Expected whole area

		double intialCurveArea=getSurfaceArea(initialCurveWire);
		double ExpectedTrapeziumArea = initialTrapaziumArea*(degree_0_percentage/100);//user percentage of the 0 point area
		double newTrapeziumHeight= getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);//calculating the height
		double ExpectedCurveArea = intialCurveArea*(degree_0_percentage/100);
		TopoDS_Wire newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);
		//pcoloredshapeList->Add(Quantity_NOC_RED,initialwholeWire);
		//pcoloredshapeList->Add(Quantity_NOC_RED,newWholeWire);



		TopoDS_Wire scaledNewWire;
		TopoDS_Shape rotatedCurve;
		TopoDS_Shape rotatedShape;
		TopoDS_Shape rotatedBase;
		TopoDS_Face curveFace;
		TopoDS_Wire trapaziumBaseWire;
		gp_Trsf transfer;
		BRepOffsetAPI_ThruSections sections;
		BRepOffsetAPI_ThruSections scaledSections;
		BRepOffsetAPI_ThruSections curveSections;
		BRepOffsetAPI_ThruSections baseSections;

		gp_Ax1 axis(gp_Pnt(0,-15,0),gp_Dir(1,0,0));
		//initial wire rotation
		scaledNewWire=createOuterShell(height,width,ang1Tan,ang2Tan,thickness);
		//pcoloredshapeList->Add(Quantity_NOC_GREEN,scaledNewWire);


		double ang=360;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated(newWholeWire,transfer);
		rotatedShape=rotated.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);


		scaledSections.AddWire(scaledNewWire);
		double areaIncreasingFactor=(100-degree_0_percentage)/18;
		CString areaIncreacingFactor;


		for(int i=2;i<=18;i++){

			double ang=20*(i-1);
			transfer.SetRotation(axis,ang*PI/180);


			double expectedWholeArea = initailWholeArea*((areaIncreasingFactor/100)*i);
			double ExpectedTrapeziumArea = initialTrapaziumArea*((areaIncreasingFactor/100)*i);
			double newTrapeziumHeight=getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);
			newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);
			trapaziumBaseWire=createNewTrapazium(newTrapeziumHeight,width,ang1Tan,ang2Tan);

			TopoDS_Wire curveWire=createCurveEdge(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea);
			//TopoDS_Wire curveWire =TopoDS::Wire(newCurve);			

			BRepBuilderAPI_Transform baseRotated(trapaziumBaseWire,transfer);
			rotatedBase=baseRotated.Shape();
			trapaziumBaseWire=TopoDS::Wire(rotatedBase);
			baseSections.AddWire(trapaziumBaseWire);

			BRepBuilderAPI_Transform curveRotaded(curveWire,transfer);
			rotatedCurve =curveRotaded.Shape();
			curveWire =TopoDS::Wire(rotatedCurve);
			curveSections.AddWire(curveWire);


			BRepBuilderAPI_Transform rotated(newWholeWire,transfer);
			rotatedShape=rotated.Shape();
			newWholeWire=TopoDS::Wire(rotatedShape);
			sections.AddWire(newWholeWire);

			scaledSections.AddWire(scaledNewWire);
		}





		ang=360;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated1(initialwholeWire,transfer);
		rotatedShape=rotated1.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);

		scaledSections.AddWire(scaledNewWire);
		double newWholArea=getSurfaceArea(newWholeWire);

		CString initialWholeAreaString;
		CString expectedWholeAreaString;
		CString newWholeAreaString;
		CString msg;
		double exitAreaFraction=0.8;
		double exitPipeLength=40;

		msg=initialWholeAreaString+expectedWholeAreaString+newWholeAreaString;
		AfxMessageBox(msg);

		sections.Build();
		scaledSections.Build();

		TopoDS_Shape sh1=sections.Shape();
		TopoDS_Shape completeVolute=sections.Shape();
		//TopoDS_Shape scaledVolute=scaledSections.Shape();
		//circle
		gp_Dir dir(0,0,1); 
		double exitPipeRadius=sqrt((initailWholeArea*exitAreaFraction)/PI);


		gp_Pnt point(width/2,height*0.7,exitPipeLength);


		gp_Circ circle1(gp_Ax2( point, dir), exitPipeRadius);
		BRepBuilderAPI_MakeEdge circle(circle1);
		BRepBuilderAPI_MakeWire exitPipeCircleWire(circle);
		BRepBuilderAPI_MakeFace faceToBeRemoved(exitPipeCircleWire); 
		BRepOffsetAPI_ThruSections exitSection;		
		gp_Circ circleT(gp_Ax2( point, dir), exitPipeRadius+thickness*10);
		BRepBuilderAPI_MakeEdge circleT1(circleT);
		BRepBuilderAPI_MakeWire exitPipeCircleWireT(circleT1);
		BRepBuilderAPI_MakeFace exitFace(exitPipeCircleWireT);
		BRepOffsetAPI_ThruSections thicknessExitSection;


		thicknessExitSection.AddWire(scaledNewWire);
		thicknessExitSection.AddWire(exitPipeCircleWireT);
		thicknessExitSection.Build();
		TopoDS_Shape thicknessAddedPipe=thicknessExitSection.Shape();

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

		GeomAbs_JoinType joinType = GeomAbs_Arc;
		curveSections.Build();
		TopoDS_Shape sh4=curveSections.Shape();
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, curveSections.Shape());

		BRepBuilderAPI_Sewing connectedFaces(1.0e-06,true,true,true,true);

		for(TopExp_Explorer aFaceExplorer(sh4, TopAbs_FACE);aFaceExplorer.More();aFaceExplorer.Next()){
			TopoDS_Face aFace =TopoDS::Face(aFaceExplorer.Current());
			curveFace=aFace;


		}

		TopoDS_Face faces [20];
		int i=0;
		for(TopExp_Explorer aFaceExploree(voluteAndExit,TopAbs_FACE);aFaceExploree.More();aFaceExploree.Next()){

			TopoDS_Face aFace= TopoDS::Face(aFaceExploree.Current());
			faces[i]=aFace;
			i++;

		}


		for(int j=0;j<10;j++){
			BRepOffsetAPI_MakeOffsetShape curveFaceOffset(faces[j],3,1e-6,BRepOffset_Skin,false,false,joinType);
			connectedFaces.Add(curveFaceOffset);

			//	BRepBuilderAPI_Sewing 
		}

		BRepPrimAPI_MakeCylinder cylinder(12,5);
		BRepOffsetAPI_MakeOffsetShape voluteOffSet(completeVolute,3,1e-6,BRepOffset_Skin,false,false,joinType);
		//BRepOffsetAPI_MakeOffsetShape exitOffSet(,3,1e-6,BRepOffset_Skin,false,false,joinType);

		connectedFaces.Perform();
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN, completeVolute);
		m_pcoloredshapeList->Add(Quantity_NOC_RED, voluteOffSet);

		//->Add(Quantity_NOC_CORAL2,exitOffSet);
		TopoDS_Shape shVolute=connectedFaces.SewedShape();
		BRepTools::Write(shVolute,"D:/Breps/volute.brep");
		TopoDS_Compound completeThicknessVolute;
		BRep_Builder aBuilder1;
		aBuilder1.MakeCompound (completeThicknessVolute);
		aBuilder1.Add(completeThicknessVolute,exitPipe);

		BRepTools::Write(completeThicknessVolute,"D:/Breps/voluteThickness.brep");

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

	}


	void CImportExportDoc::OnThickness()
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
		double thickness=1;


		gp_Pnt r1(width/2.,height*2.,0);
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);
		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);
		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);


		BRepBuilderAPI_MakeEdge line1(p1,p2);
		BRepBuilderAPI_MakeEdge line2(q1,q2);
		BRepBuilderAPI_MakeEdge line3(q1,p1);
		BRepBuilderAPI_MakeEdge line4(q2,p2);



		BRepBuilderAPI_MakeWire voluteBaseWire1(line1,line3,line2);
		BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2,line4);

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

		TopoDS_Shape sh=trapeziumWire;
		Fit();

		double initialTrapaziumArea = getSurfaceArea(trapeziumWire);
		double initailWholeArea= getSurfaceArea(initialwholeWire);
		double expectedWholeArea = initailWholeArea*(degree_0_percentage/100);//Expected whole area

		double intialCurveArea=getSurfaceArea(initialCurveWire);
		double ExpectedTrapeziumArea = initialTrapaziumArea*(degree_0_percentage/100);//user percentage of the 0 point area
		double newTrapeziumHeight= getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);//calculating the height
		double ExpectedCurveArea = intialCurveArea*(degree_0_percentage/100);
		TopoDS_Wire newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,initialwholeWire);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,newWholeWire);

		TopoDS_Wire scaledNewWire;

		TopoDS_Shape rotatedShape;
		gp_Trsf transfer;
		BRepOffsetAPI_ThruSections sections;
		BRepOffsetAPI_ThruSections scaledSections;
		gp_Ax1 axis(gp_Pnt(0,-15,0),gp_Dir(1,0,0));

		scaledNewWire=createOuterShell(height,width,ang1Tan,ang2Tan,thickness);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,scaledNewWire);
		double ang=360;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated(newWholeWire,transfer);
		rotatedShape=rotated.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);


		scaledSections.AddWire(scaledNewWire);
		double areaIncreasingFactor=(100-degree_0_percentage)/18;
		CString areaIncreacingFactor;


		for(int i=2;i<=18;i++){

			double ang=20*(i-1);
			transfer.SetRotation(axis,ang*PI/180);
			double expectedWholeArea = initailWholeArea*((areaIncreasingFactor/100)*i);
			double ExpectedTrapeziumArea = initialTrapaziumArea*((areaIncreasingFactor/100)*i);
			double newTrapeziumHeight=getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);
			newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);
			m_pcoloredshapeList->Add(Quantity_NOC_RED,newWholeWire);
			scaledNewWire=createOuterShell(newTrapeziumHeight,width,ang1Tan,ang2Tan,thickness);
			m_pcoloredshapeList->Add(Quantity_NOC_GREEN,scaledNewWire);
			scaledSections.AddWire(scaledNewWire);
			sections.AddWire(newWholeWire);
			BRepBuilderAPI_Transform rotated(newWholeWire,transfer);

			rotatedShape=rotated.Shape();
			newWholeWire=TopoDS::Wire(rotatedShape);

		}


		ang=360;
		transfer.SetRotation(axis,ang*PI/180);
		BRepBuilderAPI_Transform rotated1(initialwholeWire,transfer);
		rotatedShape=rotated1.Shape();
		newWholeWire=TopoDS::Wire(rotatedShape);
		sections.AddWire(newWholeWire);
		scaledSections.AddWire(scaledNewWire);
		double newWholArea=getSurfaceArea(newWholeWire);

		CString initialWholeAreaString;
		CString expectedWholeAreaString;
		CString newWholeAreaString;
		CString msg;
		double exitAreaFraction=0.8;
		double exitPipeLength=40;

		msg=initialWholeAreaString+expectedWholeAreaString+newWholeAreaString;
		AfxMessageBox(msg);

		sections.Build();
		scaledSections.Build();

		TopoDS_Shape sh1=sections.Shape();
		TopoDS_Shape completeVolute=sections.Shape();

		gp_Dir dir(0,0,1); 
		double exitPipeRadius=sqrt((initailWholeArea*exitAreaFraction)/PI);


		gp_Pnt point(width/2,height*0.7,exitPipeLength);


		gp_Circ circle1(gp_Ax2( point, dir), exitPipeRadius);
		BRepBuilderAPI_MakeEdge circle(circle1);
		BRepBuilderAPI_MakeWire exitPipeCircleWire(circle);
		BRepBuilderAPI_MakeFace faceToBeRemoved(exitPipeCircleWire); 
		BRepOffsetAPI_ThruSections exitSection;		
		gp_Circ circleT(gp_Ax2( point, dir), exitPipeRadius+thickness*10);
		BRepBuilderAPI_MakeEdge circleT1(circleT);
		BRepBuilderAPI_MakeWire exitPipeCircleWireT(circleT1);
		BRepBuilderAPI_MakeFace exitFace(exitPipeCircleWireT);
		BRepOffsetAPI_ThruSections thicknessExitSection;


		thicknessExitSection.AddWire(scaledNewWire);
		thicknessExitSection.AddWire(exitPipeCircleWireT);
		thicknessExitSection.Build();
		TopoDS_Shape thicknessAddedPipe=thicknessExitSection.Shape();

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

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

	}



	gp_Pnt CImportExportDoc::getCentrePoint(TopoDS_Wire wire)
	{
		BRepBuilderAPI_MakeFace wireFace(wire);
		TopoDS_Face face =wireFace;
		BRepGProp gprop; 
		GProp_GProps surfaceProps;

		gprop.SurfaceProperties(wireFace,surfaceProps);

		gp_Pnt point=surfaceProps.CentreOfMass();

		return point;



	}

	TopoDS_Wire CImportExportDoc:: createNewShapeWithRightArea(double height,double width,double ang1Tan,double ang2Tan,double expectedArea,double initialHieght){


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


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);


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
		double Tolerance =  0.00001*expectedArea; 
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


	TopoDS_Wire CImportExportDoc::createOuterShell(double height, double width, double ang1Tan,double ang2Tan,double thickness)
	{
		gp_Pnt r1(width/2.,height*2+thickness,0);

		gp_Pnt q1(0-thickness,0-thickness,0); 
		gp_Pnt p1(width+thickness,0-thickness,0);

		gp_Pnt q2(-(height/ang2Tan)-thickness*2,height+thickness,0);
		gp_Pnt p2(((height/ang1Tan)+(width+thickness*2)),height+thickness,0);

		gp_Pnt q4(-((height-thickness+height/4)/ang2Tan),(height+height/4)+thickness,0);
		gp_Pnt p4(width+thickness+((height+height/4)/ang2Tan),(height+height/4)+thickness,0);

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



	TopoDS_Wire CImportExportDoc::createNewShapeAccordingToR1height(double height,double width,double ang1Tan,double ang2Tan,double r1Height){



		gp_Pnt r1(width/2.,r1Height,0);
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);
		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);
		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);



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



		/*TColgp_Array1OfPnt qCurvePoles(1,7); 
		qCurvePoles(1) = q2; 
		qCurvePoles(2) = q4; 
		qCurvePoles(3) = r1;
		qCurvePoles(4) = p4; 
		qCurvePoles(5) = p2; */

		return wholeWire;


	}

	/*TopoDS_Wire CImportExportDoc::createNewTrapazium(double height,double width, double ang1Tan,double ang2Tan,double expectedArea,double initialHeight){

	gp_Pnt q1(0,0,0); 
	gp_Pnt q2(-(height/ang2Tan),height,0);

	gp_Pnt p1(width,0,0);
	gp_Pnt p2(((height/ang1Tan)+width),height,0);

	BRepBuilderAPI_MakeEdge line1(p1,p2);
	BRepBuilderAPI_MakeEdge line2(q1,q2);
	BRepBuilderAPI_MakeEdge line3(q1,p1);
	BRepBuilderAPI_MakeEdge line4(q2,p2);
	BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2,line4);

	return trapeziumWire;


	}

	*/
	TopoDS_Wire CImportExportDoc:: createNewTrapazium(double height,double width, double ang1Tan,double ang2Tan)
	{
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);
		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);

		BRepBuilderAPI_MakeEdge line1(p1,p2);
		BRepBuilderAPI_MakeEdge line2(q1,q2);
		BRepBuilderAPI_MakeEdge line3(q1,p1);
		BRepBuilderAPI_MakeEdge line4(q2,p2);

		BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2);

		return trapeziumWire;


	} 

	TopoDS_Wire CImportExportDoc::createTrapazium(double height,double width,double ang1Tan,double ang2Tan){

		gp_Pnt r1(width/2.,height*2.,0);

		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);

		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);


		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);

		BRepBuilderAPI_MakeEdge line1(p1,p2);
		BRepBuilderAPI_MakeEdge line2(q1,q2);
		BRepBuilderAPI_MakeEdge line3(q1,p1);
		BRepBuilderAPI_MakeEdge line4(q2,p2);

		BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2);

		TColgp_Array1OfPnt qCurvePoles(1,5); 
		qCurvePoles(1) = q2; 
		qCurvePoles(2) = q4; 
		qCurvePoles(3) = r1;
		qCurvePoles(4) = p4; 
		qCurvePoles(5) = p2; 

		Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
		TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 

		BRepBuilderAPI_MakeWire wholeWire ( trapeziumWire,curveEdge);

		return wholeWire;

	}




	TopoDS_Wire CImportExportDoc::createCurveEdge(double height,double width,double ang1Tan,double ang2Tan,double expectedArea){

		double variationTolerance=expectedArea*0.01;
		double heightVariation=height*0.5;
		double r1Height;
		gp_Pnt r1(width/2.,height*2.,0);
		r1Height=height*2;
		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);

		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);

		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);

		BRepBuilderAPI_MakeEdge line4(q2,p2);

		TColgp_Array1OfPnt qCurvePoles(1,5); 

		qCurvePoles(1) = q2; 
		qCurvePoles(2) = q4; 
		qCurvePoles(3) = r1;
		qCurvePoles(4) = p4; 
		qCurvePoles(5) = p2; 
		Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
		TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 
		TopoDS_Wire curveWire= BRepBuilderAPI_MakeWire(curveEdge);
		return curveWire;




	}




	/*while(left<surfaceArea<right){

	if(surfaceArea<expectedArea){

	changedHeight=height+heightVariation;
	newCurveWire=curveConstructor(height,width,ang1Tan,ang2Tan,changedHeight);
	surfaceArea=getSurfaceArea(newCurveWire);

	}
	if(surfaceArea>expectedArea){

	changedHeight=height+heightVariation;
	newCurveWire=curveConstructor(height,width,ang1Tan,ang2Tan,changedHeight);
	surfaceArea=getSurfaceArea(newCurveWire);

	}*/







	TopoDS_Wire CImportExportDoc::curveConstructor(double height,double width,double ang1Tan,double ang2Tan,double r1Ycordinate){





		gp_Pnt r1(width/2.,r1Ycordinate,0);

		gp_Pnt q1(0,0,0); 
		gp_Pnt q2(-(height/ang2Tan),height,0);

		gp_Pnt p1(width,0,0);
		gp_Pnt p2(((height/ang1Tan)+width),height,0);


		gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);

		BRepBuilderAPI_MakeEdge line4(q2,p2);

		TColgp_Array1OfPnt qCurvePoles(1,5); 

		qCurvePoles(1) = q2; 

		qCurvePoles(2) = q4; 

		qCurvePoles(3) = r1;

		qCurvePoles(4) = p4; 

		qCurvePoles(5) = p2; 


		Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
		TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 

		BRepBuilderAPI_MakeWire CurveWire(curveEdge,line4);
		return CurveWire;

	}






	double CImportExportDoc::numberTesting(){

		double  number =5;

		CString fivs;

		fivs.Format(_T("root %g \n"),number);
		AfxMessageBox(fivs);

		return number;

	}


	void CImportExportDoc::OnCut(){
		//

		double bx,by,bz;

		bx=5;
		by=5;
		bz=5;
		BRepPrimAPI_MakeBox box1(by,by,bz);





		////gp_Ax2 cyAx2(cyLocation, cyAxis);
		//
		//
		//BRepPrimAPI_MakeCylinder cylinder(cyAx2,1,5);
		TopoDS_Shape s1, s2,s3;
		s1 = box1.Shape();
		//s2 = cylinder.Shape();
		//
		//BRepAlgoAPI_Cut drilled(s1,s2);	





		gp_Pnt cyLocation(5/2, 5/2, 0);
		gp_Dir cyAxis = gp::DZ();




		//second method

		Standard_Real feature_diameter = 0.2;
		gp_Ax1 feature_origin = gp_Ax1(gp_Pnt(bx/2,0,bz/2), gp_Dir(0,1,0));

		BRepFeat_MakeCylindricalHole feature_maker;
		feature_maker.Init(s1,feature_origin );

		feature_maker.Perform(feature_diameter);

		s3=feature_maker.Shape();



		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, s3);
		m_pcoloredshapeList->Display(myAISContext);

		Fit();

	}



	void CImportExportDoc::OnMakeBoxDrill(){


		double length;
		double width;
		double height;

		length =10;
		width =10;
		height =10;


		gp_Pnt pnt1(length,width,0);
		gp_Pnt pnt2(-length,width,0);
		gp_Pnt pnt3(-length,-width,0);
		gp_Pnt pnt4(length,-width,0);


		Handle(Geom_TrimmedCurve) segment1= GC_MakeSegment(pnt1,pnt2);
		Handle(Geom_TrimmedCurve) segment2= GC_MakeSegment(pnt2,pnt3);
		Handle(Geom_TrimmedCurve) segment3= GC_MakeSegment(pnt3,pnt4);
		Handle(Geom_TrimmedCurve) segment4= GC_MakeSegment(pnt4,pnt1);

		TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(segment1);
		TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(segment2);
		TopoDS_Edge aEdge3 = BRepBuilderAPI_MakeEdge(segment3);
		TopoDS_Edge aEdge4 = BRepBuilderAPI_MakeEdge(segment4);


		TopoDS_Wire wire1 = BRepBuilderAPI_MakeWire(aEdge1, aEdge2, aEdge3,aEdge4);

		gp_Ax1 xAxis = gp::OX();

		TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(wire1);

		gp_Vec aPrismVec(0, 0, height);

		TopoDS_Shape myBody = BRepPrimAPI_MakePrism(myFaceProfile, aPrismVec);



		m_pcoloredshapeList->Add(Quantity_NOC_STEELBLUE1,myBody);

		m_pcoloredshapeList->Display(myAISContext);
		Fit();





	}



	void CImportExportDoc::OnBREPFile(){

		BRepBuilderAPI_Sewing sewer(1.0e-7);
		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadBREPNew();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			sewer.Add(aSeqOfShape->Value(i));


		}
		sewer.Perform();
		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, sewer.SewedShape());
		m_pcoloredshapeList->Display(myAISContext);
		Fit();



	}

	void CImportExportDoc::OnFileImportIges() 
	{   
		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadIGES();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, aSeqOfShape->Value(i));
			m_pcoloredshapeList->Display(myAISContext);
		}
		Fit();
	}
	void CImportExportDoc::OnFileExportIges() 
	{   CImportExport::SaveIGES(myAISContext);}
	//-------------------------------------------------------------------------


	void CImportExportDoc::OnFileImportStep() 
	{   
		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadSTEP();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, aSeqOfShape->Value(i));
			m_pcoloredshapeList->Display(myAISContext);
		}
		Fit();
	}


	//--------------------------------------------------------------------------


	void CImportExportDoc::OnSTEPFile() 
	{   
		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadSTEPNew();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, aSeqOfShape->Value(i));
			m_pcoloredshapeList->Display(myAISContext);
		}
		Fit();
	}

	void CImportExportDoc::OnFileExportStep() 
	{   CImportExport::SaveSTEP(myAISContext);}


	void CImportExportDoc::OnFileExportVrml() 
	{   CImportExport::SaveVRML(myAISContext);}

	void CImportExportDoc::OnFileExportStl() 
	{   CImportExport::SaveSTL(myAISContext);}

	void  CImportExportDoc::Popup(const Standard_Integer  x,
		const Standard_Integer  y ,
		const Handle(V3d_View)& aView   ) 
	{
		Standard_Integer PopupMenuNumber=0;
		myAISContext->InitCurrent();
		if (myAISContext->MoreCurrent())
			PopupMenuNumber=1;

		CMenu menu;
		VERIFY(menu.LoadMenu(IDR_Popup3D));
		CMenu* pPopup = menu.GetSubMenu(PopupMenuNumber);

		ASSERT(pPopup != NULL);
		if (PopupMenuNumber == 1) // more than 1 object.
		{
			bool OneOrMoreInShading = false;
			for (myAISContext->InitCurrent();myAISContext->MoreCurrent ();myAISContext->NextCurrent ())
				if (myAISContext->IsDisplayed(myAISContext->Current(),1)) OneOrMoreInShading=true;
			if(!OneOrMoreInShading)
				pPopup->EnableMenuItem(5, MF_BYPOSITION | MF_DISABLED | MF_GRAYED);
		}

		POINT winCoord = { x , y };
		Handle(WNT_Window) aWNTWindow=
			Handle(WNT_Window)::DownCast(aView->Window());
		ClientToScreen ( (HWND)(aWNTWindow->HWindow()),&winCoord);
		pPopup->TrackPopupMenu(TPM_LEFTALIGN | TPM_RIGHTBUTTON , winCoord.x, winCoord.y , 
			AfxGetMainWnd());


	}



	void CImportExportDoc::OnBox() 
	{
		AIS_ListOfInteractive aList;
		myAISContext->DisplayedObjects(aList);
		AIS_ListIteratorOfListOfInteractive aListIterator;
		for(aListIterator.Initialize(aList);aListIterator.More();aListIterator.Next()){
			myAISContext->Remove(aListIterator.Value());
		}

		BRepPrimAPI_MakeBox B(200.,150.,100.);

		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, B.Shape());

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

		// document has been modified
		SetModifiedFlag(TRUE);
	}


	void CImportExportDoc::OnBox1() 
	{
		AIS_ListOfInteractive aList;
		myAISContext->DisplayedObjects(aList);
		AIS_ListIteratorOfListOfInteractive aListIterator;
		for(aListIterator.Initialize(aList);aListIterator.More();aListIterator.Next()){
			myAISContext->Remove(aListIterator.Value());
		}



		// B->~BRepPrimAPI_MakeBox(200.,150.,100.);

		gp_Ax1 origin = gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1));

		TopoDS_Shape drillBox =box.Shape();


		Standard_Real feature_diameter = 20.;






		Fit();

		// document has been modified
		SetModifiedFlag(TRUE);
	}


	void CImportExportDoc::OnFillet(){


		AIS_ListOfInteractive aList;
		myAISContext->DisplayedObjects(aList);
		AIS_ListIteratorOfListOfInteractive aListIterator;
		for(aListIterator.Initialize(aList);aListIterator.More();aListIterator.Next()){
			myAISContext->Remove(aListIterator.Value());
		}
		//myAISContext->RemoveAll();
		m_pcoloredshapeList->Remove( box.Shape());



		BRepFilletAPI_MakeFillet mkFillet(box.Shape());

		TopoDS_Edge anEdge ;

		TopExp_Explorer anEdgeExplorer(box.Shape(), TopAbs_EDGE);

		while(anEdgeExplorer.More()){

			anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			mkFillet.Add(30., anEdge);
			anEdgeExplorer.Next();
		}



		TopoDS_Shape B1 = mkFillet.Shape();


		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, B1);

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

		// document has been modified
		SetModifiedFlag(TRUE);




	}



	void CImportExportDoc::OnFilletDialog(){





		TopoDS_Shape importedShape;

		Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadBREP();
		for(int i=1;i<= aSeqOfShape->Length();i++)
		{
			importedShape=aSeqOfShape->Value(i);

		}



		TopoDS_Wire wires[150];

		TopExp_Explorer anWireExplorer(importedShape, TopAbs_WIRE);
		int i=0;
		while(anWireExplorer.More()){
			//TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			TopoDS_Wire anEdge = TopoDS::Wire(anWireExplorer.Current());
			wires[i]=anEdge;
			//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,wires[i]);
			anWireExplorer.Next();
			i++;
		}

		TopoDS_Edge edges[4];
		TopExp_Explorer anEdgeExplorer(wires[0],TopAbs_EDGE);
		int j=0;
		while(anEdgeExplorer.More()){
			TopoDS_Edge anEdge=TopoDS::Edge(anEdgeExplorer.Current());
			edges[j]=anEdge;

			anEdgeExplorer.Next();
			j++;
		}

		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW,edges[0]);
		m_pcoloredshapeList->Add(Quantity_NOC_GREEN,edges[1]);

		double length=sizeof(wires);

		//CString str;
		//str.Format(_T("length g% :"),length);

		//AfxMessageBox(str);

		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, curve);
		//m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, Elps);
		m_pcoloredshapeList->Add(Quantity_NOC_RED,wires[0]);



		m_pcoloredshapeList->Display(myAISContext);
		Fit();


	}


	void CImportExportDoc::OnBottle() 
	{



		double myWidth=50;
		double myThickness=30;
		double myHeight=70;

		gp_Pnt aPnt1(-myWidth / 2., 0, 0);
		gp_Pnt aPnt2(-myWidth / 2., -myThickness / 4., 0);	
		gp_Pnt aPnt3(0, -myThickness / 2., 0);
		gp_Pnt aPnt4(myWidth / 2., -myThickness / 4., 0);
		gp_Pnt aPnt5(myWidth / 2., 0, 0);


		Standard_Real xValue1 = aPnt1.X();

		Handle(Geom_TrimmedCurve) aArcOfCircle = GC_MakeArcOfCircle(aPnt2,aPnt3,aPnt4);
		Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
		Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(aPnt4, aPnt5);


		TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(aSegment1);
		TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(aArcOfCircle);
		TopoDS_Edge aEdge3 = BRepBuilderAPI_MakeEdge(aSegment2);	
		TopoDS_Edge aEdge4 = BRepBuilderAPI_MakeEdge(aPnt1, aPnt3);
		TopoDS_Edge aEdge5 = BRepBuilderAPI_MakeEdge(aPnt4, aPnt5);


		TopoDS_Wire aWire = BRepBuilderAPI_MakeWire(aEdge1, aEdge2, aEdge3);

		gp_Pnt aOrigin(0, 0, 0);
		gp_Dir xDir(1, 0, 0);
		gp_Ax1 xAxis(aOrigin, xDir);

		xAxis = gp::OX();

		gp_Trsf aTrsf;
		aTrsf.SetMirror(xAxis);

		BRepBuilderAPI_Transform aBRepTrsf(aWire, aTrsf);
		TopoDS_Shape aMirroredShape = aBRepTrsf.Shape();
		TopoDS_Wire aMirroredWire = TopoDS::Wire(aMirroredShape);

		BRepBuilderAPI_MakeWire mkWire;
		mkWire.Add(aWire);
		mkWire.Add(aMirroredWire);
		TopoDS_Wire myWireProfile = mkWire.Wire();
		TopoDS_Shape sh  =  myWireProfile;
		TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(myWireProfile);
		gp_Vec aPrismVec(0, 0, myHeight);

		TopoDS_Shape myBody = BRepPrimAPI_MakePrism(myFaceProfile, aPrismVec);
		BRepFilletAPI_MakeFillet mkFillet(myBody);

		TopoDS_Edge anEdge ;

		TopExp_Explorer anEdgeExplorer(myBody, TopAbs_EDGE);

		while(anEdgeExplorer.More()){

			anEdge = TopoDS::Edge(anEdgeExplorer.Current());
			mkFillet.Add(myThickness / 12., anEdge);
			anEdgeExplorer.Next();
		}



		myBody = mkFillet.Shape();


		gp_Pnt neckLocation(0, 0, myHeight);
		gp_Dir neckAxis = gp::DZ();
		gp_Ax2 neckAx2(neckLocation, neckAxis);

		Standard_Real myNeckRadius = myThickness / 4.;
		Standard_Real myNeckHeight = myHeight /10 ;
		BRepPrimAPI_MakeCylinder MKCylinder(neckAx2, myNeckRadius, myNeckHeight);
		TopoDS_Shape myNeck = MKCylinder.Shape();

		BRepTools::Write(myNeck,"myNeck.brep");
		BRepTools::Write(myBody,"myBody.brep");

		myBody =BRepAlgo_Fuse(myBody, myNeck);

		TopoDS_Face faceToRemove;
		Standard_Real zMax = -1;




		for(TopExp_Explorer aFaceExplorer(myBody, TopAbs_FACE); aFaceExplorer.More(); aFaceExplorer.Next()){
			TopoDS_Face aFace = TopoDS::Face(aFaceExplorer.Current());
			// Check if <aFace> is the top face of the bottle's neck 
			Handle(Geom_Surface) aSurface = BRep_Tool::Surface(aFace);
			if(aSurface->DynamicType() == STANDARD_TYPE(Geom_Plane)){
				Handle(Geom_Plane) aPlane = Handle(Geom_Plane)::DownCast(aSurface);
				gp_Pnt aPnt = aPlane->Location();
				Standard_Real aZ   = aPnt.Z();
				if(aZ > zMax){
					zMax = aZ;
					faceToRemove = aFace;
				}
			}
		}

		TopTools_ListOfShape facesToRemove;
		facesToRemove.Append(faceToRemove);
		myBody = BRepOffsetAPI_MakeThickSolid(myBody, facesToRemove, -myThickness / 50, 1.e-3);

		Handle(Geom_CylindricalSurface) aCyl1 = new Geom_CylindricalSurface(neckAx2, myNeckRadius * 0.99);
		Handle(Geom_CylindricalSurface) aCyl2 = new Geom_CylindricalSurface(neckAx2, myNeckRadius * 1.05);

		gp_Pnt2d aPnt(2. * M_PI, myNeckHeight / 2.);
		gp_Dir2d aDir(2. * M_PI, myNeckHeight / 4.);
		gp_Ax2d anAx2d(aPnt, aDir);
		Standard_Real aMajor = 2. * M_PI;
		Standard_Real aMinor = myNeckHeight / 10;



		Handle(Geom2d_Ellipse) anEllipse1 = new Geom2d_Ellipse(anAx2d, aMajor, aMinor);
		Handle(Geom2d_Ellipse) anEllipse2 = new Geom2d_Ellipse(anAx2d, aMajor, aMinor / 4);
		Handle(Geom2d_TrimmedCurve) anArc1 = new Geom2d_TrimmedCurve(anEllipse1, 0, M_PI);
		Handle(Geom2d_TrimmedCurve) anArc2 = new Geom2d_TrimmedCurve(anEllipse2, 0, M_PI);
		gp_Pnt2d anEllipsePnt1 = anEllipse1->Value(0);
		gp_Pnt2d anEllipsePnt2 = anEllipse1->Value(M_PI);

		Handle(Geom2d_TrimmedCurve) aSegment = GCE2d_MakeSegment(anEllipsePnt1, anEllipsePnt2);

		TopoDS_Edge anEdge1OnSurf1 = BRepBuilderAPI_MakeEdge(anArc1, aCyl1);
		TopoDS_Edge anEdge2OnSurf1 = BRepBuilderAPI_MakeEdge(aSegment, aCyl1);
		TopoDS_Edge anEdge1OnSurf2 = BRepBuilderAPI_MakeEdge(anArc2, aCyl2);
		TopoDS_Edge anEdge2OnSurf2 = BRepBuilderAPI_MakeEdge(aSegment, aCyl2);
		TopoDS_Wire threadingWire1 = BRepBuilderAPI_MakeWire(anEdge1OnSurf1, anEdge2OnSurf1);
		TopoDS_Wire threadingWire2 = BRepBuilderAPI_MakeWire(anEdge1OnSurf2, anEdge2OnSurf2);
		BRepLib::BuildCurves3d(threadingWire1);
		BRepLib::BuildCurves3d(threadingWire2);

		BRepOffsetAPI_ThruSections aTool(Standard_True);
		aTool.AddWire(threadingWire1);
		aTool.AddWire(threadingWire2);
		aTool.CheckCompatibility(Standard_False);

		TopoDS_Shape myThreading = aTool.Shape();

		TopoDS_Compound aRes;
		BRep_Builder aBuilder;
		aBuilder.MakeCompound (aRes);
		aBuilder.Add (aRes, myBody);
		aBuilder.Add (aRes, myThreading);


		TopoDS_Shape All=aRes;





		//bubu



		m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, All);

		m_pcoloredshapeList->Display(myAISContext);
		Fit();


	}




	void CImportExportDoc::OnTest(){







	}


	void CImportExportDoc::OnIGESFile(){

		{  
			Handle(TopTools_HSequenceOfShape) aSeqOfShape = CImportExport::ReadIGESNew();
			for(int i=1;i<= aSeqOfShape->Length();i++)
			{
				m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, aSeqOfShape->Value(i));
				m_pcoloredshapeList->Display(myAISContext);
			}



			Fit();
		}
	}

	void CImportExportDoc::OnCylinder() 
	{
		AIS_ListOfInteractive aList;
		myAISContext->DisplayedObjects(aList);
		AIS_ListIteratorOfListOfInteractive aListIterator;
		for(aListIterator.Initialize(aList);aListIterator.More();aListIterator.Next()){
			myAISContext->Remove(aListIterator.Value());
		}

		BRepPrimAPI_MakeCylinder C(50.,200.);

		m_pcoloredshapeList->Add(Quantity_NOC_GREEN, C.Shape());

		m_pcoloredshapeList->Display(myAISContext);
		Fit();

		// document has been modified
		SetModifiedFlag(TRUE);
	}
	void CImportExportDoc::OnObjectRemove() 

	{
		for(GetAISContext()->InitCurrent();GetAISContext()->MoreCurrent();GetAISContext()->NextCurrent()) {
			Handle(AIS_Shape) aShape = Handle(AIS_Shape)::DownCast(GetAISContext()->Current());
			if(!aShape.IsNull()) {
				m_pcoloredshapeList->Remove(aShape->Shape());
			}
		}
		OCC_3dBaseDoc::OnObjectRemove();
	}

	void CImportExportDoc::OnObjectErase() 

	{
		for(GetAISContext()->InitCurrent();GetAISContext()->MoreCurrent();GetAISContext()->NextCurrent()) {
			Handle(AIS_Shape) aShape = Handle(AIS_Shape)::DownCast(GetAISContext()->Current());
			if(!aShape.IsNull()) {
				m_pcoloredshapeList->Remove(aShape->Shape());
			}
		}
		OCC_3dBaseDoc::OnObjectErase(); 
	}

	void CImportExportDoc::OnObjectDisplayall() 

	{
		OCC_3dBaseDoc::OnObjectDisplayall(); 
	}