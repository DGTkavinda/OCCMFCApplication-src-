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


	double CImportExportDoc::getSurfaceArea(TopoDS_Wire wireShape){

		BRepBuilderAPI_MakeFace wireFace(wireShape);
		TopoDS_Face face =wireFace;
		BRepGProp gprop; 
		GProp_GProps surfaceProps;



		gprop.SurfaceProperties(wireFace,surfaceProps);

		double surfaceArea=surfaceProps.Mass();
		CString fivs;
		fivs.Format(_T("surface %g\n",surfaceArea));
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

		//double angle1=45.00;
		//double angle2=45.00;
		//double height=25+5;
		//double initialHeight=height;
		//double width=10;
		//double startPointx=0;
		//double ang1Tan= tan(angle1*PI/180.0);
		//double ang2Tan= tan(angle2*PI/180.0);
		//double area=150;
		//double degree_0_percentage=5;


		//gp_Pnt r1(width/2.,height*2.,0);
		//gp_Pnt q1(0,0,0); 
		//gp_Pnt q2(-(height/ang2Tan),height,0);
		//gp_Pnt p1(width,0,0);
		//gp_Pnt p2(((height/ang1Tan)+width),height,0);
		//gp_Pnt q3(-((height+height/10)/ang2Tan),(height+height/10),0);
		//gp_Pnt q4(-((height+height/4)/ang2Tan),(height+height/4),0);
		//gp_Pnt q5(-(height/ang2Tan),(height+height/3),0);


		//gp_Pnt p3(width+((height+height/10)/ang2Tan),(height+height/10),0);
		//gp_Pnt p4(width+((height+height/4)/ang2Tan),(height+height/4),0);
		//gp_Pnt p5((width+height/ang2Tan),(height+height/3),0);


		//BRepBuilderAPI_MakeEdge line1(p1,p2);
		//BRepBuilderAPI_MakeEdge line2(q1,q2);
		//BRepBuilderAPI_MakeEdge line3(q1,p1);
		//BRepBuilderAPI_MakeEdge line4(q2,p2);



		//BRepBuilderAPI_MakeWire voluteBaseWire1(line1,line3,line2);
		//BRepBuilderAPI_MakeWire trapeziumWire(line1,line3,line2,line4);

		////----------------------------------

		//TColgp_Array1OfPnt qCurvePoles(1,5); 
		//qCurvePoles(1) = q2; 
		//qCurvePoles(2) = q4; 
		//qCurvePoles(3) = r1;
		//qCurvePoles(4) = p4; 
		//qCurvePoles(5) = p2; 
		//Handle(Geom_BezierCurve) qCurve = new Geom_BezierCurve(qCurvePoles); 
		//TopoDS_Edge curveEdge = BRepBuilderAPI_MakeEdge(qCurve); 
		//BRepBuilderAPI_MakeWire initialCurveWire(curveEdge,line4);
		//BRepBuilderAPI_MakeWire initialwholeWire (voluteBaseWire1,curveEdge);

		//BRepTools::Write( initialwholeWire,"D:/Breps/tRapAxis.brep");
		////BRepBuilderAPI_MakeWire wire2(wire1,mainE);
		////BRepBuilderAPI_MakeWire voluteCross(wire2,pE)
		//TopoDS_Shape sh=trapeziumWire;
		//Fit();

		////for 50% piece







		//double initialTrapaziumArea = getSurfaceArea(trapeziumWire);
		//double initailWholeArea= getSurfaceArea(initialwholeWire);
		//double expectedWholeArea = initailWholeArea*(degree_0_percentage/100);//Expected whole area

		//double intialCurveArea=getSurfaceArea(initialCurveWire);
		//double ExpectedTrapeziumArea = initialTrapaziumArea*(degree_0_percentage/100);//user percentage of the 0 point area
		//double newTrapeziumHeight= getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);//calculating the height
		//double ExpectedCurveArea = intialCurveArea*(degree_0_percentage/100);
		//TopoDS_Wire newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);



		//TopoDS_Shape rotatedShape;
		//gp_Trsf transfer;
		//BRepOffsetAPI_ThruSections sections;

		//gp_Ax1 axis(gp_Pnt(0,-15,0),gp_Dir(1,0,0));
		////initial wire rotation
		//

		//double ang=360;
		//transfer.SetRotation(axis,ang*PI/180);
		//BRepBuilderAPI_Transform rotated(newWholeWire,transfer);
		//rotatedShape=rotated.Shape();
		//newWholeWire=TopoDS::Wire(rotatedShape);
		//sections.AddWire(newWholeWire);


		//double areaIncreasingFactor=(100-degree_0_percentage)/18;

		//CString areaIncreacingFactor;





		//for(int i=2;i<=18;i++){

		//	double ang=20*(i-1);
		//	transfer.SetRotation(axis,ang*PI/180);


		//	double expectedWholeArea = initailWholeArea*((areaIncreasingFactor/100)*i);
		//	double ExpectedTrapeziumArea = initialTrapaziumArea*((areaIncreasingFactor/100)*i);
		//	double newTrapeziumHeight=getTrapezuimHeight(ExpectedTrapeziumArea,width,ang1Tan,ang2Tan);
		//	newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);


		//	BRepBuilderAPI_Transform rotated(newWholeWire,transfer);

		//	rotatedShape=rotated.Shape();
		//	newWholeWire=TopoDS::Wire(rotatedShape);


		//	sections.AddWire(newWholeWire);


		//}


		//ang=360;
		//transfer.SetRotation(axis,ang*PI/180);
		//BRepBuilderAPI_Transform rotated1(initialwholeWire,transfer);
		//rotatedShape=rotated1.Shape();
		//newWholeWire=TopoDS::Wire(rotatedShape);
		//sections.AddWire(newWholeWire);




		////TopoDS_Wire newWholeWire=createNewShapeWithRightArea(newTrapeziumHeight,width,ang1Tan,ang2Tan,expectedWholeArea,initialHeight);

		//double newWholArea=getSurfaceArea(newWholeWire);



		////TopoDS_Wire newCurveWire=createNewCurve(newTrapeziumHeight,width,ang1Tan,ang2Tan,ExpectedCurveArea);
		////double newlyCreatedCurveWireArea=getSurfaceArea(newCurveWire);
		///*TopoDS_Wire newWhaoleshape=createTrapazium(newTrapeziumHeight,width,ang1Tan,ang2Tan);
		//double areaOfNewWholeShape=getSurfaceArea(newWhaoleshape);*/
		////TopoDS_Wire newTrapazium=createNewTrapazium(newTrapeziumHeight,width,ang1Tan,ang2Tan);
		////double newTrapazimArea=getSurfaceArea(newTrapazium);

		////sections.Build();

		////double number=getTrapezuimHeight(initialTrapaziumArea,width,ang1Tan,ang2Tan);




		//CString initialWholeAreaString;
		//CString expectedWholeAreaString;
		//CString newWholeAreaString;
		//CString msg;
		//initialWholeAreaString.Format(_T("initialWholeAreaString %g \n" ),initailWholeArea);
		//expectedWholeAreaString.Format(_T("expectedWholeAreaString %g \n"),expectedWholeArea);
		//newWholeAreaString.Format(_T("newWholeArea %g \n"),newWholArea);





		////areaOfNewTrapziumString.Format(_T("area of new Trapazium %g \n"),newTrapazimArea);
		//msg=initialWholeAreaString+expectedWholeAreaString+newWholeAreaString;
		//AfxMessageBox(msg);

		//sections.Build();

		//TopoDS_Shape sh1=sections.Shape();

		//BRepTools::Write(sh1,"D:/Breps/whoaleShape.brep");	


		////m_pcoloredshapeList->Add(Quantity_NOC_YELLOW, sh1);
		//m_pcoloredshapeList->Add(Quantity_NOC_STEELBLUE,sh1);
		//m_pcoloredshapeList->Display(myAISContext);




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

		double Tolerance =  0.01*expectedArea; 


		double topMiddle=up;
		double bottomMiddle=bottom; 

		if(fabs(expectedArea-surfaceArea)<Tolerance){
			/*CString  msg="if condition";
			AfxMessageBox(msg);*/
			return wholeWire;

		}else{


			/*CString  msg="else";
			AfxMessageBox(msg);*/

			while(fabs(expectedArea-surfaceArea)>Tolerance){

				//CString  msg="while";
				//AfxMessageBox(msg);
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




	TopoDS_Edge CImportExportDoc::createCurveEdge(double height,double width,double ang1Tan,double ang2Tan,double expectedArea){

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
		return curveEdge;




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


		//filletDialog = new CFilletDialog();
		std::string s ="Test";

		// filletDialog->Create( IDD_DIALOG2);



		// filletDialog->ShowWindow(SW_SHOW);


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