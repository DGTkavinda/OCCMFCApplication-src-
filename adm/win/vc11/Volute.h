#pragma once
#include <vector>
#include "Geom_BSplineCurve.hxx"
#include "Geom_BezierCurve.hxx"
#include <ColoredShapes.h>
#include <OCC_3dDoc.h>
#include "BRepPrimAPI_MakeBox.hxx"
#include "Geom_BezierCurve.hxx"
#include <BRepGProp.hxx>
#include "GProp_GProps.hxx"
#include <BRepOffsetAPI_ThruSections.hxx>
#include "gp_Trsf.hxx"
#include "VoluteDialog.h"
#include "GeomAbs_JoinType.hxx"
#include "BRepOffsetAPI_MakeOffsetShape.hxx"
#include "BRepOffset_Mode.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_MakePipeShell.hxx"

class Geom_BSplineSurface;
class TopoDS_Wire;
//class Handle_Geom_BSplineCurve;
class Geom_Ellipse;
//class Handle_Geom_BezierCurve;
class TopoDS_Shape;
class TopoDS_Solid;
class TopoDS_Face;
class TopoDS_Wire;
class TopoDS_Edge;

__declspec(dllexport) class OCCVolute
{
public:
    __declspec(dllexport) void makePumpVolute(bool isClockWise,double hub_inlet_z,double shroud_inlet_z, double hub_inlet_radius,
                         double shroud_inlet_radius, double r_offset, double inlet_area, double throat_area,
                         double alpha1,double alpha2,double topRatio,
                         //wrap region data
                         double wrapLength, double m_TailThickness, double m_le_aspect_ratio,
                         //exit pipe data
                         double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v,
                         int exit_curve_order, double* pX_exitCrv, 
                         double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
                         double exit_pipe_aspect1, double exit_pipe_aspect2, double exit_pipe_lean1,
                         double exit_pipe_lean2);

	__declspec(dllexport) TopoDS_Shape getScrollShape();
	__declspec(dllexport) TopoDS_Shape getTongueShape();
	__declspec(dllexport) TopoDS_Shape getVoluteTotal();
	__declspec(dllexport) TopoDS_Shape getVoluteTotalThick();

    __declspec(dllexport) void getScrollSurf(Handle_Geom_BSplineSurface& hBsp);
    __declspec(dllexport) void gettongueSurf(Handle_Geom_BSplineSurface& hBsp);
    __declspec(dllexport) void getwrapRegionSurf(Handle_Geom_BSplineSurface& hBsp);
    __declspec(dllexport) void getexitPipeSurf0(Handle_Geom_BSplineSurface& hBsp);
    __declspec(dllexport) void getexitPipeSurf1(Handle_Geom_BSplineSurface& hBsp);
    __declspec(dllexport) void clearVolute();
    __declspec(dllexport) double getExitPipe0SectionArea(Standard_Real Vsect);
    __declspec(dllexport) double getWrapSectionArea(Standard_Real Vsect);
    __declspec(dllexport) void getTestCurve(Handle_Geom_Curve& Crv);



	__declspec(dllexport) struct asymSection
	{
		int npt;
		std::vector<double> x, y, z;
		std::vector<int> landmark;
		double thetaSect;
	};
	__declspec(dllexport) int makePumpVolute2(bool doBoolean, bool doFillets, std::vector<asymSection>& asymSects, double angTran, bool IsClockwise,
		bool isTurbine,double throat_area,
		bool firstRun, bool overRideTransPos, int& runFlag, double& calcedSecondTransPos, double earlyTurnAngle, double alpha1,double alpha2,
		//exit pipe data
         double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h,
         int exit_curve_order, double* pX_exitCrv, 
         double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
         double exit_pipe_aspect1, double exit_pipe_aspect2, int firstTransPosOption,double firstTransPos, int transSectOption,
		 bool exit_enable_extension, double exit_extension_length, double exit_extension_reduction,
		 //tongue data
		 int filletMethod,double tongueCurvControlParam,double tongueSizeControlParam,int filletOption, double filletRmin, double filletRmax,
		 double tongue_le_radius, double tongue_aspect_ratio, double tongue_z_aspect_ratio,bool Istongueless, double TongueStAngle,
		 //splitter data
		 bool isSplitter, double splitterStartAng, double splitterEndDist, std::vector<double>& thicknessAng, std::vector<double>& thicknessVal,
		 bool isConstAreaMV, int areaOptionPV, std::vector<double>& areaAng, std::vector<double>& areaRatioPV, std::vector<double>& areaAbsPV);


	__declspec(dllexport) int makeThickVolute(double wallThickness, bool isTurbine, bool exit_straight, bool isClockWise);

protected:
    Handle_Geom_BSplineSurface pScrollSurf;
    Handle_Geom_BSplineSurface pTongueSurf;
    Handle_Geom_BSplineSurface pTipBottomSurf;
    Handle_Geom_BSplineSurface pTipTopSurf;
    Handle_Geom_BSplineSurface pwrapRegionSurf;
    Handle_Geom_BSplineSurface pexitPipeSurf0;
    Handle_Geom_BSplineSurface pexitPipeSurf1;

    gp_Dir dirWrapEnd;
    gp_Dir dirTransEnd;
    Handle_Geom_BezierCurve pTongueStartCurve;
    Handle_Geom_BSplineCurve pTongueEndCurve;
    Handle_Geom_BSplineCurve pTongueSpineCurve;
    Standard_Real swirlSectU1;
    Standard_Real swirlSectU2;
	Handle_Geom_BezierCurve pexitPipeSpine_Seg1;
	Handle_Geom_BezierCurve pexitPipeSpine_Seg2;


    Handle_Geom_BSplineSurface copySurfaceToNew(Handle_Geom_BSplineSurface& pBSurf);
    Handle_Geom_BSplineCurve copyCurveToNew(Handle_Geom_BSplineCurve& pCurve);
    Handle_Geom_BSplineCurve makeSingleCurve(TopoDS_Wire wire, int numPts = 100);
    gp_Pnt findPointofTheta(Handle_Geom_Curve& pCurve, double theta, double& U);
    gp_Pnt findPointofTheta(Handle_Geom_BSplineSurface& pSurf, double U,double theta, double& V);
    gp_Pnt findPointofRadius(Handle_Geom_Curve& pCurve, double radius, bool fromEnd, double& U);


    bool makeWrapRegionSurf(double wrapLength, double m_TailThickness);

    bool makeExitTransitionSurf(bool isClockWise,double throat_area, double exit_area, double exit_length, bool exit_straight,
        double exit_pipe_angle_v, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv,
        bool exit_enable_transition, double transitionLength, double exit_pipe_aspect1, double exit_pipe_lean1);

    bool makeExitPipeSurf(bool isClockWise,double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v,
                         int exit_curve_order, double* pX_exitCrv, 
                         double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
                         double exit_pipe_aspect2, double exit_pipe_lean2);

    //returns an ellipse major/minor radius for a given target area
    bool getTargetEllipse(double targetArea, double aspectRatio, double& R_major, double& R_minor);
    //returns a curve from the ellipse with U=0 set at the minimum radius point
    Handle_Geom_BSplineCurve resetEllipsetoCurve(Geom_Ellipse& E, double fraction);
    //returns the wrap section shape at 360 degrees
    TopoDS_Wire GetWrapRegionBaseSection();

    Handle_Geom_Curve testCurve;

    //new functions and variables
	TopoDS_Shape cutScroll;	
    TopoDS_Shape tongueShape;
	TopoDS_Shape tonguePrism;   
	TopoDS_Shape tongueFinal;
    TopoDS_Shape tongueFinal1;
    TopoDS_Shape tongueFinal2;
	TopoDS_Shape wrapRegionShape;
	TopoDS_Shape exitPipe_Seg1;
	TopoDS_Shape exitPipe_Seg2;

	TopoDS_Wire scrollStartSideWire1;
	TopoDS_Wire scrollStartSideWire2;
    TopoDS_Wire tongueSweepWire;
    TopoDS_Wire tongueSectWire1;
    TopoDS_Wire tongueSectWire2;

	TopoDS_Shape debugShape;
	TopoDS_Shape debugShape2;
	TopoDS_Shape debugShape3;
	TopoDS_Shape voluteTotal;
    
    bool makeScroll(bool isClockWise, double hub_inlet_z,double shroud_inlet_z, double hub_inlet_radius,
                         double shroud_inlet_radius, double r_offset, double inlet_area, double throat_area,
                         double alpha1,double alpha2,double topRatio);

    bool makeTongueRegionPrism(double hub_inlet_z,double shroud_inlet_z);
    bool makeTongue(bool isClockWise, double hub_inlet_z,double shroud_inlet_z, double m_TailThickness,
        double m_le_aspect_ratio);
    bool makeTongueSolid();

	bool makeTongue_SurfaceFill();
	Handle_Geom_BSplineSurface FillTongueSide(Handle_Geom_Curve pTongueUIso,Handle_Geom_Curve pCrv_top, Handle_Geom_Curve pCrv_bottom, Handle_Geom_Curve pCrv_ScrollStart);
	TopoDS_Shape makeSingleFace(Handle_Geom_BSplineSurface pTongueSideSurf1, Handle_Geom_BSplineSurface pTongueSideSurf2, Handle_Geom_Curve pE2aCrv, 
		Handle_Geom_Curve pE2bCrv);

	TopoDS_Shape makeSingleFaceNew(Handle_Geom_BSplineSurface pTongueSideSurf1, Handle_Geom_BSplineSurface pTongueSideSurf2, Handle_Geom_Curve pE2aCrv, 
		Handle_Geom_Curve pE2bCrv);

    bool Operations();

	Handle_Geom_BSplineCurve makeApproxCurve(Handle_Geom_Curve curve);




	///new functions
	TopoDS_Shape scrollShape1;
	TopoDS_Shape scrollShape2;
	TopoDS_Shape scrollShape3;
	TopoDS_Shape scrollShape4;
	TopoDS_Shape transPipeShape1;
	TopoDS_Shape exitPipeShape1;
	TopoDS_Shape exitPipeShape2;
	TopoDS_Shape fusedJoint;
	TopoDS_Shape trimFace1;
	TopoDS_Shape trimFace2;
	TopoDS_Shape totalVolute; //final product
	TopoDS_Shape totalVoluteThick; //volute will wall thickness
	TopoDS_Shape m_totalVoluteWithoutInOut; //final product
	TopoDS_Shape m_exitPlane;
	TopoDS_Shape m_inputPlane;
	TopoDS_Shape TongueTransitionPatch;//New transition patch method
	TopoDS_Shape VoluteWithOutTransPatch;//New transition patch method

	std::vector<TopoDS_Shape> scrollFaces;
	std::vector<TopoDS_Shape> exitPipeFaces;
	std::vector<TopoDS_Shape> tongueRegionScrollFaces;
	std::vector<TopoDS_Shape> tongueRegionPipeFaces;

	std::vector<TopoDS_Shape> debugList;
	std::vector<double> thetaSects;
	std::vector<TopoDS_Wire> scrollSectWires;
	std::vector<TopoDS_Wire> seg0Wires,seg1Wires, seg2Wires,seg3Wires,seg4Wires,seg5Wires,seg6Wires;
	std::vector<TopoDS_Wire> exitTranSectWires;
	std::vector<TopoDS_Wire> exitTransSpineWires;

	gp_Dir dirScrollEnd;



	bool makeScroll(bool IsClockWise,std::vector<asymSection>& asymSects, double angTran);
	bool modifyAsymStartSegments(bool IsClockWise,double angTran, std::vector<asymSection>& asymSects);
	bool makeScrollSectWires(std::vector<asymSection> asymSects);
	bool makeScrollSolid();
	bool updateScrollEndDirection();

	bool makeExitPipeTransition(bool firstRun, bool overRideTransPos,bool isClockWise,std::vector<asymSection> asymSects,double throat_area,double exit_area, double exit_length,
		bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition,
		double transitionLength, double exit_pipe_aspect1, double firstTransPos, int transSectOption, int firstTransPosOption, int& runFlag, double& calcedSecondTransPos,
		double earlyTurnAngle);


	bool makeExitStartSection(asymSection scrollStartSect, asymSection scrollEndSect, double earlyTurnAngle);
	bool makeExitPipeSpine(bool isClockWise, double exit_length, bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h,
		int exit_curve_order, double* pX_exitCrv,
		double* pY_exitCrv,bool exit_enable_transition, double transitionLength, double firstTransPos, int firstTransPosOption, double earlyTurnAngle);
	bool makeTransitionXsections(bool isClockWise, double throat_area, bool exit_enable_transition, double transitionLength, double exit_area, double exit_pipe_aspect1, double exit_length);
	bool sweepTransitionPipe(bool exit_straight, bool exit_enable_transition);


	bool makeExitPipe(bool isClockWise,double exit_area, bool exit_enable_transition, double exit_pipe_aspect2,bool exit_straight, bool exit_enable_extension, 
		double exit_extension_length, double exit_extension_reduction);
	bool makeExitEndSection(bool isClockWise, bool exit_enable_transition, double exit_pipe_aspect2, double exit_area);
	bool makeExitExtensionSection(double exit_extension_length, double exit_extension_reduction);
	bool sweepExitPipe(bool exit_straight, bool exit_enable_transition, bool exit_enable_extension);


	int tongueBoolean(int filletMethod,bool exit_straight, double exit_length, double transitionLength,double firstTransPos, int firstTransPosOption, double earlyTurnAngle);
	int JoinTransitionPipe2();
	int JoinTransitionPipe2B(int filletMethod,TopoDS_Face face1, TopoDS_Face face2,double exit_length, double transitionLength,double firstTransPos, int firstTransPosOption, double earlyTurnAngle);
	bool sewAllFaces(bool IsClockWise, bool doFillets,int filletMethod,double tongueCurvControlParam,double tongueSizeControlParam, int filletOption, double filletRmin, double filletRmax, bool exit_enable_extension);
	bool makeTongueFillet(TopoDS_Face face1, TopoDS_Face face2, int filletProfile, double R1, double R2, double R3, TopoDS_Shape& filleted, int edgeNum);

	//functions to add wall thickness
	int makeOuterShell(double wallThickness, bool isTurbine, bool exit_straight, bool isClockWise);
	TopoDS_Shape makeTransitionExp(double Offset1,bool exit_straight, TopoDS_Face scrollFace1Exp, TopoDS_Shape exitPipeExp, Handle_Geom_Curve pScrollSectLastExp);
	TopoDS_Wire makeExitStartSectionExp(Handle_Geom_Curve pscrollStartCrv, Handle_Geom_Curve pscrollEndCrv);
	TopoDS_Wire makeExitEndSectionExp(TopoDS_Shape exitPipeExp);
	TopoDS_Shape makeBooleanExp(TopoDS_Face scroll1Face, TopoDS_Face trans1Face, TopoDS_Wire planeWire);
	Handle_Geom_Curve makeProjCurve(Handle_Geom_Curve pCrv, double baseRadius);
	TopoDS_Shape makeFlangePipe(TopoDS_Shape exitPipeExp, double Offset1);
	TopoDS_Shape makeInnerShellNoSplitter(Handle_Geom_Curve pScrollSectLastExp);


	std::vector<TopoDS_Shape> tongueRegionPipeFacesExp;
	TopoDS_Wire transStartWireExp;
	TopoDS_Shape trimFace1Exp;
	TopoDS_Shape trimFace2Exp;
	TopoDS_Shape filletedShapeSaved;


	void writeStepFile();


	public:
		__declspec(dllexport) TopoDS_Shape getScrollShape1();
		__declspec(dllexport) TopoDS_Shape getScrollFace(int faceNum);
		__declspec(dllexport) int getNumScrollFaces();
		__declspec(dllexport) TopoDS_Shape getexitTranShape1();
		__declspec(dllexport) TopoDS_Shape getexitPipeShape1();
		__declspec(dllexport) TopoDS_Shape getexitPipeShape2();
		__declspec(dllexport) TopoDS_Shape gettotalVolute();
		__declspec(dllexport) TopoDS_Shape gettotalVoluteWithoutInOut();
		__declspec(dllexport) TopoDS_Shape getExitPlane();
		__declspec(dllexport) TopoDS_Shape getInputPlane();
	

//volute splitters
		bool makeScrollSplitter(bool IsClockWise,std::vector<asymSection>& asymSects, double angTran);	  
								
		bool makeSplitterVolute(bool IsClockWise,std::vector<asymSection>& asymSects, double angTran,double splitterStartAng,double alpha1,
								double alpha2, double splitterEndDist, std::vector<double>& thicknessAng, std::vector<double>& thicknessVal,
								bool isConstAreaMV, int areaOptionPV, std::vector<double>& areaAng, std::vector<double>& areaRatioPV,
								std::vector<double>& areaAbsPV);
									  
		bool makeStartVoluteMV(bool IsClockWise,std::vector<asymSection>& asymSects, double angTran,double splitterStartAng);
		bool OCCVolute::makeExitPipeTransitionPV(bool firstRun, bool overRideTransPos, bool isClockWise,bool isTurbine,std::vector<asymSection> asymSects,
		double throat_area,double exit_area, double exit_length,bool exit_straight, 
						double exit_pipe_angle_v, double exit_pipe_angle_h, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
                        double exit_pipe_aspect1, double firstTransPos, int transSectOption, int firstTransPosOption, int& runFlag, double& calcedSecondTransPos,
						double earlyTurnAngle);
		bool SplitExitStartSectionPV();	
		bool SplitTransitionXsectionsPVMV(TopoDS_Wire exitTranSectWiresPVMV,TopoDS_Wire exitTransSpineWiresPVMV,
												TopoDS_Wire& WirePV,TopoDS_Wire& WireMV);
		bool ExitStartSectionMV(bool isTurbine);
		bool sweepTransitionPipePVMV(bool exit_straight, bool exit_enable_transition);
		bool sweepExitPipePVMV(bool exit_straight, bool exit_enable_transition, bool exit_enable_extension);
		bool makeExitPipePVMV(bool isClockWise,double exit_area, bool exit_enable_transition, double exit_pipe_aspect2,
				bool exit_straight, bool exit_enable_extension, double exit_extension_length, double exit_extension_reduction, 
				double splitterEndDist, bool isTurbine);
						 
		int tongueBooleanPVMV(int filletMethod,bool exit_straight, double exit_length, double transitionLength, double firstTransPos,
						 int firstTransPosOption, double earlyTurnAngle);
		bool sewAllFacesPVMV(bool IsClockWise, bool doFillets,int filletMethod,double tongueCurvControlParam,double tongueSizeControlParam, int filletOption, double filletRmin, double filletRmax, bool exit_enable_extension);
		//! Separate faces into inlet, exit and mid segments
		//! Assign m_inlet, m_exit, m_hollowVolute
		bool separateFaces();
		bool FindPVThickness(std::vector<asymSection>& asymSects,int index,std::vector<double>& thicknessAng, std::vector<double>& thicknessVal);	
		bool FindAreaRatioPVMV(double splitterStartAng,std::vector<asymSection>& asymSects,int index,bool isConstAreaMV,int areaOptionPV, std::vector<double>& areaAng, 
									std::vector<double>& areaRatioPV,std::vector<double>& areaAbsPV);
		bool ExitSplitterEnd(bool exit_straight, bool exit_enable_transition, bool exit_enable_extension,double splitterEndDist, bool isTurbine);
		TopoDS_Shape makeFilletOnSharpEdge(Handle_Geom_BSplineCurve pCrv1,Handle_Geom_BSplineCurve pCrv2,Handle_Geom_BSplineCurve pCrv3,Handle_Geom_BSplineCurve pCrv4);
		bool MakeSplitterFillet();
		TopoDS_Shape makeInnerShellWithSplitter(Handle_Geom_Curve pScrollSectLastExp);
		TopoDS_Shape makePipeFlangeouterFacePVMV();  //Make exit pipe flance outer face
		bool FindInterSecParametersofSecondEdge(TopoDS_Edge cuttingEdge,TopoDS_Edge TempEdge1,double& Parameter1,double& Parameter2);
		bool makeTongueFilletUsingThrughSections(TopoDS_Face face1, TopoDS_Face face2,double tongueCurvControlParam,double tongueSizeControlParam,TopoDS_Shape& filleted);
		double FindParameterAtaDistanceAlongACurve(Handle_Geom_Curve& CurveUsed,double& Uin,double Distance);
		//New variables for Partial Volute Function
		double PVStartAngle,PVEndAngle ;
		int i_mid;
		bool hasSplitter; //indicates whether volute has a splitter
		bool Spliterendtrue; //Splitter end distance has a non zero value OK
		double ARatioPVMV;	//Area Ratio of the volute sections (Partial Volute)/(Main Volute)
		bool filletOK;	//to check wether Tongue fillet of the Splitter(PV) is ok
		TopoDS_Wire ExitStartSecWirePV;	//Exit Start Section PV wire
		TopoDS_Wire ExitStartSecWireMV;	//Exit Start Section MV wire
		TopoDS_Wire TransExitWirePV;	//Exit Transition section end PV wire
		TopoDS_Wire TransExitWireMV;	//Exit Transition section end MV wire
		TopoDS_Wire ExitWirePV;	//Exit end section end PV wire
		TopoDS_Wire ExitWireMV;	//Exit end section end MV wire
		TopoDS_Shape transPipeShapePV;	//Transition(Exit) pipe shape PV with fillet or without
		TopoDS_Shape transPipeShapeMV;	//Transition(Exit) pipe shape MV
		TopoDS_Shape ExitPipeShapePV;	//Exit pipe shape PV
		TopoDS_Shape ExitPipeShapeMV;	//Exit pipe shape MV
		TopoDS_Shape ExitPipePV;	//Exit pipe Filleted finish shape PV
		TopoDS_Shape ExitPipeMV;	//Exit pipe Filleted finish shape MV
		
		TopoDS_Shape PVMVTong;  //Splitter and Main Volute part
		TopoDS_Shape MVFinal,PVFinal,NobottomPVFinal;		//Main Volute and Partial Volute Parts
		TopoDS_Shape VoluteStartAMV; //Start section Main Volute from 0 deg to 60 deg
		TopoDS_Shape VoluteStartBMV; ////Start section Main Volute from 60 deg to "PVStartAngle"(Partial volute start angle)
		TopoDS_Shape fillettransPipePV;
		TopoDS_Shape fillettransPipeMV;
		TopoDS_Edge InterSecEdge; // this is the edge received after tongueBoolean funtion (cut edge of Face1 and Face2 intersection)
		double Ufirst1,Ulast1;		//Parameters of Tongue curve of PVMVTong need to do the filleting
		std::vector<double> AreaRatioArray;
		std::vector<double> PVThicknessArray;  //"Splitter"(Partial Volute Thicknee:PVThickness) Array
		
		//std::vector<TopoDS_Wire> scrollSectWires;
		std::vector<Handle_Geom_BSplineCurve> PVTsecCurveArray;
		std::vector<Handle_Geom_BSplineCurve> MVsecCurveArray;
		std::vector<Handle_Geom_BSplineCurve> PVsecCurveArray1;
		std::vector<Handle_Geom_BSplineCurve> PVsecCurveArray3;

		//new pump method with more user controllable tonge shape
		bool makeTransitionPatch(bool firstRun, bool overRideTransPos,bool isClockWise,std::vector<asymSection> asymSects, double tongue_le_radius, double tongue_R_aspect_ratio, double tongue_z_aspect_ratio, int & TongueTransEndIndex, double TongueStAngle,
				double throat_area,double exit_area, double exit_length,
				bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition,
				double transitionLength, double exit_pipe_aspect1, double firstTransPos, int transSectOption, int firstTransPosOption, int& runFlag, double& calcedSecondTransPos,
				double earlyTurnAngle);

	//	void GetFaceAandFilletEdge(TopoDS_Face FaceA, TopoDS_Face &NewFaceA, gp_Pnt PntBottom, gp_Pnt PntTop, TopoDS_Edge FilletEdgeSideA, double TongueR);
	//	void GetFaceBandFilletEdge(TopoDS_Face FaceB, TopoDS_Face &NewFaceB, gp_Pnt PntBottom, gp_Pnt PntTop, TopoDS_Edge FilletEdgeSideB);
		TopoDS_Edge TranslateVoluteEndToVotStart(TopoDS_Shape VoluteEndFace, TopoDS_Wire VoluteStartEdge, TopoDS_Wire VoluteEndWire, bool IsPipeShell, gp_Vec & VoluteExitSt_Vec);
		double GetTopPointParameterOfTongueFillet(double TongueR, gp_Pnt Pnt2, gp_Pnt Pnt_P1new, gp_Pnt Pnt4, gp_Pnt Pnt_P2new);
		TopoDS_Shape GetTongueFilletAndNewTransShape(TopoDS_Shape ThroughTransShape, TopoDS_Face & TongueFilletShape, double TongueR, double AspectRatio, double TranslateFactor, 
																					int TongueStartindex, int TongueMidStartindex,  bool Tonguless,
																					gp_Pnt SegABot, gp_Pnt SegATopMid, gp_Pnt SegATopSt, gp_Pnt SegBBot, gp_Pnt SegBTopMid , gp_Pnt SegBTopSt,
																					gp_Pnt TongueBotMidPnt, gp_Pnt TongueTopMidPnt, gp_Vec VoluteExitSt_Vec,
																					TopoDS_Edge SlantSegATop, TopoDS_Edge SlantSegBTop);

		void ModifyFaceA0AndFaceB0(TopoDS_Face & FaceA0, TopoDS_Face & FaceB0, int TongueStartIndex);
		TopoDS_Wire GetVolwireAtMidDistanceToStart(TopoDS_Edge TranslatedVolEdgeMidle, TopoDS_Wire StartMidWire);
		TopoDS_Edge mergeTwoEdges(TopoDS_Edge edge1, TopoDS_Edge edge2);
		bool makeVoluteWithoutPatch(int TongueTransEndIndex);
		bool sewAllNewPatchMethod();
		//! These two methods are coppid from cadgeometrybuilder
		//! Method compare two faces
		//! @param theFace1 first edge to compare
		//! @param theFace2 second edge to compare
		//! @param tol tolerence
		//! @return is equla or not
		bool isSameFace(const TopoDS_Face& theFace1, const TopoDS_Face& theFace2, const double tol = 1e-5);
		//! Method compare two edges
		//! @param theEdge1 first edge to compare
		//! @param theEdge2 second edge to compare
		//! @param tol tolerence
		//! @return is equla or not
		bool isSameEdge(const TopoDS_Edge& theEdge1, const TopoDS_Edge& theEdge2, const double tol = 1e-5);
};