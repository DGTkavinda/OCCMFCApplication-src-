#include "stdafx.h"
#include "Volute.h"
#include "ShapeAnalysis_Shell.hxx"
#include "ShapeAnalysis_ShapeContents.hxx"
#include "ShapeUpgrade_ShellSewing.hxx"
#include "GC_MakeCircle.hxx"
#include "GeomAPI_IntCS.hxx"
#include "GeomFill.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "BRepExtrema_ExtPC.hxx"
#include "BRepOffsetAPI_NormalProjection.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "Geom_Line.hxx"
#include "BRepExtrema_ExtCC.hxx"
#include "GeomLib_Tool.hxx"
#include "GC_MakeEllipse.hxx"
#include "GC_MakeArcOfEllipse.hxx"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "TopTools_MapOfOrientedShape.hxx"
#include "TopExp.hxx"


void OCCVolute::makePumpVolute(bool isClockWise,double hub_inlet_z,double shroud_inlet_z, double hub_inlet_radius,
                         double shroud_inlet_radius, double r_offset, double inlet_area, double throat_area,
                         double alpha1,double alpha2,double topRatio,
                         //wrap region data
                         double wrapLength, double m_TailThickness, double m_le_aspect_ratio,
                         //exit pipe data
                         double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v,
                         int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition,
                         double transitionLength, double exit_pipe_aspect1, double exit_pipe_aspect2,
                         double exit_pipe_lean1, double exit_pipe_lean2)
{
    bool success = makeScroll(isClockWise,hub_inlet_z, shroud_inlet_z, hub_inlet_radius, shroud_inlet_radius,
        r_offset, inlet_area, throat_area, alpha1, alpha2, topRatio);

    success = makeTongue(isClockWise,hub_inlet_z, shroud_inlet_z, m_TailThickness, m_le_aspect_ratio);

    success = makeTongueRegionPrism(hub_inlet_z, shroud_inlet_z);

    success = makeTongueSolid();

	success = makeWrapRegionSurf(wrapLength, m_TailThickness);

	success = makeExitTransitionSurf(isClockWise,throat_area, exit_area, exit_length, exit_straight, exit_pipe_angle_v,
        exit_curve_order, pX_exitCrv, pY_exitCrv, exit_enable_transition, transitionLength, exit_pipe_aspect1,
             exit_pipe_lean1);

	if(exit_enable_transition)
        success = makeExitPipeSurf(isClockWise,exit_area, exit_length, exit_straight, exit_pipe_angle_v,
                         exit_curve_order,pX_exitCrv, pY_exitCrv, exit_enable_transition,
                         transitionLength, exit_pipe_aspect2, exit_pipe_lean2);

    success = Operations();

    return;
}

void OCCVolute::getScrollSurf(Handle_Geom_BSplineSurface& hBsp)
{
    hBsp = pScrollSurf;
    return;
}

void OCCVolute::getexitPipeSurf0(Handle_Geom_BSplineSurface& hBsp)
{
    hBsp = pexitPipeSurf0;
    return;
}

void OCCVolute::getexitPipeSurf1(Handle_Geom_BSplineSurface& hBsp)
{
    hBsp = pexitPipeSurf1;
    return;
}

void OCCVolute::getwrapRegionSurf(Handle_Geom_BSplineSurface& hBsp)
{
    hBsp = pwrapRegionSurf;
    return;
}

void OCCVolute::gettongueSurf(Handle_Geom_BSplineSurface& hBsp)
{
    hBsp = pTongueSurf;
    return;
}


void OCCVolute::clearVolute()
{
    /*delete pScrollSurf;
    delete pTongueSurf;
    delete pwrapRegionSurf;
    delete pexitPipeSurf0;
    delete pexitPipeSurf1;*/
    
}

double OCCVolute::getExitPipe0SectionArea(Standard_Real Vsect)
{
    if(!pexitPipeSurf0)
        return 0.0;

    Handle_Geom_Curve pSectCurve = pexitPipeSurf0->VIso(Vsect);

    Standard_Real U2 = pSectCurve->LastParameter();
    BRepBuilderAPI_MakeWire mkWire;

    gp_Pnt pt0(0.0,0.0,0.0);
    for(int i=0; i< 11; i++)
    {
        Standard_Real Uiter = i*U2/10.0;
        gp_Pnt pt = pSectCurve->Value(Uiter);
        if(i>0)
        {
            BRepBuilderAPI_MakeEdge mkEdge(pt0,pt);
            mkWire.Add(mkEdge);
        }
        pt0 = pt;
    }

    TopoDS_Wire wire = mkWire.Wire();

    BRepBuilderAPI_MakeFace myFace(wire);
    Standard_Boolean isDone = myFace.IsDone();

    double area = 0.0;
    if(isDone)
    {
        GProp_GProps SProps;
        BRepGProp::SurfaceProperties(myFace.Face(),SProps);
        area = SProps.Mass();
    }

    return area;
}

double OCCVolute::getWrapSectionArea(Standard_Real Vsect)
{
    if(!pwrapRegionSurf)
        return 0.0;

    Handle_Geom_Curve pSectCurve = pwrapRegionSurf->VIso(Vsect);

    Standard_Real U2 = pSectCurve->LastParameter();
    gp_Pnt Pt_start = pSectCurve->Value(0.0);
    gp_Pnt Pt_end = pSectCurve->Value(U2);

    BRepBuilderAPI_MakeWire mkWire;

    gp_Pnt pt0(0.0,0.0,0.0);
    for(int i=0; i< 11; i++)
    {
        Standard_Real Uiter = i*U2/10.0;
        gp_Pnt pt = pSectCurve->Value(Uiter);
        if(i>0)
        {
            BRepBuilderAPI_MakeEdge mkEdge(pt0,pt);
            mkWire.Add(mkEdge);
        }
        pt0 = pt;
    }

    TopoDS_Wire wire = mkWire.Wire();

    BRepBuilderAPI_MakeFace myFace(wire);
    Standard_Boolean isDone = myFace.IsDone();

    GProp_GProps SProps;
    BRepGProp::SurfaceProperties(myFace.Face(),SProps);
    double area = SProps.Mass();

    return area;
}

Handle_Geom_BSplineSurface OCCVolute::copySurfaceToNew(Handle_Geom_BSplineSurface& pBSurf)
{
    //make the same surface again
    int NbUPoles = pBSurf->NbUPoles();
    int NbVPoles = pBSurf->NbVPoles();
    TColgp_Array2OfPnt Poles(1,NbUPoles,1,NbVPoles);

    int NbUKnots = pBSurf->NbUKnots();
    TColStd_Array1OfReal UKnots(1,NbUKnots);

    int NbVKnots = pBSurf->NbVKnots();
    TColStd_Array1OfReal VKnots(1,NbVKnots);

    TColStd_Array1OfInteger UMults(1,NbUKnots);
    TColStd_Array1OfInteger VMults(1,NbVKnots);

    

    pBSurf->Poles(Poles);
    pBSurf->UKnots(UKnots);
    pBSurf->VKnots(VKnots);
    pBSurf->UMultiplicities(UMults);
    pBSurf->VMultiplicities(VMults);   

    return new Geom_BSplineSurface(Poles,UKnots,VKnots,UMults,VMults,pBSurf->UDegree(),
        pBSurf->VDegree(),pBSurf->IsUPeriodic(),pBSurf->IsVPeriodic());
}

Handle_Geom_BSplineCurve OCCVolute::copyCurveToNew(Handle_Geom_BSplineCurve& pBCurve)
{
    //make the same surface again
    int NbPoles = pBCurve->NbPoles();
    TColgp_Array1OfPnt Poles(1,NbPoles);

    int NbKnots = pBCurve->NbKnots();
    TColStd_Array1OfReal Knots(1,NbKnots);

    TColStd_Array1OfInteger Mults(1,NbKnots);
 
    pBCurve->Poles(Poles);
    pBCurve->Knots(Knots);
    pBCurve->Multiplicities(Mults);

    return new Geom_BSplineCurve(Poles,Knots,Mults,pBCurve->Degree(),pBCurve->IsPeriodic());
}

Handle_Geom_BSplineCurve OCCVolute::makeSingleCurve(TopoDS_Wire wire,int numPts)
{
    BRepAdaptor_CompCurve BCC(wire);
    Standard_Real U0 = BCC.FirstParameter();
    Standard_Real U1 = BCC.LastParameter();

    int numPt = U1*numPts+1;

    Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,numPt);
	std::vector<gp_Pnt> goodPts;
	gp_Pnt Pnt0(0.0,0.0,0.0);
	int eqlPts = 0;
    for(int k=1; k<=numPt; k++)
    {
        Standard_Real Uiter = (k-1.0)/float(numPts);
        pPoints->SetValue(k,BCC.Value(Uiter));

		if(Pnt0.IsEqual(BCC.Value(Uiter),1.0e-8))
			eqlPts++;
		else
			goodPts.push_back(BCC.Value(Uiter));
		
		Pnt0 = BCC.Value(Uiter);
    }

	if(eqlPts > 0)
	{
		pPoints = new TColgp_HArray1OfPnt(1,goodPts.size());
		for(int k=0; k<goodPts.size(); k++)
			pPoints->SetValue(k+1,goodPts[k]);
	}


    GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-8);
    interp.Perform();
    Handle_Geom_BSplineCurve curve;
    if(interp.IsDone())
            curve = interp.Curve();

    return curve;
}






bool OCCVolute::makeWrapRegionSurf(double wrapLength, double m_TailThickness)
{
	if(wrapLength <= 0.0)
	{
		//set the wrap section end direction to scroll end direction to use in pipe transition
		Handle_Geom_Curve ptempCrv = pScrollSurf->UIso(0.5);
		gp_Pnt tempPt;
		gp_Vec endVec;
		ptempCrv-> D1(0.0,tempPt,endVec);
		dirWrapEnd = gp_Dir(endVec);

		pwrapRegionSurf.Nullify();
		return false;
	}
    //wrap region starting X section is same as tonguetiptopsurf end section
/*    Standard_Real U1,U2,V1,V2;
    pTipTopSurf->Bounds(U1,U2,V1,V2);
    Handle_Geom_Curve pWrapStartSect = pTipTopSurf->VIso(V2);
    BRepBuilderAPI_MakeEdge mkEdge(pWrapStartSect);
    BRepBuilderAPI_MakeWire mkWire(mkEdge);
*/
    TopoDS_Wire wrapSectTotalWire = GetWrapRegionBaseSection();//mkWire.Wire();

    //get the centroid of the starting wrap section
    
    GProp_GProps SProps;
    BRepBuilderAPI_MakeFace mkFace(wrapSectTotalWire);
    BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
    gp_Pnt startSectionCentroid = SProps.CentreOfMass();

    //get the spine for wrap sweeping
    Handle_Geom_Curve pUIsoCurve = pScrollSurf->UIso(0.5);

    GeomAdaptor_Curve GAC(pUIsoCurve);
    Standard_Real Uend = 0.0;
    GCPnts_AbscissaPoint GCPA(GAC,wrapLength,0.0);
    if(GCPA.IsDone())
        Uend = GCPA.Parameter();

    Handle(Geom_BSplineCurve) pbs = Handle(Geom_BSplineCurve)::DownCast(pUIsoCurve);
    pbs->Segment(0.0,Uend);

    //update the end direction of wrap spine
    gp_Vec endVec = pbs-> DN(Uend,1);
    dirWrapEnd = gp_Dir(endVec);

    BRepBuilderAPI_MakeEdge wrapEdge(pbs);
    BRepBuilderAPI_MakeWire wrapWire(wrapEdge);

    //make the wrap region sweep
    BRepOffsetAPI_MakePipeShell wrapRegionShell(wrapWire.Wire());
    wrapRegionShell.Add(wrapSectTotalWire,Standard_False,Standard_False) ;
    Standard_Boolean isReady =wrapRegionShell.IsReady();
    wrapRegionShell.Build();
    BRepBuilderAPI_PipeError wrapStatus = wrapRegionShell.GetStatus();

    wrapRegionShape = wrapRegionShell.Shape();

    //extract the geometric surface
    BRepBuilderAPI_NurbsConvert BNC_wrapRegion(wrapRegionShape);

    TopExp_Explorer Ex;
    TopoDS_Face wrapFace;
    int i=0;
    for (Ex.Init(BNC_wrapRegion.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        wrapFace = TopoDS::Face(Ex.Current());
    }

    Handle_Geom_Surface wrapSurf = BRep_Tool::Surface(wrapFace);
    GeomAdaptor_Surface GASurfwrap(wrapSurf);
    GeomAbs_SurfaceType surfType = GASurfwrap.GetType();
    Handle_Geom_BSplineSurface bwrapSurf = GASurfwrap.BSpline();

    pwrapRegionSurf = copySurfaceToNew(bwrapSurf);

    return true;
}


bool OCCVolute::makeExitTransitionSurf(bool isClockWise,double throat_area, double exit_area, double exit_length,
            bool exit_straight, double exit_pipe_angle_v, int exit_curve_order, double* pX_exitCrv, 
            double* pY_exitCrv, bool exit_enable_transition, double transitionLength, double exit_pipe_aspect1,
            double exit_pipe_lean1)
{
    //MAKE THE EXIT PIPE
    //section up to transition point
    //starting section is the last section of wrap region surface
    Standard_Real U1,U2,V1,V2;
	TopoDS_Wire startSectionWire;
	if(pwrapRegionSurf.IsNull()) //wrap region length is zero. Start from scroll exit
	{
		startSectionWire = GetWrapRegionBaseSection();
	}
	else
	{
		pwrapRegionSurf->Bounds(U1,U2,V1,V2);
		Handle_Geom_Curve pstartSectionCurve = pwrapRegionSurf->VIso(V2);
		gp_Pnt tempPt;
		gp_Vec startSectVec;
		pstartSectionCurve->D1(U1,tempPt,startSectVec);

		BRepBuilderAPI_MakeEdge mkstartSectionEdge(pstartSectionCurve);
		BRepBuilderAPI_MakeWire mkstartSectionWire(mkstartSectionEdge);
		startSectionWire = mkstartSectionWire.Wire();
	}

	//test whether turbine or compressor
	Standard_Boolean isTurbine = false;
	pScrollSurf->Bounds(U1,U2,V1,V2);
    Handle_Geom_Curve pVIsoCurve = pScrollSurf->VIso(V1);
	double Zshrd = pVIsoCurve->Value(U1).Z();
	double Zhub = pVIsoCurve->Value(U2).Z();
	if(Zhub < Zshrd)
		isTurbine = true;



    //make the spine of exit pipe starting from the centroid of start section
	double pipe0_length = exit_length;
    if(exit_enable_transition)
		pipe0_length = exit_length*transitionLength/100.0;
    TopoDS_Edge spineEdge;
    gp_Pnt startSectionCentroid;
    if(exit_straight)
    {
        GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(startSectionWire);
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();

        //find the normal to the stat section to get the spine direction
        BRepAdaptor_Surface tempSurf(mkFace.Face());
        BRepLProp_SLProps slProps(tempSurf,0.0,0.0,1,1.0e-6);
        gp_Dir surfaceNormal = slProps.Normal();
        gp_Vec normalVec(surfaceNormal);
        normalVec.Normalize();

        //test to see the normal vector should be in the +ve direction
        double angleBetween =dirWrapEnd.Angle(gp_Dir(normalVec));
        if(angleBetween > M_PI/2.0) // the normalVec direction has to be reversed
            normalVec.Multiply(-1.0);
               
        normalVec.Multiply(pipe0_length); 
        gp_Pnt end = startSectionCentroid;
        end.Translate(normalVec); 

        BRepBuilderAPI_MakeEdge mkEdge(startSectionCentroid,end);
        spineEdge = mkEdge.Edge(); 
    }
    else
    {
        TColgp_Array1OfPnt CurvePoles(1,exit_curve_order);
        for(int i=1;i<=exit_curve_order; i++)
        {
            gp_Pnt spinePt(pX_exitCrv[i-1], pY_exitCrv[i-1],0.0);
            CurvePoles.SetValue(i,spinePt); 
        }
        Handle(Geom_BezierCurve) pexitSpineCurve = new Geom_BezierCurve(CurvePoles);
 
        GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(startSectionWire);
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();

		//move the spine curve to the centroid of the start section
        gp_Trsf T;
        gp_Pnt P1;
        gp_Vec V1;
        pexitSpineCurve->D1(0.0,P1,V1);//StartPoint();
        T.SetTranslation(P1,startSectionCentroid);
		pexitSpineCurve->Transform(T);

		//align the spine curve with swrirl section spine curve if the curve is towards opposite direction		
		gp_Trsf  Rot; 
		gp_Vec V2(dirWrapEnd);
		Standard_Real rotAngle = V2.Angle(V1);
		if(rotAngle > M_PI/2.0)
		{
			//first rotate perpendicular to OZ direction if necessary
			rotAngle = M_PI;
			gp_Pnt P2(0,0,startSectionCentroid.Z());
			gp_Ax1 rotAxis(P2,gp_Dir(gp_Vec(P2,startSectionCentroid)));		
			Rot.SetRotation(rotAxis,rotAngle);
			pexitSpineCurve->Transform(Rot);

			//now rotate around OZ direction
			rotAxis = gp::OZ();
			T.SetTranslation(gp_Pnt(0,0,0),startSectionCentroid);
			rotAxis = rotAxis.Transformed(T);
			pexitSpineCurve->D1(0.0,P1,V1);
			rotAngle = V2.Angle(V1);		
			Rot.SetRotation(rotAxis,rotAngle);
			pexitSpineCurve->Transform(Rot);
		}

		//get the requierd length of this curve
		GeomAdaptor_Curve GAC(pexitSpineCurve);
		Standard_Real Uend = 0.0;
		U1 = GAC.FirstParameter();
		U2 = GAC.LastParameter();
		GCPnts_AbscissaPoint GCPA(GAC,pipe0_length,U1);
		if(GCPA.IsDone())
			Uend = GCPA.Parameter();

		pexitPipeSpine_Seg2 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		pexitPipeSpine_Seg2->Segment(Uend,U2);

		pexitPipeSpine_Seg1 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		pexitPipeSpine_Seg1->Segment(U1,Uend);

		pexitSpineCurve = pexitPipeSpine_Seg1;

      
        BRepBuilderAPI_MakeEdge mkEdge(pexitSpineCurve);
        spineEdge = mkEdge.Edge(); 
    }
    //spine up to transition end
    BRepBuilderAPI_MakeWire spineWire(spineEdge);

    //if transition is requested make the ellipse shape
    TopoDS_Wire endSectionWire;
    if(exit_enable_transition)
    {
        double targetArea = throat_area+transitionLength*(exit_area-throat_area)/100.0;
        gp_Pnt endSectionCentroid;
        gp_Vec Vec1;
        Standard_Real U1,U2;
        BRep_Tool::Range(spineEdge,U1,U2);
        Handle_Geom_Curve spineCurve = BRep_Tool::Curve(spineEdge,U1,U2);
        spineCurve->D1(U2,endSectionCentroid,Vec1);

        //update the end direction of wrap spine
        dirTransEnd = gp_Dir(Vec1);

        gp_Dir dir(Vec1);
        gp_Ax2 majorAxis(endSectionCentroid,dir);
		gp_Dir Vy;
		if(isClockWise)
		{
			Vy = gp_Dir(gp_Vec(0,0,1));
			if(exit_pipe_aspect1 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}
		else
		{
			Vy = gp_Dir(gp_Vec(0,0,-1));
			if(exit_pipe_aspect1 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}

        double majR = 0.0;
        double minR = 0.0;
        bool found = getTargetEllipse(targetArea, exit_pipe_aspect1, majR, minR);

        //we have the final ellipse shape now
        Geom_Ellipse E_converged(majorAxis,majR,minR);       
        gp_Elips ellipse = E_converged.Elips();

		Handle_Geom_BSplineCurve bspEllipse;
		if(isClockWise)
		{
			if(exit_pipe_aspect1 >= 1.0)
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,0.0); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,-1.0);
			}
			else
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,0.75); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,-0.25); 
			}
		}
		else
		{
			if(exit_pipe_aspect1 >= 1.0) 
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,-1.0); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,0.0); 
			}
			else
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,-0.25);  
				else
					bspEllipse = resetEllipsetoCurve(E_converged,0.75);
			}
		}

        BRepBuilderAPI_MakeEdge mkEdge(bspEllipse/*ellipse*/);
        BRepBuilderAPI_MakeWire mkWire(mkEdge);
        endSectionWire = mkWire.Wire(); 
    }


    //make the exit pipe sweep
    // if exit pipe is straight use the ThruSections algorithm whic is more robust
    // if not straight, use the makePipeShell algorithm which will sweep in a non linear path

    if(!exit_straight) 
    {
        BRepOffsetAPI_MakePipeShell exitPipeShell(spineWire.Wire());
        exitPipeShell.Add(startSectionWire,Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
            exitPipeShell.Add(endSectionWire,Standard_False,Standard_True) ;
        
        Standard_Boolean  isReady =exitPipeShell.IsReady();
        //return true;
        exitPipeShell.Build();
        BRepBuilderAPI_PipeError buildStatus = exitPipeShell.GetStatus();

        exitPipe_Seg1 = exitPipeShell.Shape();
    }
    else
    {
        BRepOffsetAPI_ThruSections thruSections;
        thruSections.AddWire(startSectionWire);
        thruSections.AddWire(endSectionWire);
        thruSections.Build();
        exitPipe_Seg1 = thruSections.Shape();
    }

    //extract the geometric surface
    BRepBuilderAPI_NurbsConvert BNC_exitPipe(exitPipe_Seg1);

    TopExp_Explorer Ex;
    TopoDS_Face exitFace;
    int i=0;
    for (Ex.Init(BNC_exitPipe.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        exitFace = TopoDS::Face(Ex.Current());
    }

    Handle_Geom_Surface exitSurf = BRep_Tool::Surface(exitFace);
    GeomAdaptor_Surface GASurfexit(exitSurf);
    GeomAbs_SurfaceType surfType = GASurfexit.GetType();
    Handle_Geom_BSplineSurface bexitSurf = GASurfexit.BSpline();

    pexitPipeSurf0 = copySurfaceToNew(bexitSurf);


    return true;
}


bool OCCVolute::makeExitPipeSurf(bool isClockWise,double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v,
        int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition, 
        double transitionLength, double exit_pipe_aspect2, double exit_pipe_lean2)
{
    //section from transition point to end

    Standard_Real U1,U2,V1,V2;
    pexitPipeSurf0->Bounds(U1,U2,V1,V2);
    Handle_Geom_Curve pstartSectionCurve = pexitPipeSurf0->VIso(V2);
    gp_Pnt tempPt;
    gp_Vec startSectVec;
    pstartSectionCurve->D1(U1,tempPt,startSectVec);

    BRepBuilderAPI_MakeEdge mkstartSectionEdge(pstartSectionCurve);
    BRepBuilderAPI_MakeWire mkstartSectionWire(mkstartSectionEdge);

	//test whether turbine or compressor
	Standard_Boolean isTurbine = false;
	pScrollSurf->Bounds(U1,U2,V1,V2);
	Handle_Geom_Curve pVIsoCurve = pScrollSurf->VIso(V1);
	double Zshrd = pVIsoCurve->Value(U1).Z();
	double Zhub = pVIsoCurve->Value(U2).Z();
	if(Zhub < Zshrd)
		isTurbine = true;


    //make the spine of exit pipe
    TopoDS_Edge spineEdge;
    gp_Pnt startSectionCentroid;
    if(exit_straight)
    {
        GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(mkstartSectionWire.Wire());
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();
        gp_Pnt base0(0.0,0.0,startSectionCentroid.Z());
        //get the surface to find a normal to it
        gp_Pnt junkPt;
        gp_Vec uVec, vVec;
        pstartSectionCurve->D1(0.0,junkPt,uVec);
        pstartSectionCurve->D1(0.2,junkPt,vVec);
        gp_Vec normalVec = vVec.Crossed(uVec);
        normalVec.Normalize();

        //test to see the normal vector should be in the +ve direction
        double angleBetween =dirTransEnd.Angle(gp_Dir(normalVec));
        if(angleBetween > M_PI/2.0) // the normalVec direction has to be reversed
            normalVec.Multiply(-1.0);

        double pipe1_length = 0.0;
        if(exit_enable_transition)
            pipe1_length = exit_length*(1.0-transitionLength/100.0);
        
        normalVec.Multiply(pipe1_length); 
        gp_Pnt end = startSectionCentroid;
        end.Translate(normalVec); 
 
        BRepBuilderAPI_MakeEdge mkEdge(startSectionCentroid,end);
        spineEdge = mkEdge.Edge();   
    }
    else
    {
		BRepBuilderAPI_MakeEdge mkEdge(pexitPipeSpine_Seg2);
        spineEdge = mkEdge.Edge();        
    }
    //spine from transition up to end
    BRepBuilderAPI_MakeWire spineWire(spineEdge);

    //if transition is requested make the ellipse shape
    TopoDS_Wire endSectionWire;
    if(exit_enable_transition)
    {
        gp_Pnt endSectioncentroid;
        gp_Vec V1;
        Standard_Real U1,U2;
        Handle_Geom_Curve spineCurve = BRep_Tool::Curve(spineEdge,U1,U2);
        spineCurve->D1(U2,endSectioncentroid,V1);

        gp_Dir dir(V1);
        gp_Ax2 majorAxis(endSectioncentroid,dir);
		gp_Dir Vy;
		if(isClockWise)
		{
			Vy = gp_Dir(gp_Vec(0,0,1));
			if(exit_pipe_aspect2 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}
		else
		{
			Vy = gp_Dir(gp_Vec(0,0,-1));
			if(exit_pipe_aspect2 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}

        
        double majR = 0.0;
        double minR = 0.0;
        bool found = getTargetEllipse(exit_area, exit_pipe_aspect2, majR, minR);

        //we have the final ellipse shape now
        Geom_Ellipse E_converged(majorAxis,majR,minR);
        gp_Elips ellipse = E_converged.Elips();

		Handle_Geom_BSplineCurve bspEllipse;
		if(isClockWise)
		{
			if(exit_pipe_aspect2 >= 1.0)
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,0.0); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,-1.0); 
			}
			else
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,0.75); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,-0.25); 
			}
		}
		else
		{
			if(exit_pipe_aspect2 >= 1.0)
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,-1.0); 
				else
					bspEllipse = resetEllipsetoCurve(E_converged,0.0);
			}
			else
			{
				if(isTurbine)
					bspEllipse = resetEllipsetoCurve(E_converged,-0.25);  
				else
					bspEllipse = resetEllipsetoCurve(E_converged,0.75);  
			}

		}

        BRepBuilderAPI_MakeEdge mkEdge(bspEllipse/*ellipse*/);
        BRepBuilderAPI_MakeWire mkWire(mkEdge);
        endSectionWire = mkWire.Wire();
        
    }

    //make the exit pipe sweep
	if(!exit_straight) 
    {
        BRepOffsetAPI_MakePipeShell exitPipeShell(spineWire.Wire());
        exitPipeShell.Add(mkstartSectionWire.Wire(),Standard_False,Standard_False) ;
        exitPipeShell.Add(endSectionWire,Standard_False,Standard_True) ;
        exitPipeShell.Build();
        BRepBuilderAPI_PipeError buildStatus = exitPipeShell.GetStatus();

        exitPipe_Seg2 = exitPipeShell.Shape();
    }
    else
    {
        BRepOffsetAPI_ThruSections thruSections;
		thruSections.AddWire(mkstartSectionWire.Wire());
		thruSections.AddWire(endSectionWire);
		thruSections.Build();
		exitPipe_Seg2 = thruSections.Shape();
    }

    //extract the geometric surface
    BRepBuilderAPI_NurbsConvert BNC_exitPipe(exitPipe_Seg2);

    TopExp_Explorer Ex;
    TopoDS_Face exitFace;
    int i=0;
    for (Ex.Init(BNC_exitPipe.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        exitFace = TopoDS::Face(Ex.Current());
    }

    Handle_Geom_Surface exitSurf = BRep_Tool::Surface(exitFace);
    GeomAdaptor_Surface GASurfexit(exitSurf);
    GeomAbs_SurfaceType surfType = GASurfexit.GetType();
    Handle_Geom_BSplineSurface bexitSurf = GASurfexit.BSpline();

    pexitPipeSurf1 = copySurfaceToNew(bexitSurf);

    return true;
}


bool OCCVolute::getTargetEllipse(double targetArea, double aspectRatio, double& R_major, double& R_minor)
{
    double tolerance = 1.0e-5;

    //make the starting ellipse. note geom_ellipse works only for major radius >= minor radius
    double majR = sqrt(targetArea/M_PI);
    double minR = majR/aspectRatio;
    if(aspectRatio < 1.0)
        minR = majR*aspectRatio;
    gp_Ax2 majorAxis = gp_Ax2(gp_Pnt(0,0,0),gp_Dir(gp_Vec(0,0,1)));
    Geom_Ellipse E0(majorAxis,majR,minR);
    BRepBuilderAPI_MakeEdge mkEdge(E0.Elips());
    BRepBuilderAPI_MakeWire mkWire(mkEdge.Edge());
    BRepBuilderAPI_MakeFace mkFace(mkWire.Wire());

    if(!mkFace.IsDone())
        return false;

    GProp_GProps SProps;
    BRepGProp::SurfaceProperties(mkFace.Shape(),SProps);
    double area = SProps.Mass();

        
    //iterate on the requd area to fix the X section
    bool foundLowLimit = false;
    bool foundHighLimit = false;
    double majR_low = 0.0;
    double majR_high = 0.0;
 
    int j;
    for(j=0;j<100;j++)
    {
        if(foundLowLimit && foundHighLimit)
        {
            if(area < targetArea)
                majR_low = majR;
            else
                majR_high = majR;

            majR = (majR_low+majR_high)/2.0;
        }
        else
        {
            if(area < targetArea)
            {
                foundLowLimit = true;
                majR_low = majR;
                majR = majR*1.1;
            }
            else
            {
                foundHighLimit = true;
                majR_high = majR;
                majR = majR*0.9;
            }
        }

        //recalculate the points
        minR = majR/aspectRatio;
        if(aspectRatio < 1.0)
            minR = majR*aspectRatio;

        //make the section face
        Geom_Ellipse E_iter(majorAxis,majR,minR);
        BRepBuilderAPI_MakeEdge mkEllipseEdge(E_iter.Elips());

        BRepBuilderAPI_MakeWire mkEllipseWire(mkEllipseEdge.Edge());
        BRepBuilderAPI_MakeFace mkEllipseFace(mkEllipseWire.Wire());

        BRepGProp::SurfaceProperties(mkEllipseFace.Shape(),SProps);
        area = SProps.Mass();

        if(fabs(area-targetArea)<tolerance)
            break;

    }

    if(j>= 100) //didn't find the ellipse
        return false;
    else
    {
        R_major = majR;
        R_minor = minR;
    }

    return true;
}

Handle_Geom_BSplineCurve OCCVolute::resetEllipsetoCurve(Geom_Ellipse& E, double fraction)
{
	//first interpolate this to a BSpline curve
	const int numPt = 100;
	Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,numPt+1);

	Standard_Real U1,U2;
	U1 = E.FirstParameter();
    U2 = E.LastParameter();

	for(int i=0; i<=numPt; i++)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/numPt;
		gp_Pnt P;
		E.D0(Uiter,P);
		pPoints->SetValue(i+1,P);	
	}

	GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
    interp.Perform();
    Handle_Geom_BSplineCurve curve;
    if(interp.IsDone())
            curve = interp.Curve();

	if(fraction == 0)
		return curve;
	else if(fraction == -1.0)
	{
		curve->Reverse();
		return curve;
	}
	else if(fraction < 0.0)
	{
		curve->Reverse();
		fraction = -1.0*fraction;
	}

	Handle_Geom_BSplineCurve bSpCrv1 = Handle(Geom_BSplineCurve)::DownCast(curve->Copy()); 
	Handle_Geom_BSplineCurve bSpCrv2 = Handle(Geom_BSplineCurve)::DownCast(curve->Copy());
	

	U1 = bSpCrv1->FirstParameter();
	U2 = bSpCrv1->LastParameter();
	Standard_Real Ubreak = U1+fraction*(U2-U1);
	bSpCrv1->Segment(U1,Ubreak);
	bSpCrv2->Segment(Ubreak,U2);

	BRepBuilderAPI_MakeEdge mkEdge1(bSpCrv1);
	BRepBuilderAPI_MakeEdge mkEdge2(bSpCrv2);
	BRepBuilderAPI_MakeWire mkWire(mkEdge2,mkEdge1);

	curve = makeSingleCurve(mkWire.Wire());

	return curve;
}




//returns a point on the curve corresponding to a given theta value on absolute coordinate
//the corresponding U parameter is updated
gp_Pnt OCCVolute::findPointofTheta(Handle_Geom_Curve& pCurve, double theta, double& U)
{
    gp_Pnt tempPt;
    Standard_Real U1 = pCurve->FirstParameter();
    Standard_Real U2 = pCurve->LastParameter();

    Standard_Real Uup = -1.0;
    Standard_Real Udown = -1.0;
    bool isBracket = false;
    
    gp_Pnt pt_minTheta, pt_maxTheta;
    double theta_min = 1.0e6;
    double theta_max = -1.0e6;

    for(int i=0; i<=10; i++)
    {
        Standard_Real Uiter = U1+i*(U2-U1)/10.0;
        tempPt = pCurve->Value(Uiter);
        double tempTheta = atan2(tempPt.Y(),tempPt.X());

        if(fabs(theta-tempTheta) < 1.0e-4)
        {
            U = Uiter;
            return tempPt;
        }
        if(tempTheta < theta)
            Udown = Uiter;
        else
            Uup = Uiter;

        if(tempTheta > theta_max)
        {
            theta_max = tempTheta;
            pt_maxTheta = tempPt;
        }
        if(tempTheta < theta_min)
        {
            theta_min = tempTheta;
            pt_minTheta = tempPt;
        }


        if(Udown>= 0.0 && Uup >= 0.0)
        {
            isBracket = true;
            break;
        }
    }
    if(!isBracket)
    {
        if(theta < theta_min)
            return pt_minTheta;
        else if(theta > theta_max)
            return pt_maxTheta;
        else
            return gp_Pnt();
    }

    for(int i=0; i<100; i++)
    {
        double Uiter = 0.5*(Uup+Udown);
        tempPt = pCurve->Value(Uiter);
        double tempTheta = atan2(tempPt.Y(),tempPt.X());

        if(fabs(theta-tempTheta) < 1.0e-6)
        {
            U = Uiter;
            break;
        }
        else if(theta > tempTheta)
            Udown = Uiter;
        else
            Uup = Uiter;
    }
        
    return tempPt;
}

gp_Pnt OCCVolute::findPointofRadius(Handle_Geom_Curve& pCurve, double radius, bool fromEnd, double& U)
{
    gp_Pnt tempPt;
    Standard_Real U1 = pCurve->FirstParameter();
    Standard_Real U2 = pCurve->LastParameter();

    Standard_Real Uup = -1.0;
    Standard_Real Udown = -1.0;
    bool isBracket = false;

    for(int i=0; i<=10; i++)
    {
        Standard_Real Uiter;
        if(fromEnd)
            Uiter = U2-i*(U2-U1)/10.0;
        else
            Uiter = U1+i*(U2-U1)/10.0;

        tempPt = pCurve->Value(Uiter);
        double tempRadius = sqrt(tempPt.Y()*tempPt.Y()+tempPt.X()*tempPt.X());

        if(tempRadius < radius)
            Udown = Uiter;
        else
            Uup = Uiter;

        if(Udown>= 0.0 && Uup >= 0.0)
        {
            isBracket = true;
            break;
        }
    }
    if(!isBracket)
        return gp_Pnt();

    for(int i=0; i<100; i++)
    {
        double Uiter = 0.5*(Uup+Udown);
        tempPt = pCurve->Value(Uiter);
        double tempRadius = sqrt(tempPt.Y()*tempPt.Y()+tempPt.X()*tempPt.X());

        if(fabs(radius-tempRadius) < 1.0e-6)
        {
            U = Uiter;
            break;
        }
        else if(radius > tempRadius)
            Udown = Uiter;
        else
            Uup = Uiter;
    }
        
    return tempPt;
}




gp_Pnt OCCVolute::findPointofTheta(Handle_Geom_BSplineSurface& pSurf, double U,double theta, double& V)
{
    gp_Pnt tempPt;
    Standard_Real U1,U2,V1,V2;
    pSurf->Bounds(U1,U2,V1,V2);

    if(U < U1 || U > U2)
        return gp_Pnt();
    
    Standard_Real Vup = -1.0;
    Standard_Real Vdown = -1.0;
    bool isBracket = false;

    gp_Pnt pt_minTheta, pt_maxTheta;
    double theta_min = 1.0e6;
    double theta_max = -1.0e6;

    for(int i=0; i<=10; i++)
    {
        Standard_Real Viter = V2-i*(V2-V1)/10.0;
        pSurf->D0(U,Viter,tempPt);
        double tempTheta = atan2(tempPt.Y(),tempPt.X());

        if(fabs(theta-tempTheta) < 1.0e-6)
        {
            V = Viter;
            return tempPt;
        }

        if(tempTheta < theta)
            Vdown = Viter;
        else
            Vup = Viter;

        if(tempTheta > theta_max)
        {
            theta_max = tempTheta;
            pt_maxTheta = tempPt;
        }
        if(tempTheta < theta_min)
        {
            theta_min = tempTheta;
            pt_minTheta = tempPt;
        }

        if(Vdown>= 0.0 && Vup >= 0.0)
        {
            isBracket = true;
            break;
        }
    }
    if(!isBracket)
    {
        if(theta < theta_min)
            return pt_minTheta;
        else if(theta > theta_max)
            return pt_maxTheta;
        else
            return gp_Pnt();
    }

    for(int i=0; i<100; i++)
    {
        double Viter = 0.5*(Vup+Vdown);
        pSurf->D0(U,Viter,tempPt);
        double tempTheta = atan2(tempPt.Y(),tempPt.X());

        if(fabs(theta-tempTheta) < 1.0e-6)
        {
            V=Viter;
            break;
        }
        else if(theta > tempTheta)
            Vdown = Viter;
        else
            Vup = Viter;
    }
        
    return tempPt;
}



TopoDS_Wire OCCVolute::GetWrapRegionBaseSection()
{
    //MAKE THE WRAP REGION
    Standard_Real U1,U2,V1,V2;
    pScrollSurf->Bounds(U1,U2,V1,V2);
    Handle_Geom_Curve pVIsoCurve = pScrollSurf->VIso(V2);

	Handle(Geom_BSplineCurve) pwrapSect = Handle(Geom_BSplineCurve)::DownCast(pVIsoCurve->Copy());
	pwrapSect->Segment(swirlSectU1,swirlSectU2);
	BRepBuilderAPI_MakeEdge wrapSectEdge(pwrapSect);
	
    
    //get the bottom closing part. this is same as tongue surface end curve
    pTongueSurf->Bounds(U1,U2,V1,V2);
    Handle(Geom_BSplineCurve) pwrapBottomCurve = Handle(Geom_BSplineCurve)::DownCast(pTongueEndCurve->Copy()/*pTongueSurf->VIso(V2)*/);
	U1 = pwrapBottomCurve->FirstParameter(); 
	U2 = pwrapBottomCurve->LastParameter(); 
    Handle_Geom_BSplineCurve pwrapBottom1 = GeomConvert::SplitBSplineCurve(pwrapBottomCurve,U2/2.0, U2,1.0e-3);
    Handle_Geom_BSplineCurve pwrapBottom2 = GeomConvert::SplitBSplineCurve(pwrapBottomCurve,U1, U2/2.0,1.0e-3);

    BRepBuilderAPI_MakeEdge wrapSectBottom1(pwrapBottom1);
    BRepBuilderAPI_MakeEdge wrapSectBottom2(pwrapBottom2);
        

    //the bottom bezier curve is broken in to two and added in two parts so that U=0 point of the 
    // X section will be at the bottom of this curve. This is needed for shape extrusion in next segments
    BRepBuilderAPI_MakeWire wrapSectWire(wrapSectBottom1,wrapSectEdge,wrapSectBottom2);
    
    

    TopoDS_Wire wrapSectTotalWire;
    if(wrapSectWire.IsDone())
    {
        Handle_Geom_BSplineCurve BSPcurve = makeSingleCurve(wrapSectWire.Wire());
        BRepBuilderAPI_MakeEdge mkwrapSectEdge(BSPcurve);
        BRepBuilderAPI_MakeWire mkwrapSectWire(mkwrapSectEdge);
        wrapSectTotalWire = mkwrapSectWire.Wire();
    }

    return wrapSectTotalWire;
}






void OCCVolute::getTestCurve(Handle_Geom_Curve& Crv)
{
    Crv = testCurve;
}


bool OCCVolute::makeScroll(bool isClockWise,double hub_inlet_z,double shroud_inlet_z, double hub_inlet_radius,
                         double shroud_inlet_radius, double r_offset, double inlet_area, double throat_area,
                         double alpha1,double alpha2,double topRatio)
{
    //MAKE THE SCROLL
    //make the inlet shape
    const double D2R = M_PI/180.0;
    const double tolerance = 1.0e-6;
    //this fix is to avoid a bug in occ thru sections algorithm
    if(alpha1 == alpha2)
        alpha2 = 0.999*alpha1;

    const int numSections = 11;
    double startTheta = 0.0*D2R;
    double endTheta = 360.0*D2R;
    double ThetaInterval = (endTheta-startTheta)/(numSections-1.0);
    double areaInterval = (throat_area-inlet_area)/(numSections-1.0);
    
    double ThetaArray[numSections];
    double AreaArray[numSections];
    for(int i=0; i<numSections; i++)
    {
        if(isClockWise)
            ThetaArray[i] = i*ThetaInterval;
        else
            ThetaArray[i] = endTheta-(i*ThetaInterval);

        AreaArray[i] = inlet_area+i*areaInterval;
    }
    

    //get the base shape at zero cross section
    gp_Pnt pt1(shroud_inlet_radius, 0.0, shroud_inlet_z);
    gp_Pnt pt4(hub_inlet_radius, 0.0, hub_inlet_z);
    gp_Pnt pt0(0.5*(shroud_inlet_radius+hub_inlet_radius),0.0,0.5*(shroud_inlet_z+hub_inlet_z));
    double h0start = 0.3*(inlet_area/fabs(hub_inlet_z-shroud_inlet_z));

	gp_Pnt pt2,pt3;
	if(hub_inlet_z > shroud_inlet_z) //compressor and pump cases
	{
		pt2 = gp_Pnt(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z-h0start*tan(alpha1*D2R));
		pt3 = gp_Pnt(hub_inlet_radius+h0start, 0.0,hub_inlet_z+h0start*tan(alpha2*D2R));
	}
	else
	{
		pt2 = gp_Pnt(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z+h0start*tan(alpha1*D2R));
		pt3 = gp_Pnt(hub_inlet_radius+h0start, 0.0,hub_inlet_z-h0start*tan(alpha2*D2R));
	}
//    gp_Pnt pt2(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z-h0start*tan(alpha1*D2R));
//    gp_Pnt pt3(hub_inlet_radius+h0start, 0.0,hub_inlet_z+h0start*tan(alpha2*D2R));
    double r0start1 = (h0start/cos(alpha1*D2R))*topRatio;
    double r0start2 = (h0start/cos(alpha2*D2R))*topRatio;

	gp_Pnt pt2a,pt3a;
	if(hub_inlet_z > shroud_inlet_z) //compressor and pump cases
	{
		pt2a = gp_Pnt(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z-r0start1*sin(alpha1*D2R));
		pt3a = gp_Pnt(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z+r0start2*sin(alpha2*D2R));
	}
	else
	{
		pt2a = gp_Pnt(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z+r0start1*sin(alpha1*D2R));
		pt3a = gp_Pnt(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z-r0start2*sin(alpha2*D2R));
	}

//    gp_Pnt pt2a(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z-r0start1*sin(alpha1*D2R));
//    gp_Pnt pt3a(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z+r0start2*sin(alpha2*D2R));


    TColgp_Array1OfPnt topPts(1,4);
    topPts.SetValue(1,pt2); 
    topPts.SetValue(2,pt2a);
    topPts.SetValue(3,pt3a);
    topPts.SetValue(4,pt3);
    Handle(Geom_BezierCurve) ptopCurve = new Geom_BezierCurve(topPts);

    BRepBuilderAPI_MakeEdge side1Edge(pt1,pt2);
    BRepBuilderAPI_MakeEdge side2Edge(pt3,pt4);
    BRepBuilderAPI_MakeEdge bottomEdge(pt4,pt1);
    BRepBuilderAPI_MakeEdge topEdge(ptopCurve);

    BRepBuilderAPI_MakeWire outline(side1Edge,topEdge,side2Edge,bottomEdge);
    BRepBuilderAPI_MakeFace section(outline);

    GProp_GProps SProps;
    BRepGProp::SurfaceProperties(section.Shape(),SProps);
    double area = SProps.Mass();

    TopoDS_Edge leftSide;
    TopoDS_Edge rightSide;
    TopoDS_Edge topSide;
    TopoDS_Edge bottomSide1;
    TopoDS_Edge bottomSide2;

    BRepOffsetAPI_ThruSections BTS(Standard_False,Standard_False,1.0e-06);
    BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
    
    TopoDS_Wire scrollEndClosedWire;

    for(int i=0; i<numSections; i++)
    {
        double targetArea = AreaArray[i];

        //iterate on the requd area to fix the X section
        bool foundLowLimit = false;
        bool foundHighLimit = false;
        double h_low = 0.0;
        double h_high = 0.0;
        for(int j=0;j<100;j++)
        {
            if(foundLowLimit && foundHighLimit)
            {
                if(area < targetArea)
                    h_low = h0start;
                else
                    h_high = h0start;

                h0start = (h_low+h_high)/2.0;
            }
            else
            {
                if(area < targetArea)
                {
                    foundLowLimit = true;
                    h_low = h0start;
                    h0start = h0start*1.1;
                }
                else
                {
                    foundHighLimit = true;
                    h_high = h0start;
                    h0start = h0start*0.9;
                }
            }

            //recalculate the points
       //     gp_Pnt pt2(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z-h0start*tan(alpha1*D2R));
       //     gp_Pnt pt3(hub_inlet_radius+h0start, 0.0,hub_inlet_z+h0start*tan(alpha2*D2R));
			if(hub_inlet_z > shroud_inlet_z) //compressor and pump cases
			{
				pt2 = gp_Pnt(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z-h0start*tan(alpha1*D2R));
			    pt3 = gp_Pnt(hub_inlet_radius+h0start, 0.0,hub_inlet_z+h0start*tan(alpha2*D2R));
			}
			else
			{
				pt2 = gp_Pnt(shroud_inlet_radius+h0start, 0.0,shroud_inlet_z+h0start*tan(alpha1*D2R));
			    pt3 = gp_Pnt(hub_inlet_radius+h0start, 0.0,hub_inlet_z-h0start*tan(alpha2*D2R));
			}
            r0start1 = (h0start/cos(alpha1*D2R))*topRatio;
            r0start2 = (h0start/cos(alpha2*D2R))*topRatio;
        //    gp_Pnt pt2a(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z-r0start1*sin(alpha1*D2R));
        //    gp_Pnt pt3a(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z+r0start2*sin(alpha2*D2R));
			if(hub_inlet_z > shroud_inlet_z) //compressor and pump cases
			{
				pt2a = gp_Pnt(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z-r0start1*sin(alpha1*D2R));
			    pt3a = gp_Pnt(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z+r0start2*sin(alpha2*D2R));
			}
			else
			{
				pt2a = gp_Pnt(shroud_inlet_radius+r0start1*cos(alpha1*D2R), 0.0, shroud_inlet_z+r0start1*sin(alpha1*D2R));
			    pt3a = gp_Pnt(hub_inlet_radius+r0start2*cos(alpha2*D2R), 0.0, hub_inlet_z-r0start2*sin(alpha2*D2R));
			}

            //make the section face
            TColgp_Array1OfPnt topPts(1,4);
            topPts.SetValue(1,pt2); 
            topPts.SetValue(2,pt2a);
            topPts.SetValue(3,pt3a);
            topPts.SetValue(4,pt3);
            Handle(Geom_BezierCurve) ptopCurve = new Geom_BezierCurve(topPts);
            
            BRepBuilderAPI_MakeEdge side1Edge(pt1,pt2);
            BRepBuilderAPI_MakeEdge side2Edge(pt3,pt4);
            BRepBuilderAPI_MakeEdge bottomEdge(pt4,pt1);
            BRepBuilderAPI_MakeEdge topEdge(ptopCurve);

            BRepBuilderAPI_MakeWire outline(side1Edge,topEdge,side2Edge,bottomEdge);
            BRepBuilderAPI_MakeFace section(outline);

            GProp_GProps SProps;
            BRepGProp::SurfaceProperties(section.Shape(),SProps);
            area = SProps.Mass();

            if(fabs(area-targetArea)<tolerance)
            {
                leftSide = TopoDS::Edge(side1Edge.Shape());
                rightSide = TopoDS::Edge(side2Edge.Shape());
                topSide = TopoDS::Edge(topEdge.Shape());
                //bottomSide = TopoDS::Edge(bottomEdge.Shape());

                if(i == numSections-1)
                    scrollEndClosedWire = outline.Wire(); //save this last section wire to start exit pipe

                if(i== 0)
				{
                    pTongueStartCurve = ptopCurve; //save this first section top curve to start tongue surface
					scrollStartSideWire1 = TopoDS::Wire(BRepBuilderAPI_MakeWire(side1Edge));
					scrollStartSideWire2 = TopoDS::Wire(BRepBuilderAPI_MakeWire(side2Edge));
				}


                break;
            }

        }

        //we have the target area. So make the X section
        //need to split the bottomSide at the middle to align the curve at the bottom
        bottomSide1 = BRepBuilderAPI_MakeEdge(pt0,pt1);
        bottomSide2 = BRepBuilderAPI_MakeEdge(pt4,pt0);

        //BRepBuilderAPI_MakeWire Xsection(bottomSide1,leftSide,topSide,rightSide);
        //Xsection.Add(bottomSide2);
        BRepBuilderAPI_MakeWire Xsection(leftSide,topSide,rightSide);


        gp_Trsf  rotate; 
        rotate.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Vec(0,0,1)),ThetaArray[i]);
        BRepBuilderAPI_Transform rotated(Xsection,rotate); 
        TopoDS_Wire rotatedWire = TopoDS::Wire(rotated.Shape());

        Handle_Geom_BSplineCurve curve = makeSingleCurve(rotatedWire);
        
        BRepBuilderAPI_MakeEdge totalEdge(curve);
        BRepBuilderAPI_MakeWire totalWire(totalEdge);
        
		BTS.AddWire(totalWire); 
    }

    //make the scroll shell
    BTS.Build();
    TopoDS_Shape scrollShape= BTS.Shape();

    BRepBuilderAPI_NurbsConvert BNC(scrollShape);

    TopExp_Explorer Ex; 
    TopoDS_Face scrollFace;
    i=0;
    for (Ex.Init(BNC.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        if(i==1)
            scrollFace = TopoDS::Face(Ex.Current());
    }

    Handle_Geom_Surface scrollSurf = BRep_Tool::Surface(scrollFace);
    GeomAdaptor_Surface GASurf(scrollSurf);
    GeomAbs_SurfaceType surfType = GASurf.GetType();
    Handle_Geom_BSplineSurface bScrollSurf = GASurf.BSpline();

     pScrollSurf = copySurfaceToNew(bScrollSurf);

    return true;

}


bool OCCVolute::makeTongueRegionPrism(double hub_inlet_z,double shroud_inlet_z)
{
    Standard_Boolean status;
    //base circle arc
    Handle_Geom_Curve pVIso = pScrollSurf->VIso(0.98);
    gp_Pnt Pt1;
    pVIso->D0(0.0,Pt1);
    pVIso = pScrollSurf->VIso(0.99);
    gp_Pnt Pt2;
    pVIso->D0(0.0,Pt2);
    pVIso = pScrollSurf->VIso(1.0);
    gp_Pnt Pt3;
    pVIso->D0(0.0,Pt3);
    
    GC_MakeArcOfCircle mkArc(Pt1,Pt2,Pt3);
    BRepBuilderAPI_MakeEdge mkEdge1(mkArc.Value());
    testCurve = mkArc.Value();
    

    //vertial line
    gp_Pnt tempPt;
    pVIso->D0(0.5,tempPt);
    gp_Pnt Pt4(tempPt.X(),tempPt.Y(),Pt3.Z());
    BRepBuilderAPI_MakeEdge mkEdge2(Pt3,Pt4);
    status = mkEdge2.IsDone();

    Standard_Real U1,U2;
    BRep_Tool::Range(mkEdge2.Edge(),U1,U2);
    Handle_Geom_Curve verLine = BRep_Tool::Curve(mkEdge2.Edge(),U1,U2);
//    gp_Pnt Pt5;
//    verLine->D0(0.2*U2,Pt5);
//    BRepBuilderAPI_MakeEdge mkEdge3(Pt3,Pt5);

    pTongueEndCurve->D0(0.0,tempPt);
    gp_Pnt Pt5(tempPt.X(),tempPt.Y(),Pt3.Z());
    BRepBuilderAPI_MakeEdge mkEdge3(Pt3,Pt5);

    

    //top arc
    Handle_Geom_Curve tempCrv = mkArc.Value();
    tempCrv = tempCrv->Reversed();
    U1 = tempCrv->FirstParameter();
    gp_Vec Vec1;
    tempCrv->D1(U1,tempPt,Vec1);
    GC_MakeArcOfCircle mkArc2(Pt5,Vec1,Pt1);
    BRepBuilderAPI_MakeEdge mkEdge4(mkArc2.Value());

    gp_Pnt Center(0,0,Pt3.Z());
    BRepBuilderAPI_MakeEdge tedge1(Center,Pt5);
    BRepBuilderAPI_MakeEdge tedge2(Center,Pt1);
    BRepBuilderAPI_MakeWire mkWire(tedge1,tedge2,mkEdge4);
    //BRepBuilderAPI_MakeWire mkWire(mkEdge1,mkEdge3,mkEdge4);


    //end test
    BRepBuilderAPI_MakeFace mkFace(mkWire.Wire());
    TopoDS_Face prismFace = mkFace.Face();
    status = mkFace.IsDone();


    //now we have the shape. make the surface
    gp_Trsf  t1; 
    t1.SetTranslation(gp_Vec(0,0,-3.0*(hub_inlet_z-shroud_inlet_z)));
    BRepBuilderAPI_Transform trf1(t1);
    trf1.Perform(prismFace,false);
    TopoDS_Shape trfShape = trf1.Shape();

    //make a surface 
    gp_Vec		aPrismVec(0 , 0 , 6.0*(hub_inlet_z-shroud_inlet_z)); 
    TopoDS_Shape prism = BRepPrimAPI_MakePrism(trfShape,aPrismVec);
    tonguePrism = prism;

    return true;
}


bool OCCVolute::makeTongue(bool isClockWise, double hub_inlet_z,double shroud_inlet_z, double m_TailThickness,
        double m_le_aspect_ratio)
{
    Standard_Real U1,U2,V1,V2;
    //get the z center of starting curve
    double zCenter = 0.5*(hub_inlet_z+shroud_inlet_z);
    //locate the point in starting curve
    gp_Pnt startCenter;
    U1 = pTongueStartCurve->FirstParameter();
    U2 = pTongueStartCurve->LastParameter();
    pTongueStartCurve->D0(0.5*(U1+U2),startCenter);

    
    double lastSectBottomRadi = sqrt(startCenter.X()*startCenter.X()+startCenter.Y()*startCenter.Y()) + m_TailThickness;

    //locate this radi point on the last section
    //make the last curve section of the tongue using section at 360 degree
    pScrollSurf->Bounds(U1,U2,V1,V2);
    Handle_Geom_Curve pVIsoCurve = pScrollSurf->VIso(V2);
    double Ur = 0.0;
    double Us = 0.0;
    gp_Pnt lastSectPt1;
    gp_Pnt lastSectPt3;
    lastSectPt1 = findPointofRadius(pVIsoCurve, 1.05*lastSectBottomRadi, false, Ur);
    lastSectPt3 = findPointofRadius(pVIsoCurve, 1.05*lastSectBottomRadi, true, Us);
    //update variables to use in the swirl section
    swirlSectU1 = Ur;
    swirlSectU2 = Us;

    //get a top Z center point
    
    gp_Pnt tempPt;
    pVIsoCurve->D0(0.5*(U1+U2),tempPt);
    gp_Pnt zCenterTop(tempPt.X(),tempPt.Y(),startCenter.Z());
    
    //make a line from the Z center points. the tongue bottom end section center lies on this line
    Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,2);
    pPoints->SetValue(1,startCenter);
    pPoints->SetValue(2,zCenterTop);

    GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
    interp.Perform();
    Handle_Geom_BSplineCurve centerLine;
    if(interp.IsDone())
            centerLine = interp.Curve();

    GeomAdaptor_Curve GAC(centerLine);
    Standard_Real Uend = 0.0;
    U1 = GAC.FirstParameter();
    GCPnts_AbscissaPoint GCPA(GAC,m_TailThickness,U1);
    if(GCPA.IsDone())
        Uend = GCPA.Parameter();

    gp_Pnt endCenter;
    gp_Vec endVector;
    centerLine->Segment(U1,Uend);
    centerLine->D1(Uend,endCenter, endVector); 
    gp_Pnt spineEllipseCenter; //center of the ellipse that make up the spine for tongue sweep
    centerLine->D0(0.5*(U1+Uend),spineEllipseCenter);


    //we have 3 points on the last tongue curve. make the curve
	pPoints = new TColgp_HArray1OfPnt(1,3);
	pPoints->SetValue(1,lastSectPt3);
	pPoints->SetValue(2,endCenter);
	pPoints->SetValue(3,lastSectPt1);

	GeomAPI_Interpolate interpLastSect(pPoints,Standard_False,1.0e-6);
    interpLastSect.Perform();
    Handle_Geom_BSplineCurve pCrv;
    if(interpLastSect.IsDone())
		pCrv = interpLastSect.Curve();
    pTongueEndCurve = copyCurveToNew(pCrv);
   

    //make the spine for tongue surface sweep
    gp_Ax2 Axis = gp_Ax2(spineEllipseCenter,gp_Dir(gp_Vec(0,0,1)));

    gp_Dir Vy(endVector);
    if(m_le_aspect_ratio >= 1.0) //depending on where the major radius would be
        Axis.SetYDirection(Vy);
    else
        Axis.SetXDirection(Vy);

    double R1 = spineEllipseCenter.Distance(startCenter);
    double R2 = R1*m_le_aspect_ratio;


    pPoints = new TColgp_HArray1OfPnt(1,51);
    if(m_le_aspect_ratio >= 1.0)
    {
        Geom_Ellipse spineEllipse(Axis,R2,R1);
        for(int i=1; i<=51; i++)
        {
            gp_Pnt Pt;
            double Uiter;
            if(isClockWise)
                Uiter = 1.5*M_PI+(i-1)*M_PI/50.0;
            else
                Uiter = 1.5*M_PI-(i-1)*M_PI/50.0;

            spineEllipse.D0(Uiter,Pt);
            pPoints->SetValue(i,Pt);
        }

        GeomAPI_Interpolate interpSpine(pPoints,Standard_False,1.0e-6);
        interpSpine.Perform();
        Handle_Geom_BSplineCurve hCrv;
        if(interpSpine.IsDone())
            hCrv = interpSpine.Curve();
        pTongueSpineCurve = copyCurveToNew(hCrv);
        
    }
    else
    {
        Geom_Ellipse spineEllipse(Axis,R1,R2);
        for(int i=1; i<=51; i++)
        {
            gp_Pnt Pt;
            double Uiter;
            if(isClockWise)
                Uiter = M_PI+(i-1)*M_PI/50.0;//2.0*M_PI-(i-1)*M_PI/50.0;//(i-1)*M_PI/50.0;
            else
                Uiter = M_PI-(i-1)*M_PI/50.0; //2.0*M_PI-(i-1)*M_PI/50.0;

            spineEllipse.D0(Uiter,Pt);
            pPoints->SetValue(i,Pt);
        }

        GeomAPI_Interpolate interpSpine(pPoints,Standard_False,1.0e-6);
        interpSpine.Perform();
        Handle_Geom_BSplineCurve hCrv;
        if(interpSpine.IsDone())
            hCrv = interpSpine.Curve();
        pTongueSpineCurve = copyCurveToNew(hCrv);

    }

    

	//generate the tongue center surface. We take segments of the pTongueStartCurve & pTongueEndCurve for this
	double Startfraction = 0.1;
	double Endfraction = 0.2;
	U1 = pTongueStartCurve->FirstParameter();
	U2 = pTongueStartCurve->LastParameter();
	Handle_Geom_BezierCurve pTongueStart_mini = Handle(Geom_BezierCurve)::DownCast(pTongueStartCurve->Copy()); 
	Standard_Real Ustart = U1+ Startfraction*(U2-U1);
	Uend = U1+(1.0-Startfraction)*(U2-U1);
	pTongueStart_mini->Segment(Ustart,Uend);

	U1 = pTongueEndCurve->FirstParameter();
	U2 = pTongueEndCurve->LastParameter();
	Handle_Geom_BSplineCurve pTongueEnd_mini = Handle(Geom_BSplineCurve)::DownCast(pTongueEndCurve->Copy()); 
	Ustart = U1+ Endfraction*(U2-U1);
	Uend = U1+(1.0-Endfraction)*(U2-U1);
	pTongueEnd_mini->Segment(Ustart,Uend);


    //make the wires and do the tongue surface sweep
    BRepBuilderAPI_MakeEdge startSectEdge(pTongueStart_mini);
    BRepBuilderAPI_MakeWire startSectWire(startSectEdge);

    BRepBuilderAPI_MakeEdge endSectEdge(pTongueEnd_mini);
    BRepBuilderAPI_MakeWire endSectWire(endSectEdge);

    BRepBuilderAPI_MakeEdge spineEdge(pTongueSpineCurve);
    BRepBuilderAPI_MakeWire spineWire(spineEdge);

    //make the tongue center face sweep
    BRepOffsetAPI_MakePipeShell tongueSurfSweep(spineWire.Wire());
    tongueSurfSweep.Add(startSectWire.Wire(),Standard_False,Standard_False) ;
    tongueSurfSweep.Add(endSectWire.Wire(),Standard_False,Standard_False) ;
    
    tongueSurfSweep.Build();
    BRepBuilderAPI_PipeError tongueStatus = tongueSurfSweep.GetStatus();

    tongueShape = tongueSurfSweep.Shape();
    tongueSweepWire = spineWire.Wire();
    tongueSectWire1 = startSectWire.Wire();
    tongueSectWire2 = endSectWire.Wire();

    //extract the geometric surface
    BRepBuilderAPI_NurbsConvert BNC_tongueSurf(tongueShape);

    TopExp_Explorer Ex;
    TopoDS_Face tongueFace;
    int i=0;
    for (Ex.Init(BNC_tongueSurf.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        tongueFace = TopoDS::Face(Ex.Current());
    }

    Handle_Geom_Surface tongueSurf = BRep_Tool::Surface(tongueFace);
    GeomAdaptor_Surface GASurftongue(tongueSurf);
    GeomAbs_SurfaceType surfType = GASurftongue.GetType();
    Handle_Geom_BSplineSurface btongueSurf = GASurftongue.BSpline();

    pTongueSurf = copySurfaceToNew(btongueSurf);




    return true;
}

bool OCCVolute::makeTongueSolid()
{
	makeTongue_SurfaceFill();
    
    return true;
}

Handle_Geom_BSplineSurface OCCVolute::FillTongueSide(Handle_Geom_Curve pCrv1,Handle_Geom_Curve pCrv2, Handle_Geom_Curve pCrv3, Handle_Geom_Curve pCrv4)
{
	Standard_Real Tol3d = 1.0e-6;
    Standard_Real Tolang = 1.0e-3;

	Handle_GeomAdaptor_HSurface H0= new GeomAdaptor_HSurface(pTongueSurf);
    Adaptor3d_CurveOnSurface CurveOnSurf(H0);
    Handle_GeomFill_BoundWithSurf B0 = new GeomFill_BoundWithSurf(CurveOnSurf,Tol3d,Tolang);

	Handle_GeomAdaptor_HCurve H1 = new GeomAdaptor_HCurve(pCrv1);
    Handle_GeomFill_SimpleBound B1 = new GeomFill_SimpleBound(H1,Tol3d,Tolang);

    Handle_GeomAdaptor_HCurve H2 = new GeomAdaptor_HCurve(pCrv2);
    Handle_GeomFill_SimpleBound B2 = new GeomFill_SimpleBound(H2,Tol3d,Tolang);

    Handle_GeomAdaptor_HCurve H3 = new GeomAdaptor_HCurve(pCrv3);
    Handle_GeomFill_SimpleBound B3 = new GeomFill_SimpleBound(H3,Tol3d,Tolang);

    Handle_GeomAdaptor_HCurve H4 = new GeomAdaptor_HCurve(pCrv4);
    Handle_GeomFill_SimpleBound B4 = new GeomFill_SimpleBound(H4,Tol3d,Tolang);

    Standard_Integer MaxDeg = 12;
    Standard_Integer MaxSeg = 12;
    GeomFill_ConstrainedFilling GCF(MaxDeg,MaxSeg);
    GCF.Init(B1,B2,B3,B4) ;

	GCF.SetDomain(1.0,B0);
    GCF.ReBuild();

	Handle_Geom_BSplineSurface pSurf = GCF.Surface();
    
	return pSurf;
}


bool OCCVolute::makeTongue_SurfaceFill()
{
	//do the boolean operations
	Standard_Real U1,U2,V1,V2;
    Standard_Boolean success;
    BRepBuilderAPI_MakeFace face1(pScrollSurf,1.0e-6);
	
 /*   BOPTools_DSFiller df;
    df.SetShapes (face1.Shape(),tonguePrism);
    df.Perform();
*/ 
    BRepAlgoAPI_Cut makeBool(face1.Shape(),tonguePrism/*,df*/);
    makeBool.Build();
    success = makeBool.IsDone();
    TopoDS_Shape shape = makeBool.Shape();
	cutScroll = shape;


    BRepAlgoAPI_Common makeBool2(face1.Shape(),tonguePrism/*,df*/);
    makeBool2.Build();
    success = makeBool2.IsDone();
    TopoDS_Shape tongueSide = makeBool2.Shape();

    //get the side surfaces
    TopoDS_Face sideFace1;
    TopoDS_Face sideFace2;
    int i=0;
	TopExp_Explorer Ex;
    for (Ex.Init(tongueSide,TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
        if(i==1)
            sideFace1 = TopoDS::Face(Ex.Current());
        if(i==2)
            sideFace2 = TopoDS::Face(Ex.Current());
    }

	//get the curves that make up the side surfaces
    i = 0;
	TopoDS_Edge E1a,E2a;
	TopoDS_Shape sh1,sh2;
    for (Ex.Init(sideFace1,TopAbs_EDGE); Ex.More(); Ex.Next())
    { 
        i++;
		if(i==1) sh1 = (Ex.Current());
		if(i==2) sh2 = (Ex.Current());
    }
	E1a = TopoDS::Edge(sh1);
	E2a = TopoDS::Edge(sh2);

    BRep_Tool::Range(E1a,U1,U2);
    Handle_Geom_Curve pTempCrv = BRep_Tool::Curve(E1a,U1,U2);
    Handle(Geom_BSplineCurve) pE1aCrv = Handle(Geom_BSplineCurve)::DownCast(pTempCrv->Copy());
    pE1aCrv->Segment(U1,U2);
    pE1aCrv->Reverse();

    BRep_Tool::Range(E2a,U1,U2);
    pTempCrv = BRep_Tool::Curve(E2a,U1,U2);
    Handle(Geom_BSplineCurve) pE2aCrv = Handle(Geom_BSplineCurve)::DownCast(pTempCrv->Copy());
    pE2aCrv->Segment(U1,U2);
    pE2aCrv->Reverse();


    i = 0;
	TopoDS_Edge E1b,E2b;
    for (Ex.Init(sideFace2,TopAbs_EDGE); Ex.More(); Ex.Next())
    { 
        i++;
        if(i==1) sh1 = (Ex.Current());
        if(i==2) sh2 = (Ex.Current());
    }
	E1b = TopoDS::Edge(sh1);
	E2b = TopoDS::Edge(sh2);

    BRep_Tool::Range(E1b,U1,U2);
    pTempCrv = BRep_Tool::Curve(E1b,U1,U2);
    Handle(Geom_BSplineCurve) pE1bCrv = Handle(Geom_BSplineCurve)::DownCast(pTempCrv->Copy());
    pE1bCrv->Segment(U1,U2);
    pE1bCrv->Reverse();

    BRep_Tool::Range(E2b,U1,U2);
    pTempCrv = BRep_Tool::Curve(E2b,U1,U2);
    Handle(Geom_BSplineCurve) pE2bCrv = Handle(Geom_BSplineCurve)::DownCast(pTempCrv->Copy());
    pE2bCrv->Segment(U1,U2);
 //   pE2bCrv->Reverse(); //do not reverse this curve it is already in the correct orientation



	Standard_Real Ustart, Uend;
	double Startfraction = 0.1;
	double Endfraction = 0.2;
	BRepBuilderAPI_MakeEdge mke1(pE1aCrv);
	BRepBuilderAPI_MakeEdge mke2(pE2aCrv);
	BRepBuilderAPI_MakeWire mkW1(mke1,mke2);
	Handle_Geom_BSplineCurve bsp1 = makeSingleCurve(mkW1.Wire());

	U1 = pTongueStartCurve->FirstParameter();
	U2 = pTongueStartCurve->LastParameter();
	Handle_Geom_BezierCurve ptempCrv1 = Handle(Geom_BezierCurve)::DownCast(pTongueStartCurve->Copy()); 
	Ustart = U1;
	Uend = U1+Startfraction*(U2-U1);
	ptempCrv1->Segment(Ustart,Uend);
	BRepBuilderAPI_MakeEdge mkEd1(ptempCrv1);
	BRepBuilderAPI_MakeWire mkWire1(scrollStartSideWire1);
	mkWire1.Add(mkEd1.Edge());
	Handle_Geom_BSplineCurve ptempCrv2 = makeSingleCurve(mkWire1.Wire());

	Handle_Geom_BSplineCurve ptempCrv3 = Handle(Geom_BSplineCurve)::DownCast(pTongueEndCurve->Copy());
	U1 = ptempCrv3->FirstParameter();
	U2 = ptempCrv3->LastParameter();
	Ustart = U1+(1.0-Endfraction)*(U2-U1);
	ptempCrv3->Segment(Ustart, U2);

	pTongueSurf->Bounds(U1,U2,V1,V2);

	Handle_Geom_BSplineSurface pTongueSideSurf1 = FillTongueSide(pTongueSurf->UIso(U1),bsp1, ptempCrv2, ptempCrv3);
	BRepBuilderAPI_MakeFace mkFace1(pTongueSideSurf1,1.0e-6);	
    tongueFinal1 = mkFace1.Face();


	BRepBuilderAPI_MakeEdge mke3(pE1bCrv);
	BRepBuilderAPI_MakeEdge mke4(pE2bCrv);
	BRepBuilderAPI_MakeWire mkW2(mke3,mke4);
	Handle_Geom_BSplineCurve bsp2 = makeSingleCurve(mkW2.Wire());

	U1 = pTongueStartCurve->FirstParameter();
	U2 = pTongueStartCurve->LastParameter();
	ptempCrv1 = Handle(Geom_BezierCurve)::DownCast(pTongueStartCurve->Copy()); 
	Ustart = U1+(1.0-Startfraction)*(U2-U1);
	Uend = U2;
	ptempCrv1->Segment(Ustart,Uend);
	BRepBuilderAPI_MakeEdge mkEd2(ptempCrv1);
	BRepBuilderAPI_MakeWire mkWire2(scrollStartSideWire2);
	mkWire2.Add(mkEd2.Edge());
	ptempCrv2 = makeSingleCurve(mkWire2.Wire());

	ptempCrv3 = Handle(Geom_BSplineCurve)::DownCast(pTongueEndCurve->Copy());
	U1 = ptempCrv3->FirstParameter();
	U2 = ptempCrv3->LastParameter();
	Ustart = U1;
	Uend = U1+(Endfraction)*(U2-U1);
	ptempCrv3->Segment(Ustart, Uend);

	pTongueSurf->Bounds(U1,U2,V1,V2);

	Handle_Geom_BSplineSurface pTongueSideSurf2 = FillTongueSide(pTongueSurf->UIso(U2),bsp2, ptempCrv2, ptempCrv3);
	BRepBuilderAPI_MakeFace mkFace2(pTongueSideSurf2,1.0e-6);	
    tongueFinal2 = mkFace2.Face();

	tongueFinal = makeSingleFaceNew(pTongueSideSurf1,pTongueSideSurf2, pE2aCrv,pE2bCrv );


	return true;
}

TopoDS_Shape OCCVolute::makeSingleFace(Handle_Geom_BSplineSurface pTongueSideSurf1, Handle_Geom_BSplineSurface pTongueSideSurf2,
									   Handle_Geom_Curve pE2aCrv, Handle_Geom_Curve pE2bCrv )
{
	Standard_Real U1,U2,V1,V2;

	BRepOffsetAPI_ThruSections BTS;
	BTS.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pScrollSurf->VIso(0.0))));

	const int numSects = 7;
	pTongueSurf->Bounds(U1,U2,V1,V2);
	Handle_Geom_Curve ptempCrv1  = pTongueSurf->UIso(U1);
	Standard_Real Ustart = V1;
	Standard_Real Udiv = (V2-V1)/1000.0;
	
	for(int i=1; i<=numSects; i++)
	{
		pTongueSideSurf1->Bounds(U1,U2,V1,V2);
		Standard_Real Uiter = U1+i*(U2-U1)/(numSects+1.0);
		Handle_Geom_Curve ptempCrv2 = pTongueSideSurf1->UIso(Uiter);
		Handle_Geom_Curve ptempCrv3 = pTongueSideSurf2->UIso(Uiter);
		gp_Pnt endPt;
		ptempCrv2->D0(V1,endPt);

		//locate this point on the Crv1;
		Standard_Real Uiter2;
		for(int j=0; j<500; j++)
		{
			Uiter2 = Ustart+j*Udiv;
			gp_Pnt Pt2;
			ptempCrv1->D0(Uiter2,Pt2);
			Standard_Real dist = Pt2.Distance(endPt);
			if(fabs(dist) < 1.0e-4)
			{
				Ustart = Uiter2;
				break;
			}
		}

		Handle_Geom_Curve ptempCrv4  = pTongueSurf->VIso(Uiter2);

		BRepBuilderAPI_MakeEdge mkE1(ptempCrv2);
		BRepBuilderAPI_MakeEdge mkE2(ptempCrv3);
		BRepBuilderAPI_MakeEdge mkE3(ptempCrv4);

		//make a single curve
		int numPt = 261;
		int index = 1;
		Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,numPt);

		ptempCrv2->Reverse();
		U1 = ptempCrv2->FirstParameter();
		U2 = ptempCrv2->LastParameter();
		for(int j=0; j<90; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			ptempCrv2->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		U1 = ptempCrv4->FirstParameter();
		U2 = ptempCrv4->LastParameter();
		for(int j=10; j<91; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			ptempCrv4->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		U1 = ptempCrv3->FirstParameter();
		U2 = ptempCrv3->LastParameter();
		for(int j=11; j<=100; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			ptempCrv3->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
		interp.Perform();
		Handle_Geom_BSplineCurve aproxCrv;
		if(interp.IsDone())
            aproxCrv = interp.Curve();

		BTS.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(aproxCrv)));
	}

	//add the last section
	pTongueSurf->Bounds(U1,U2,V1,V2);
	Handle_Geom_Curve pSect2 = pTongueSurf->VIso(V2);


	//make a single curve
		int numPt = 301;
		int index = 1;
		Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,numPt);

		
		
		gp_Pnt P2;
		pE2aCrv->D0(0.0,P2);

		U1 = pE2aCrv->FirstParameter();
		U2 = pE2aCrv->LastParameter();
		for(int j=0; j<=100; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			pE2aCrv->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		gp_Pnt P1;
		pSect2->D0(0.0,P1);
		U1 = pSect2->FirstParameter();
		U2 = pSect2->LastParameter();
		for(int j=1; j<=100; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			pSect2->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		pE2bCrv->Reverse();
		gp_Pnt P0;
		pE2bCrv->D0(0.0,P0);
		U1 = pE2bCrv->FirstParameter();
		U2 = pE2bCrv->LastParameter();
		for(int j=1; j<=100; j++)
		{
			Standard_Real Ui = U1+j*(U2-U1)/100.0;
			gp_Pnt tempPt;
			pE2bCrv->D0(Ui,tempPt);
			pPoints->SetValue(index,tempPt);
			index++;
		}

		GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
		interp.Perform();
		Handle_Geom_BSplineCurve aproxCrv;
		if(interp.IsDone())
            aproxCrv = interp.Curve();

		BTS.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(aproxCrv)));


	BTS.Build();
	
	return BTS.Shape();



}





TopoDS_Shape OCCVolute::makeSingleFaceNew(Handle_Geom_BSplineSurface pTongueSideSurf1, Handle_Geom_BSplineSurface pTongueSideSurf2,
									   Handle_Geom_Curve pE2aCrv, Handle_Geom_Curve pE2bCrv )
{
	BRepBuilderAPI_MakeFace mkF1(pTongueSideSurf1,1.0e-6);
	BRepBuilderAPI_MakeFace mkF2(pTongueSideSurf2,1.0e-6);
	BRepBuilderAPI_MakeFace mkF3(pTongueSurf,1.0e-6);

	BRepBuilderAPI_Sewing sewing;
    sewing.Add(mkF3.Face());
	sewing.Add(mkF1.Face());
	sewing.Add(mkF2.Face());
    
    sewing.Perform(); 
    TopoDS_Shape shape = sewing.SewedShape(); 

	return shape;

//	BRepAlgoAPI_Fuse fuse1(mkF1.Shape(),mkF3.Shape());
//	BRepAlgoAPI_Fuse fuse2(fuse1.Shape(), mkF2.Shape());

//	return fuse2.Shape();
}


TopoDS_Shape OCCVolute::getScrollShape()
{
	return cutScroll;
}

TopoDS_Shape OCCVolute::getTongueShape()
{
	return tongueFinal;
}

TopoDS_Shape OCCVolute::getVoluteTotal()
{
	return voluteTotal;
	
}

TopoDS_Shape OCCVolute::getVoluteTotalThick()
{
	return totalVoluteThick;
}


Handle_Geom_BSplineCurve OCCVolute::makeApproxCurve(Handle_Geom_Curve curve)
{
    Standard_Real U1 = curve->FirstParameter();
    Standard_Real U2 = curve->LastParameter();

    int numPt = 11;

    Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,numPt);
    for(int k=1; k<=numPt; k++)
    {
        Standard_Real Uiter = U1+ (U2-U1)*(k-1.0)/(numPt-1.0);
		gp_Pnt tempPt;
		curve->D0(Uiter,tempPt);
        pPoints->SetValue(k,tempPt);
    }

    GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
    interp.Perform();
    Handle_Geom_BSplineCurve aproxCrv;
    if(interp.IsDone())
            aproxCrv = interp.Curve();

    return aproxCrv;
}

bool OCCVolute::Operations()
{
	//sew everything together
	BRepBuilderAPI_Sewing sewing;

	sewing.Add(cutScroll);
	sewing.Add(wrapRegionShape);
	sewing.Add(exitPipe_Seg1);
	sewing.Add(exitPipe_Seg2);
 
    sewing.Perform(); 
    voluteTotal = sewing.SewedShape();

    

    //writing a step file
/*    STEPControl_StepModelType type = STEPControl_AsIs;
    IFSelect_ReturnStatus status;
    STEPControl_Writer writer;

  //    status = writer.Transfer( debugShape , type );
 //   status = writer.Transfer( debugShape2 , type );
//    status = writer.Transfer( debugShape3 , type );
//	status = writer.Transfer( debugShape4 , type );
//	status = writer.Transfer( debugShape5 , type );
 //   status = writer.Transfer( cutScroll , type );
//    status = writer.Transfer( tongueShape , type );
//    status = writer.Transfer( tongueFinal1 , type );
    status = writer.Transfer( voluteTotal , type );
 	status = writer.Transfer( tongueFinal , type );
    

    Standard_CString filename = "test.step";
    status = writer.Write( filename );
*/

    return true;
}

/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

TopoDS_Shape OCCVolute::getScrollShape1()
{
	return scrollShape1;
}

int OCCVolute::getNumScrollFaces()
{
	return scrollFaces.size();
}

TopoDS_Shape OCCVolute::getScrollFace(int faceNum)
{
	if(scrollFaces.size()> faceNum)
		return scrollFaces[faceNum];
	else
	{
		TopoDS_Shape nullShape;
		nullShape.Nullify();
		return nullShape;
	}
}

TopoDS_Shape OCCVolute::getexitTranShape1()
{
	return transPipeShape1;
}

TopoDS_Shape OCCVolute::getexitPipeShape1()
{
	return exitPipeShape1;
}

TopoDS_Shape OCCVolute::getexitPipeShape2()
{
	return exitPipeShape2;
}

TopoDS_Shape OCCVolute::gettotalVolute()
{
	return totalVolute;
}

TopoDS_Shape OCCVolute::gettotalVoluteWithoutInOut()
{
	return m_totalVoluteWithoutInOut;
}

TopoDS_Shape OCCVolute::getExitPlane()
{
	return m_exitPlane;
}

TopoDS_Shape OCCVolute::getInputPlane()
{
	return m_inputPlane;
}

int OCCVolute::makePumpVolute2(bool doBoolean, bool doFillets, std::vector<asymSection>& asymSects, double angTran, bool IsClockWise, bool isTurbine,
							double throat_area,
						 bool firstRun, bool overRideTransPos, int& runFlag, double& calcedSecondTransPos, double earlyTurnAngle, double alpha1,double alpha2,
						 //exit pipe data
                         double exit_area, double exit_length, bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h,
                         int exit_curve_order, double* pX_exitCrv, 
                         double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
                         double exit_pipe_aspect1, double exit_pipe_aspect2, int firstTransPosOption, double firstTransPos, int transSectOption,
						 bool exit_enable_extension, double exit_extension_length, double exit_extension_reduction,
						 //tongue data
						 int filletMethod,double tongueCurvControlParam,double tongueSizeControlParam, int filletOption, double filletRmin, double filletRmax,
						 double tongue_le_radius, double tongue_Raspect_ratio, double tongue_z_aspect_ratio, bool Istongueless, double TongueStAngle,
						 //splitter data
						  bool isSplitter, double splitterStartAng, double splitterEndDist, std::vector<double>& thicknessAng, std::vector<double>& thicknessVal,
						  bool isConstAreaMV, int areaOptionPV, std::vector<double>& areaAng, std::vector<double>& areaRatioPV, std::vector<double>& areaAbsPV)
{
	int errorCode = 0;
	bool success = false;
	hasSplitter = isSplitter;
	

	if(firstRun)
	{
		if(isSplitter)
			{success = makeSplitterVolute(IsClockWise,asymSects, angTran, splitterStartAng,alpha1, alpha2, splitterEndDist, thicknessAng, thicknessVal,
									 isConstAreaMV, areaOptionPV, areaAng, areaRatioPV, areaAbsPV);
			makeStartVoluteMV(IsClockWise,asymSects, angTran,splitterStartAng);
			}
		else
			success = makeScroll(IsClockWise,asymSects, angTran );
			
		if(!success)
			errorCode = 1;
	}
	else
		success = true;
		
	


		

	if(success)
	{	if(isSplitter)
		success = makeExitPipeTransitionPV(firstRun,overRideTransPos,IsClockWise,isTurbine,asymSects,throat_area,exit_area, exit_length, exit_straight, exit_pipe_angle_v, exit_pipe_angle_h, exit_curve_order, pX_exitCrv, pY_exitCrv, 
			exit_enable_transition, transitionLength, exit_pipe_aspect1,firstTransPos, transSectOption,firstTransPosOption, runFlag, calcedSecondTransPos, 
			earlyTurnAngle);
		else 
		success = makeExitPipeTransition(firstRun,overRideTransPos,IsClockWise,asymSects,throat_area,exit_area, exit_length, exit_straight, exit_pipe_angle_v, exit_pipe_angle_h, exit_curve_order, pX_exitCrv, pY_exitCrv, 
			exit_enable_transition, transitionLength, exit_pipe_aspect1,firstTransPos, transSectOption,firstTransPosOption, runFlag, calcedSecondTransPos, 
			earlyTurnAngle);

		if(success && runFlag == 1)
			return errorCode;

		if(!success)
			errorCode = 2;
	}

	if(success)
	{	if(isSplitter)
		success = makeExitPipePVMV(IsClockWise,exit_area, exit_enable_transition, exit_pipe_aspect2,exit_straight,exit_enable_extension, 
								exit_extension_length, exit_extension_reduction,splitterEndDist, isTurbine);
		else 
		success = makeExitPipe(IsClockWise,exit_area, exit_enable_transition, exit_pipe_aspect2,exit_straight,exit_enable_extension, exit_extension_length, exit_extension_reduction);

		if(!success)
			errorCode = 3;
	}

	//work on new transition patch method
	bool useNewPatchMethod = true;
	if(useNewPatchMethod && Istongueless && !isSplitter ){
		
		
		int TongueTransEndIndex; //this indicate Last scrollsec wire index used to make Tongue Transition patch

		bool successnew = makeTransitionPatch(firstRun,overRideTransPos,IsClockWise,asymSects, tongue_le_radius, tongue_Raspect_ratio, tongue_z_aspect_ratio,
			 TongueTransEndIndex,  TongueStAngle,
			throat_area,exit_area, exit_length, exit_straight, exit_pipe_angle_v, exit_pipe_angle_h, exit_curve_order, pX_exitCrv, pY_exitCrv, 
			exit_enable_transition, transitionLength, exit_pipe_aspect1,firstTransPos, transSectOption,firstTransPosOption, runFlag, calcedSecondTransPos, 
			earlyTurnAngle);

	// Make Volute without Volute start and Tongue area 
	{
		successnew =makeVoluteWithoutPatch(TongueTransEndIndex);
		
	}

		{// sew all faces of new patch method to make totalVolute


			successnew =sewAllNewPatchMethod();

		}
		//end work on new transition patch method
		if(successnew) separateFaces();
	return 0;
	}
	
	
	//perform the final boolean operations only if the user requests
	if(doBoolean && success)
	{	if(isSplitter)
		{
		//	makeStartVoluteMV(IsClockWise,asymSects, angTran,splitterStartAng);
			errorCode =tongueBooleanPVMV(filletMethod,exit_straight, exit_length, transitionLength, firstTransPos, firstTransPosOption, earlyTurnAngle);
			if(errorCode == 0)
				sewAllFacesPVMV(IsClockWise, doFillets,filletMethod,tongueCurvControlParam,tongueSizeControlParam, filletOption, filletRmin, filletRmax,exit_enable_extension);  
		}
		else 
		{
			errorCode = tongueBoolean(filletMethod,exit_straight, exit_length, transitionLength, firstTransPos, firstTransPosOption, earlyTurnAngle);
			if(errorCode == 0)
				sewAllFaces(IsClockWise, doFillets,filletMethod,tongueCurvControlParam, tongueSizeControlParam, filletOption, filletRmin, filletRmax,exit_enable_extension);  
		}
	}
	else if(success) //if no Boolean op was done make a compound of the shape collection
	{
		TopoDS_Compound comp1_totalVolute;
		BRep_Builder builder1_totalVolute;
		builder1_totalVolute.MakeCompound( comp1_totalVolute );
		builder1_totalVolute.Add(comp1_totalVolute,scrollShape1);
		if(scrollFaces.size() > 0)
		{
			for(int i=0; i< scrollFaces.size(); i++)
				builder1_totalVolute.Add(comp1_totalVolute,scrollFaces[i]);
		}
		else if(!scrollShape2.IsNull()) 
		{
			builder1_totalVolute.Add(comp1_totalVolute,scrollShape2);
		}
		if(!transPipeShape1.IsNull())
			builder1_totalVolute.Add(comp1_totalVolute,transPipeShape1);
		if(!exitPipeShape1.IsNull())
			builder1_totalVolute.Add(comp1_totalVolute,exitPipeShape1);
		if(!exitPipeShape2.IsNull())
			builder1_totalVolute.Add(comp1_totalVolute,exitPipeShape2);

		totalVolute = comp1_totalVolute;

	}

	if((!errorCode == 0 || !doBoolean) && success)
	{
		if(isSplitter&&!doFillets)
		{
		
		sewAllFacesPVMV(IsClockWise, doFillets,filletMethod,tongueCurvControlParam,tongueSizeControlParam, filletOption, filletRmin, filletRmax,exit_enable_extension); 
		}
		else if(isSplitter&&!doBoolean)
			{
				sewAllFacesPVMV(IsClockWise, false,filletMethod,tongueCurvControlParam,tongueSizeControlParam, filletOption, filletRmin, filletRmax,exit_enable_extension); 
			}
	}

	if(!errorCode) separateFaces();
	//writeStepFile();


	return errorCode;
}


bool OCCVolute::makeScroll(bool IsClockWise, std::vector<asymSection>& asymSects, double angTran)
{
	bool success = false;
	//reset some data
	scrollShape1.Nullify();
	scrollShape2.Nullify();
	scrollShape3.Nullify();
	scrollShape4.Nullify();
/*	transPipeShape1.Nullify();
	transPipeShape2.Nullify();
	exitPipeShape1.Nullify();
	totalVolute.Nullify();
	trimFace7.Nullify();
	trimFace8.Nullify();
	trimFace9.Nullify();
	trimFace10.Nullify();
*/
// 	modifyAsymStartSegments(IsClockWise,angTran, asymSects);
	makeScrollSectWires(asymSects);
	
	

// 	modifyScrollStartSegments(angTran);
// 	redoScrollStartSegments(asymSects, angTran);

	success = makeScrollSolid();
	

 	updateScrollEndDirection();
 
	return success;
}
 
bool OCCVolute::modifyAsymStartSegments(bool IsClockWise,double angTran, std::vector<asymSection>& asymSects)
{
	//pick the cross section close to angTran angle.
	//Replace all the sections below that with that cross section

	double degree = 5.0;

 	angTran = std::max(degree*M_PI/180.0, angTran);

	int index = 0;
	for(int i=0; i<asymSects.size();i++)
	{
		if(asymSects[i].thetaSect >= angTran)
		{
			index = i;
			break;
		}
	}
	

	gp_Trsf  rotate;
	int numPt = asymSects[index].npt;
	int numLandmark = asymSects[index].landmark.size();
	for(int i=0; i<index; i++)
	{
		double rotAngle = asymSects[i].thetaSect-asymSects[index].thetaSect;
		if(!IsClockWise)
			rotAngle = -1.0*rotAngle;

		asymSects[i].npt = numPt;

		asymSects[i].landmark.resize(numLandmark);
		for(int j=0; j<numLandmark; j++)
			asymSects[i].landmark[j] = asymSects[index].landmark[j];

		for(int j=0; j<asymSects[index].npt; j++)
		{
			gp_Pnt tempPt = gp_Pnt(asymSects[index].x[j],asymSects[index].y[j],asymSects[index].z[j]);
			tempPt.Rotate(gp_Ax1(gp_Pnt(0,0,0),gp_Vec(0,0,1)),rotAngle);
		 	asymSects[i].x[j] = tempPt.X();
		 	asymSects[i].y[j] = tempPt.Y();
		 	asymSects[i].z[j] = tempPt.Z();
		}
	}

	return true;
}


bool OCCVolute::makeScrollSectWires(std::vector<asymSection> asymSects)
{
	int size = asymSects.size();

	thetaSects.resize(size);
	seg0Wires.resize(size);
	seg1Wires.resize(size);
	seg2Wires.resize(size);
	seg3Wires.resize(size);
	seg4Wires.resize(size);
	seg5Wires.resize(size);
	seg6Wires.resize(size);

	scrollSectWires.resize(size);

	int numSeg = 1;
	int segRange[7];

	for(int i=0;i<size;i++)
	{
 		segRange[0] = 0;
	//	segRange[1] = (asymSects[i].landmark[0]+segRange[0])/2; //asymSects[i].landmark[0];
	//	segRange[2] = asymSects[i].landmark[10];
	//	segRange[3] = asymSects[i].landmark[9];
	//	segRange[4] = asymSects[i].landmark[4];
	//	segRange[5] = asymSects[i].landmark[5];
		segRange[1] = asymSects[i].npt-1;


		for(int j=0; j<numSeg; j++)
		{
			int numPt = segRange[j+1]-segRange[j]+1;

			//no need to get all the points. Let us get every 3rd point
			int ptCount = 0;
			for(int k=segRange[j]; k<=segRange[j+1]; k++)
			{
				if (k == segRange[j] || k == segRange[j+1])
				{
					ptCount++;
				}
				else if(k%3 == 0)
				{
					ptCount++;
				}
			}

			Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,ptCount/*numPt*/);
		
			int index = 1;
			for(int k=segRange[j]; k<=segRange[j+1]; k++)
			{
				bool addThisPoint = false;
				if (k == segRange[j] || k == segRange[j+1])
				{
					addThisPoint = true;
				}
				else if(k%3 == 0)
				{
					addThisPoint = true;
				}
				
				if(addThisPoint)
				{
					gp_Pnt tempPt(asymSects[i].x[k],asymSects[i].y[k],asymSects[i].z[k]);
					pPoints->SetValue(index,tempPt);
					index++;
				}
			}
			
					
			GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
			interp.Perform();
			Handle_Geom_BSplineCurve pCurve;
			if(interp.IsDone())
			{
				pCurve = interp.Curve();
				BRepBuilderAPI_MakeEdge mkEdge(pCurve);
				BRepBuilderAPI_MakeWire mkWire(mkEdge);

				if(j==0) seg0Wires[i] = mkWire.Wire();
		 //		if(j==1) seg1Wires[i] = mkWire.Wire();
		 //		if(j==2) seg2Wires[i] = mkWire.Wire();
		//		if(j==3) seg3Wires[i] = mkWire.Wire();
		//		if(j==4) seg4Wires[i] = mkWire.Wire();
		//		if(j==5) seg5Wires[i] = mkWire.Wire();

				
			}
		}

		int firstIndex = segRange[0];
		gp_Pnt Pt_start(asymSects[i].x[firstIndex],asymSects[i].y[firstIndex],asymSects[i].z[firstIndex]);
		int lastIndex = segRange[1];//asymSects[i].npt-1;
		gp_Pnt Pt_end(asymSects[i].x[lastIndex],asymSects[i].y[lastIndex],asymSects[i].z[lastIndex]);
		BRepBuilderAPI_MakeEdge mkE(Pt_end,Pt_start);
		BRepBuilderAPI_MakeWire mkW(mkE);
		seg6Wires[i] = mkW.Wire();


		BRepBuilderAPI_MakeWire mkW2(seg0Wires[i]);
//	 	mkW2.Add(seg1Wires[i]);
//	 	mkW2.Add(seg2Wires[i]);
	//	mkW2.Add(seg3Wires[i]);
	//	mkW2.Add(seg4Wires[i]);
	//	mkW2.Add(seg5Wires[i]);
  	 	mkW2.Add(seg6Wires[i]);
			
		
		scrollSectWires[i] = mkW2.Wire();
		thetaSects[i] = asymSects[i].thetaSect;	

	}

	return true;
}


bool OCCVolute::makeScrollSolid()
{
	bool success = false;
// 	BRepOffsetAPI_ThruSections BTS1(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS3(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS4(Standard_False,Standard_False,1.0e-06);

	BRepOffsetAPI_ThruSections BTS2B(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS3B(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS4B(Standard_False,Standard_False,1.0e-06);



	int i_mid = 0;
	int i_break = 0;
	int i_break_st = 0;
	int divider1 = 3;
	int divider2 = 5;

	double angleBreakDegrees = 60.0;
	double angleBreak = angleBreakDegrees*M_PI/180.0; //M_PI/3.0;
	for(int i=0; i<scrollSectWires.size();i++)
	{
		if(thetaSects[i] > angleBreak)
		{
			i_mid = i-1;
			break;
		}
    }



 
	////
	// The first try block contain functionality to build scroll surfaces. This method is fast but less robust.
	//If that fails the function in corresponding catch block will try to do the scroll surfaces. This second medthod is more robust but slower.
	try
	{
		int counter = i_mid;
		Standard_Boolean isRuled = Standard_False;	// need smooth surfaces not ruled ones

		while(counter < scrollSectWires.size())
		{
			BRepOffsetAPI_ThruSections BTSscrollSegA(Standard_False,isRuled,1.0e-08);
			BRepOffsetAPI_ThruSections BTSscrollSegB(Standard_False,isRuled,1.0e-08);
			int numAdded = 0;
	 		BTSscrollSegA.SetContinuity(GeomAbs_C1);
	 		BTSscrollSegB.SetContinuity(GeomAbs_C1);
			
			while(true)
			{
				if(counter == i_mid || counter%3 == 0 || counter == scrollSectWires.size()-1)
				{
					BTSscrollSegA.AddWire(scrollSectWires[counter]);
					numAdded++;
				}
				if(counter == i_mid || counter%1 == 0 || counter == scrollSectWires.size()-1)
				{
					BTSscrollSegB.AddWire(scrollSectWires[counter]);
				}
		
	 		  //	if(numAdded >= 5 && scrollSectWires.size()-counter > 3)
	 		  //		break;

				counter++;

				if(scrollSectWires.size()-counter < 1)
					break;
			}

			try
			{
				BTSscrollSegA.Build();
				scrollFaces.push_back(BTSscrollSegA.Shape());					
			}
			catch(...)
			{
				BTSscrollSegB.Build();
				scrollFaces.push_back(BTSscrollSegB.Shape());
			}

		}
	}
	catch(...) ////if the first function fails try the second method to build scroll surfaces. This is more robust but slower
	{
		//return false;  
		scrollFaces.clear();	//remove any faces added in the previous try

		//2nd quarter
		for(int i=i_mid; (i<scrollSectWires.size() && thetaSects[i] <= M_PI) ;i++)
		{
			if(i%divider1 == 0 || i == i_mid)
			{
				BTS2.AddWire(scrollSectWires[i]);
				i_break = i;
			}				
		}

		try
		{
			BTS2.Build();
			scrollShape2 = BTS2.Shape();
		}
		catch(...) 
		{
			for(int i=i_mid; (i<scrollSectWires.size() && thetaSects[i] <= M_PI) ;i++)
			{
				if(i%divider2 == 0 || i == i_mid)
				{
					BTS2B.AddWire(scrollSectWires[i]);
					i_break = i;
				}				
			}
			
			try
			{
				BTS2B.Build();
				scrollShape2 = BTS2B.Shape();
			}
			catch(...)
			{
				return false;
			}
			
		}

 
		//3rd quarter
		i_break_st = i_break;
		for(int i=i_break; (i<scrollSectWires.size() && thetaSects[i] <= 1.5*M_PI) ;i++)
		{
			if(i%divider1 == 0 || i == i_break) //get only altrenate sections in the second half 
			{
				BTS3.AddWire(scrollSectWires[i]);
				i_break = i;
			}			
		}

		try
		{
			BTS3.Build();
			scrollShape3 = BTS3.Shape();
		}
		catch(...)
		{
			for(int i=i_break_st; (i<scrollSectWires.size() && thetaSects[i] <= 1.5*M_PI) ;i++)
			{
				if(i%divider2 == 0 || i == i_break_st) //get only altrenate sections in the second half 
				{
					BTS3B.AddWire(scrollSectWires[i]);
					i_break = i;
				}			
			}

			try
			{
				BTS3B.Build();
				scrollShape3 = BTS3B.Shape();
			}
			catch(...)
			{
				return false;
			}


		}

		//4th quarter
		i_break_st = i_break;

		for(int i=i_break; i< scrollSectWires.size();i++)
		{
			if(i%divider1 == 0 || i== scrollSectWires.size()-1 || i == i_break) //get only altrenate sections in the second half 
				BTS4.AddWire(scrollSectWires[i]);				
		}

		try
		{
			BTS4.Build();
			scrollShape4 = BTS4B.Shape();
		}
		catch(...)
		{
			for(int i=i_break_st; i< scrollSectWires.size() ;i++)
			{
				if(i%divider2 == 0 || i== scrollSectWires.size()-1 || i == i_break_st) //get only altrenate sections in the second half 
					BTS4B.AddWire(scrollSectWires[i]);
			}

			try
			{
				BTS4B.Build();
				scrollShape4 = BTS4B.Shape();
			}
			catch(...)
			{
				return false;
			}
		}

 		scrollFaces.push_back(scrollShape2);
 		scrollFaces.push_back(scrollShape3);
 		scrollFaces.push_back(scrollShape4);
		
	}
	////end of 2nd method to build scroll surface
	

	//make the first quarter
	int i_st = 0;
	int numMid = 2; //when the # of sections used in creating the surface is increased, the Boolean operation takes a longer time (numMid was originally set to 5)


	double gap = angleBreak/double(numMid);
	BRepOffsetAPI_ThruSections thruSects(Standard_False,Standard_False,1.0e-06);
	thruSects.AddWire(scrollSectWires[i_st]);
	for(int i=1; i< numMid; i++)
	{
		double ang = thetaSects[0]+i*gap;
		for(int j=0; i<scrollSectWires.size();j++)
		{
			if(ang <= thetaSects[j])
			{
				thruSects.AddWire(scrollSectWires[j]);
				break;
			}
		}
	}
	thruSects.AddWire(scrollSectWires[i_mid]);

	try
	{
		thruSects.Build();
		scrollShape1 = thruSects.Shape(); 
	}
	catch(...)
	{
		//return false;
		// try with a different number of sections
		numMid = 7;
		gap = angleBreak/double(numMid);
		BRepOffsetAPI_ThruSections thruSects2(Standard_False,Standard_False,1.0e-06);
		thruSects2.AddWire(scrollSectWires[i_st]);
		for(int i=1; i< numMid; i++)
		{
			double ang = thetaSects[0]+i*gap;
			for(int j=0; i<scrollSectWires.size();j++)
			{
				if(ang <= thetaSects[j])
				{
					thruSects2.AddWire(scrollSectWires[j]);
					break;
				}
			}
		}
		thruSects2.AddWire(scrollSectWires[i_mid]);

		try
		{
			thruSects2.Build();
			scrollShape1 = thruSects2.Shape(); 
		}
		catch(...)
		{
			return false;
		}
	}
	
	//Find the bottom faces
	BRepBuilderAPI_Sewing sewing1(0.005/*1.0e-6*/);
	//TopoDS_Face face1,face2;
	for(int j = 0; j < scrollFaces.size() ; j++)
	{
		TopoDS_Shape ashape = scrollFaces[j];
		TopoDS_Shape subShape;
		TopAbs_ShapeEnum subType;
		int facecount = 0;

		for(TopExp_Explorer explr(ashape,TopAbs_FACE);explr.More(); explr.Next())
		{
			//TopoDS_Face tempFace = TopoDS::Face(explr.Current());
			
			 if(facecount == 1)
			{
				//face1 = tempFace;
				sewing1.Add(TopoDS::Face(explr.Current()));
				break;
			}
			facecount++;
		}
	}
	int facecount = 0;
	for(TopExp_Explorer explr(scrollShape1,TopAbs_FACE);explr.More(); explr.Next())
	{
		//TopoDS_Face tempFace = TopoDS::Face(explr.Current());

		if(facecount == 1)
		{
			//face2 = tempFace;
			sewing1.Add(TopoDS::Face(explr.Current()));
			break;
		}
		facecount++;
	}

	sewing1.Perform();
	m_inputPlane = sewing1.SewedShape();

	success = true;

	return success;
}

bool OCCVolute::updateScrollEndDirection()
{
	TopoDS_Vertex V0,V1,V2;
	for(int i=0; i<3; i++)
	{
		TopExp_Explorer Ex;
		
		int index = scrollSectWires.size()-3+i;
    
		for (Ex.Init(scrollSectWires[index],TopAbs_VERTEX); Ex.More(); Ex.Next())
		{ 
			if(i==0) V0 = TopoDS::Vertex(Ex.Current());
			if(i==1) V1 = TopoDS::Vertex(Ex.Current());
			if(i==2) V2 = TopoDS::Vertex(Ex.Current());
			break;
		}

	}

	gp_Pnt P0 = BRep_Tool::Pnt(V0);
	gp_Pnt P1 = BRep_Tool::Pnt(V1);
	gp_Pnt P2 = BRep_Tool::Pnt(V2);

	Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,3);

	pPoints->SetValue(1,P0);	
	pPoints->SetValue(2,P1);
	pPoints->SetValue(3,P2);
	

	GeomAPI_Interpolate interp(pPoints,Standard_False,1.0e-6);
    interp.Perform();
    Handle_Geom_BSplineCurve curve;
    if(interp.IsDone())
            curve = interp.Curve();

	gp_Pnt tempPt;
	gp_Vec endVec;
	Standard_Real U2 = curve->LastParameter();
	curve-> D1(U2,tempPt,endVec);
	dirScrollEnd = gp_Dir(endVec);

	return true;
}




bool OCCVolute::makeExitPipeTransition(bool firstRun, bool overRideTransPos, bool isClockWise,std::vector<asymSection> asymSects,double throat_area,double exit_area, double exit_length,
						bool exit_straight, 
						double exit_pipe_angle_v, double exit_pipe_angle_h, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition, double transitionLength, 
                        double exit_pipe_aspect1, double firstTransPos, int transSectOption, int firstTransPosOption, int& runFlag, double& calcedSecondTransPos,
						double earlyTurnAngle)
{
	bool success = false;
	exitTranSectWires.resize(4);
	exitTransSpineWires.resize(3);
	exitTranSectWires[0] = scrollSectWires[scrollSectWires.size()-1];

	success = makeExitStartSection(asymSects[0],asymSects[asymSects.size()-1],earlyTurnAngle);
	if(!success) 
	
		return false;

	success = makeExitPipeSpine(isClockWise, exit_length, exit_straight, exit_pipe_angle_v, exit_pipe_angle_h, exit_curve_order, pX_exitCrv, pY_exitCrv, exit_enable_transition, transitionLength, firstTransPos,
		firstTransPosOption, earlyTurnAngle);
	if(!success)
		return false;

//	success = makeTransitionMiddleSection(asymSects[0],asymSects[asymSects.size()-1], isClockWise, throat_area, transitionLength, exit_area, exit_pipe_aspect1, exit_length, firstTransPos, transSectOption);

	success = makeTransitionXsections(isClockWise, throat_area,exit_enable_transition, transitionLength, exit_area, exit_pipe_aspect1, exit_length);
	if(!success)
		return false;
 
 	success = sweepTransitionPipe(exit_straight,exit_enable_transition);
 
	return success;
}

bool OCCVolute::makeExitStartSection(asymSection scrollStartSect, asymSection scrollEndSect, double earlyTurnAngle)
{
	//make a curve from start section
	Handle_TColgp_HArray1OfPnt pSec1Points = new TColgp_HArray1OfPnt(1,scrollStartSect.npt);

	int i=0;
	int index = 1;
	for(i=0; i< scrollStartSect.npt; i++)
	{
		gp_Pnt tempPt(scrollStartSect.x[i],scrollStartSect.y[i],scrollStartSect.z[i]);
		pSec1Points->SetValue(index,tempPt);
		index++;
	}
			
					
	GeomAPI_Interpolate interpSec1(pSec1Points,Standard_False,1.0e-6);
	interpSec1.Perform();
	Handle_Geom_BSplineCurve pSec1Curve;
	if(interpSec1.IsDone())
		pSec1Curve = interpSec1.Curve();

	Standard_Real U1,U2;
	U1 = pSec1Curve->FirstParameter();
	U2 = pSec1Curve->LastParameter();

	gp_Pnt Pnt1;
	//check for symmetry
 	gp_Vec Vst,Vend;
	gp_Vec VecX = gp_Vec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
 	pSec1Curve->D1(U1,Pnt1,Vst);
	pSec1Curve->D1(U2,Pnt1,Vend);
 	double Ang1 = fabs(Vst.Angle(VecX))*180/M_PI;
	if(Ang1 > 90.0)
		Ang1 = 180.0-Ang1;
	double Ang2 = fabs(Vend.Angle(VecX))*180/M_PI;
	if(Ang2 > 90.0)
		Ang2 = 180.0-Ang2;


	Standard_Boolean isSymmetrical = fabs(Ang1-Ang2) < 1.0;
	
	////start from middle and go towards each end to find the U param of straight line segment
	Standard_Real Ubreak1;

	double angPrev = 0.0;
	gp_Vec VecPrev(gp_Pnt(0.0,0.0,0.0),gp_Pnt(0.0,0.0,1.0));
	for(i=100; i> 0; i--)
	{
		Standard_Real Uint = (U1+0.4*(U2-U1)) -(100-i)*(U2-U1)/100.0;
		gp_Vec V2;
		pSec1Curve->D1(Uint,Pnt1,V2);
		
		double ang = fabs(VecPrev.Angle(V2));
		if(i < 95 && ang < angPrev)
		{
			Ubreak1 = Uint;
			break;
		}
		else
		{
			angPrev = ang;
			VecPrev = V2;
		}
	}

	Standard_Real Ubreak2;

	if(isSymmetrical)
		Ubreak2 = U2-(Ubreak1-U1);
	else
	{
		for(i=0; i<100; i++)
		{
			Standard_Real Uint = (U1+0.6*(U2-U1)) + i*(U2-U1)/100.0;
			gp_Vec V2;
			pSec1Curve->D1(Uint,Pnt1,V2);
			
			double ang = fabs(VecPrev.Angle(V2));
			if(i > 5 && ang < angPrev)
			{
				Ubreak2 = Uint;
				break;
			}
			else
			{
				angPrev = ang;
				VecPrev = V2;
			}
		}
	}


	//sometimes the Ubreak1 or Ubreak2 seems to locate points too early. therefore get the max of both
	pSec1Curve->Segment(Ubreak1,Ubreak2);
/*	double gap = max((Ubreak1-U1),(U2-Ubreak2));
	pSec1Curve->Segment(U1+gap,U2-gap);
*/	

	//make the curve from end section
	Handle_TColgp_HArray1OfPnt pSecLastPoints = new TColgp_HArray1OfPnt(1,scrollEndSect.npt);
	i=0;
	index = 1;
	for(i=0; i< scrollEndSect.npt; i++)
	{
		gp_Pnt tempPt(scrollEndSect.x[i],scrollEndSect.y[i],scrollEndSect.z[i]);
		pSecLastPoints->SetValue(index,tempPt);
		index++;
	}
			
					
	GeomAPI_Interpolate interpSecLast(pSecLastPoints,Standard_False,1.0e-6);
	interpSecLast.Perform();
	Handle_Geom_BSplineCurve pSecLastCurve;
	if(interpSecLast.IsDone())
		pSecLastCurve = interpSecLast.Curve();

	ShapeAnalysis_Curve SAC;
	gp_Pnt projPnt1,projPnt2;
	Standard_Real Uparam1,Uparam2;
	double dist = SAC.Project(pSecLastCurve,pSec1Curve->StartPoint(),1.0e-8,projPnt1,Uparam1);
	dist = SAC.Project(pSecLastCurve,pSec1Curve->EndPoint(),1.0e-8,projPnt2,Uparam2);

	Standard_Real Ua = min(Uparam1,Uparam2);
	Standard_Real Ub = std::max(Uparam1,Uparam2);
	Standard_Real Ugap = (Ub-Ua)/10.0;

	if(earlyTurnAngle > 1.0e-3)
	{
		Ua = pSecLastCurve->FirstParameter();
		Ub = pSecLastCurve->LastParameter();
	}
	else
	{
		try
		{
			pSecLastCurve->Segment(Ua,Ub);
		}catch(...)
		{
			//this happens verry rarly
			//have to change Ua or Ub 
			if(Ua < pSecLastCurve->FirstParameter()) Ua = pSecLastCurve->FirstParameter();
			if(Ub > pSecLastCurve->LastParameter()) Ub = pSecLastCurve->LastParameter();
			pSecLastCurve->Segment(Ua,Ub);
		}
	}
	
	//re assign the exact start and end point on pSecLastCurve
	projPnt1 = pSecLastCurve->StartPoint();
	projPnt2 = pSecLastCurve->EndPoint();

	//now make the bottom segment (we have to pick one of the two choices below)

	//choice 1 (circular arc 1 which goes tangent to the top segment)
	gp_Vec dirVec;
	gp_Pnt tempPt;
	pSecLastCurve->D1(Ua,tempPt,dirVec);
	dirVec.Reverse();
	GC_MakeArcOfCircle mkCirc(projPnt1,dirVec,projPnt2);
	
	//get the mid point of the bottom arc
	Handle_Geom_Curve pbottomArc = mkCirc.Value();
	gp_Pnt PtmidArc = pbottomArc->Value(0.5*(pbottomArc->FirstParameter()+pbottomArc->LastParameter()));


	//choice 2 (circular arc 2 which goes tangent to the bottom level)
	gp_Pnt Pt1bottom = gp_Pnt(scrollEndSect.x[0],scrollEndSect.y[0],scrollEndSect.z[0]);
	gp_Pnt Pt2bottom = gp_Pnt(scrollEndSect.x[scrollEndSect.npt-1],scrollEndSect.y[scrollEndSect.npt-1],scrollEndSect.z[scrollEndSect.npt-1]);
	BRepBuilderAPI_MakeEdge mkEdge(Pt1bottom,Pt2bottom);
	TopoDS_Edge bottomEdge = mkEdge.Edge();
	Standard_Real U1bottom,U2bottom;
	Handle_Geom_Curve pbottomCrv = BRep_Tool::Curve(bottomEdge,U1bottom,U2bottom);
	gp_Pnt Ptmidbottom = pbottomCrv->Value(0.5*(U1bottom+U2bottom));
	
	Handle_Geom_Curve pSelectCrv;
	TopoDS_Edge SelectEdge;
	
	

//	BRepBuilderAPI_MakeWire mkW(BRepBuilderAPI_MakeEdge(pSecLastCurve).Edge());
/* 	bool straightBottom = true;
	TopoDS_Edge btmEdge;
	if(straightBottom)
	{
		btmEdge = BRepBuilderAPI_MakeEdge(projPnt1,projPnt2).Edge();
	}  
*/
	if(earlyTurnAngle > 1.0e-3)
	{
		SelectEdge = BRepBuilderAPI_MakeEdge(projPnt1,projPnt2).Edge();
	}
	else if((Ptmidbottom.X()*Ptmidbottom.X()+Ptmidbottom.Y()*Ptmidbottom.Y()) > PtmidArc.X()*PtmidArc.X()+PtmidArc.Y()*PtmidArc.Y() && earlyTurnAngle < 1.0e-3 ) //arc is too low
	{
		GC_MakeArcOfCircle mkCirc2(projPnt1,Ptmidbottom,projPnt2);
//		mkW.Add(BRepBuilderAPI_MakeEdge(mkCirc2.Value()).Edge());
		pSelectCrv = mkCirc2.Value();
		SelectEdge = BRepBuilderAPI_MakeEdge(pSelectCrv).Edge();

	}
	else
	{
//		mkW.Add(BRepBuilderAPI_MakeEdge(mkCirc.Value()).Edge());
		pSelectCrv = mkCirc.Value();
		SelectEdge = BRepBuilderAPI_MakeEdge(pSelectCrv).Edge();
	}

	//take parts off of the top curve to add to the bottom curve and make one single bottom curve
	Handle_Geom_BSplineCurve pPartACrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy());  
	Handle_Geom_BSplineCurve pPartBCrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy());     
	Handle_Geom_BSplineCurve pPartCCrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy()); 

	pPartACrv->Segment(Ua,Ua+Ugap);
	pPartBCrv->Segment(Ub-Ugap,Ub);
	pPartCCrv->Segment(Ua+Ugap,Ub-Ugap); 

	BRepBuilderAPI_MakeWire mkBottomWire(SelectEdge/*BRepBuilderAPI_MakeEdge(pSelectCrv).Edge()*/);
	mkBottomWire.Add(BRepBuilderAPI_MakeEdge(pPartACrv).Edge());
	Standard_Boolean isdone = mkBottomWire.IsDone();
	mkBottomWire.Add(BRepBuilderAPI_MakeEdge(pPartBCrv).Edge());
	isdone = mkBottomWire.IsDone();
	if(!isdone)
	{
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSelectCrv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pPartACrv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pPartBCrv).Edge());
		return false;
	}

	TopoDS_Wire tempWire = mkBottomWire.Wire();
	Handle_Geom_BSplineCurve pBottomCrvFinal = makeSingleCurve(tempWire,10);

	BRepBuilderAPI_MakeWire mkStartSectWire(BRepBuilderAPI_MakeEdge(pBottomCrvFinal).Edge());
	mkStartSectWire.Add(BRepBuilderAPI_MakeEdge(pPartCCrv).Edge());


	exitTranSectWires[0] = mkStartSectWire.Wire();

	return true;
}



bool OCCVolute::makeExitPipeSpine(bool IsClockWise,double exit_length, bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h,
									  int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv,
									  bool exit_enable_transition, double transitionLength, double firstTransPos, int firstTransPosOption,
									  double earlyTurnAngle)
{
	 //make the spine of exit pipe starting from the centroid of start section
	double pipe0_length = exit_length;
    if(exit_enable_transition)
		pipe0_length = exit_length*transitionLength/100.0;
    TopoDS_Edge spineEdge;
    gp_Pnt startSectionCentroid;
    if(exit_straight)
    {
        GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(exitTranSectWires[0]);
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();

        //find the normal to the stat section to get the spine direction
        BRepAdaptor_Surface tempSurf(mkFace.Face());
        BRepLProp_SLProps slProps(tempSurf,0.0,0.0,1,1.0e-6);
        gp_Dir surfaceNormal = slProps.Normal();
        gp_Vec normalVec(surfaceNormal);
        normalVec.Normalize();
		
        //test to see the normal vector should be in the +ve direction
        double angleBetween =dirScrollEnd.Angle(gp_Dir(normalVec));
        if(angleBetween > M_PI/2.0) // the normalVec direction has to be reversed
			normalVec.Multiply(-1.0);

		//in plane and out of plane angles
		if(exit_pipe_angle_v != 0.0)
		{
			double angInPlane = -exit_pipe_angle_v*M_PI/180.0;
			gp_Trsf  Rot;
			gp_Ax1 Ax1 = gp::OZ();
			gp_Trsf T;
			T.SetTranslation(gp_Pnt(0,0,0),startSectionCentroid);
			Ax1 = Ax1.Transformed(T);
			normalVec.Rotate(Ax1,angInPlane);
		}

		if(exit_pipe_angle_h != 0.0)
		{
			double angOutPlane = -exit_pipe_angle_h*M_PI/180.0;
			gp_Trsf  Rot;
			gp_Ax1 Ax1 = gp::OX();
			gp_Trsf T;
			T.SetTranslation(gp_Pnt(0,0,0),startSectionCentroid);
			Ax1 = Ax1.Transformed(T);
			normalVec.Rotate(Ax1,angOutPlane);
		}

		//the transition region is broken in to two segments by fraction
		double fraction = 1.0;
/*		if(firstTransPosOption == 0) //program calculated
		{
			//take the scroll cirumference
			double radius = sqrt(startSectionCentroid.X()*startSectionCentroid.X()+startSectionCentroid.Y()*startSectionCentroid.Y());
			double circum = 2.0*M_PI*radius;
			double firstSegLength = circum/36.0;
			fraction = firstSegLength/pipe0_length;
		}
		else //user specified
			fraction = firstTransPos/100.0;
*/
		gp_Vec vec1 = normalVec;
		double seg1Length = fraction*pipe0_length;
		vec1.Multiply(seg1Length);
		gp_Pnt end = startSectionCentroid;
		end.Translate(vec1);

		

		BRepBuilderAPI_MakeEdge mkEdge1(startSectionCentroid,end);
        spineEdge = mkEdge1.Edge();
		BRepBuilderAPI_MakeWire spineWire(spineEdge);
		exitTransSpineWires[0] = spineWire.Wire();
/*
		gp_Vec vec2 = normalVec;
		double seg2Length = (1.0-fraction)*pipe0_length;
		vec2.Multiply(seg2Length);
		gp_Pnt end2 = end;
		end2.Translate(vec2);

		BRepBuilderAPI_MakeEdge mkEdge2(end,end2);
        spineEdge = mkEdge2.Edge();
		BRepBuilderAPI_MakeWire spineWire2(spineEdge);
		exitTransSpineWires[1] = spineWire2.Wire();
*/
		//exit pipe end segment spine
		gp_Vec vec3 = normalVec;
		double seg3Length = exit_length-pipe0_length;
		vec3.Multiply(seg3Length);
		gp_Pnt end3 = end;//end2;
		end3.Translate(vec3);

		BRepBuilderAPI_MakeEdge mkEdge3(end,end3);
        spineEdge = mkEdge3.Edge();
		BRepBuilderAPI_MakeWire spineWire3(spineEdge);
		exitTransSpineWires[1] = spineWire3.Wire();

    }
    else
    {
        TColgp_Array1OfPnt CurvePoles(1,exit_curve_order);
        for(int i=1;i<=exit_curve_order; i++)
        {
            gp_Pnt spinePt(pX_exitCrv[i-1], pY_exitCrv[i-1],0.0);
            CurvePoles.SetValue(i,spinePt); 
        }
        Handle(Geom_BezierCurve) pexitSpineCurve = new Geom_BezierCurve(CurvePoles);

		Handle_Geom_BezierCurve pBezCrv5 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
				 
        GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(exitTranSectWires[0]);
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();
		
		Standard_Real umin, umax, vmin, vmax;
		BRepTools::UVBounds(mkFace.Face(),umin, umax, vmin, vmax);
		Handle(Geom_Surface) aSurface = BRep_Tool::Surface(mkFace.Face());
		GeomLProp_SLProps props(aSurface, umin, vmin,1, 0.01);
		gp_Dir normal = props.Normal();
		//make sure the normal is out of scroll exit 
		TopoDS_Face tempFace;
		for(TopExp_Explorer explr(scrollShape1,TopAbs_FACE);explr.More();explr.Next())
		{
			tempFace = TopoDS::Face(explr.Current());
			break;
		}
		Handle_Geom_Surface pScrollShape1Surf = BRep_Tool::Surface(tempFace);
		Standard_Real Us1,Us2,Vs1,Vs2;
		pScrollShape1Surf->Bounds(Us1,Us2,Vs1,Vs2);
		Handle_Geom_Curve ptempCrv = pScrollShape1Surf->UIso(Us1);
		gp_Pnt tempPt;
		gp_Vec tempVec;
		ptempCrv->D1(Vs1,tempPt,tempVec);
		gp_Dir tempDir(tempVec);
		if(fabs(tempVec.Angle(normal)) > M_PI/2.0)
			normal.Reverse();

		//move the spine curve to the centroid of the start section
        gp_Trsf T;
        gp_Pnt P1;
        gp_Vec V1;
        pexitSpineCurve->D1(0.0,P1,V1);//StartPoint();
        T.SetTranslation(P1,startSectionCentroid);
		pexitSpineCurve->Transform(T);

		Handle_Geom_BezierCurve pBezCrv = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		
		//align the spine curve with swrirl section spine curve if the curve is towards opposite direction		
		gp_Trsf  Rot; 
		//gp_Vec V2(dirScrollEnd);
		gp_Vec V2(normal);
		Standard_Real rotAngle = V2.Angle(V1);
		if(IsClockWise || earlyTurnAngle > 1.0e-2)
			rotAngle = -rotAngle;
		if(true)
		{
			gp_Ax1 rotAxis = gp::OZ();
			T.SetTranslation(gp_Pnt(0,0,0),startSectionCentroid);
			rotAxis = rotAxis.Transformed(T);
			Rot.SetRotation(rotAxis,rotAngle);
			pexitSpineCurve->Transform(Rot);

			//second rotation around scroll exit axis
			if(!IsClockWise)
			{
				rotAxis = gp_Ax1(startSectionCentroid,normal);
				rotAngle = M_PI;
				Rot.SetRotation(rotAxis,rotAngle);
				pexitSpineCurve->Transform(Rot);
			}
		}
/*		if(rotAngle > M_PI/2.0)
		{
			//first rotate perpendicular to OZ direction if necessary
			rotAngle = M_PI;
			gp_Pnt P2(0,0,startSectionCentroid.Z());
			gp_Ax1 rotAxis(P2,gp_Dir(gp_Vec(P2,startSectionCentroid)));		
			Rot.SetRotation(rotAxis,rotAngle);
			pexitSpineCurve->Transform(Rot);
		}
*/
		gce_MakeLin mkLin(startSectionCentroid,normal);
		BRepBuilderAPI_MakeEdge mke(mkLin.Value(),0.0,0.1);
		
		Handle_Geom_BezierCurve pBezCrv2 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		
		//now rotate around OZ direction
/*		gp_Ax1 rotAxis = gp::OZ();
		T.SetTranslation(gp_Pnt(0,0,0),startSectionCentroid);
		rotAxis = rotAxis.Transformed(T);
		pexitSpineCurve->D1(0.0,P1,V1);
		rotAngle = V2.Angle(V1);
		if(IsClockWise)
			Rot.SetRotation(rotAxis,-1.0*rotAngle);
		else
			Rot.SetRotation(rotAxis,rotAngle);
		pexitSpineCurve->Transform(Rot);
*/		

		//get the requierd length of this curve
		double fraction = 1.0;
/*		if(firstTransPosOption == 0) //default
		{
			//take the scroll cirumference
			double radius = sqrt(startSectionCentroid.X()*startSectionCentroid.X()+startSectionCentroid.Y()*startSectionCentroid.Y());
			double circum = 2.0*M_PI*radius;
			double firstSegLength = circum/36.0;
			fraction = firstSegLength/pipe0_length;
		}
		else
			fraction = firstTransPos/100.0;
*/
		double seg1Length = fraction*pipe0_length;

		GeomAdaptor_Curve GAC(pexitSpineCurve);
		Standard_Real Uend = 0.0;
		Standard_Real U1 = GAC.FirstParameter();
		Standard_Real U2 = GAC.LastParameter();
		GCPnts_AbscissaPoint GCPA(GAC,seg1Length,U1);
		if(GCPA.IsDone())
			Uend = GCPA.Parameter();

	/*	Standard_Real Uend2 = 0.0;
		GCPnts_AbscissaPoint GCPA2(GAC,pipe0_length,U1);
		if(GCPA2.IsDone())
			Uend2 = GCPA2.Parameter();
*/


		Handle_Geom_BezierCurve pexitPipeSpine_Seg2 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		pexitPipeSpine_Seg2->Segment(Uend,U2);
		exitTransSpineWires[1] = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pexitPipeSpine_Seg2).Edge()).Wire();


		Handle_Geom_BezierCurve pexitPipeSpine_Seg1 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		pexitPipeSpine_Seg1->Segment(U1,Uend);
		exitTransSpineWires[0] = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pexitPipeSpine_Seg1).Edge()).Wire();

/*		Handle_Geom_BezierCurve pexitPipeSpine_Seg3 = Handle(Geom_BezierCurve)::DownCast(pexitSpineCurve->Copy());
		pexitPipeSpine_Seg3->Segment(Uend2,U2);
		exitTransSpineWires[2] = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pexitPipeSpine_Seg3).Edge()).Wire(); */
    }


	return true;
}

bool OCCVolute::makeTransitionXsections(bool isClockWise, double throat_area, bool exit_enable_transition, double transitionLength, double exit_area, double exit_pipe_aspect1,
											double exit_length)
{
	 
    TopoDS_Wire endSectionWire;
    if(exit_enable_transition)
    {
        double targetArea = throat_area+transitionLength*(exit_area-throat_area)/100.0;
        gp_Pnt endSectionCentroid;
        gp_Vec Vec1;
        Standard_Real U1,U2;

		TopExp_Explorer Ex;
		TopoDS_Edge spineEdge;
		for (Ex.Init(exitTransSpineWires[0],TopAbs_EDGE); Ex.More(); Ex.Next())
		{ 
			spineEdge = TopoDS::Edge(Ex.Current());
		}

        BRep_Tool::Range(spineEdge,U1,U2);
        Handle_Geom_Curve spineCurve = BRep_Tool::Curve(spineEdge,U1,U2);
        spineCurve->D1(U2,endSectionCentroid,Vec1);

        //update the end direction of wrap spine
       // dirTransEnd = gp_Dir(Vec1);

        gp_Dir dir(Vec1);
        gp_Ax2 majorAxis(endSectionCentroid,dir);
		gp_Dir Vy;
		if(isClockWise)
		{
			Vy = gp_Dir(gp_Vec(0,0,1));
			if(exit_pipe_aspect1 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}
		else
		{
			Vy = gp_Dir(gp_Vec(0,0,-1));
			if(exit_pipe_aspect1 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}

        double majR = 0.0;
        double minR = 0.0;
        bool found = getTargetEllipse(targetArea, exit_pipe_aspect1, majR, minR);

        //we have the final ellipse shape now
        Geom_Ellipse E_converged(majorAxis,majR,minR);       
        gp_Elips ellipse = E_converged.Elips();

		Handle_Geom_BSplineCurve bspEllipse;
		if(isClockWise)
		{
			if(exit_pipe_aspect1 >= 1.0) 
				bspEllipse = resetEllipsetoCurve(E_converged,-0.25/*-1.0*/); 
			else
				bspEllipse = resetEllipsetoCurve(E_converged,0.0/*-0.25*/); 
		}
		else
		{
			if(exit_pipe_aspect1 >= 1.0) 
				bspEllipse = resetEllipsetoCurve(E_converged,0.25/*0.0*/ ); 
			else
				bspEllipse = resetEllipsetoCurve(E_converged,0.0/*0.75*/);  
		}

        BRepBuilderAPI_MakeEdge mkEdge(bspEllipse);
        BRepBuilderAPI_MakeWire mkWire(mkEdge);
        endSectionWire = mkWire.Wire(); 

		exitTranSectWires[1] = endSectionWire;

    }




	return true;
}

bool OCCVolute::sweepTransitionPipe(bool exit_straight, bool exit_enable_transition)
{
	 //make the exit pipe sweep
    // if exit pipe is straight use the ThruSections algorithm whic is more robust
    // if not straight, use the makePipeShell algorithm which will sweep in a non linear path

    if(!exit_straight || !exit_enable_transition) 
    {
        BRepOffsetAPI_MakePipeShell exitPipeShell(exitTransSpineWires[0]);
        exitPipeShell.Add(exitTranSectWires[0],Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
            exitPipeShell.Add(exitTranSectWires[1],Standard_False,Standard_False) ;
   
		try
		{
			exitPipeShell.Build();
		//	Standard_Boolean madeSolid = exitPipeShell.MakeSolid();
		//	exitPipeShell.MakeSolid();
		}
		catch(...)
		{
			return false;
		}
        
        transPipeShape1 = exitPipeShell.Shape();

/*
		BRepOffsetAPI_MakePipeShell exitPipeShell2(exitTransSpineWires[1]);
        exitPipeShell2.Add(exitTranSectWires[1],Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
            exitPipeShell2.Add(exitTranSectWires[2],Standard_False,Standard_False) ;
       
		try
		{
			exitPipeShell2.Build();
		//	exitPipeShell2.MakeSolid();
		}
		catch(...)
		{
			return false;
		}
        
        transPipeShape2 = exitPipeShell2.Shape(); */
    }
    else
    {
        BRepOffsetAPI_ThruSections thruSections1(Standard_False,Standard_False,1.0e-6);
        thruSections1.AddWire(exitTranSectWires[0]);
        thruSections1.AddWire(exitTranSectWires[1]);
		try
		{
			thruSections1.Build();
			transPipeShape1 = thruSections1.Shape();
		}
		catch(...) 
		{
			return false;
		}      
 /*
		BRepOffsetAPI_ThruSections thruSections2(Standard_False,Standard_False,1.0e-6);
        thruSections2.AddWire(exitTranSectWires[1]);
        thruSections2.AddWire(exitTranSectWires[2]);
		try
		{
			thruSections2.Build();
			transPipeShape2 = thruSections2.Shape(); 
		}
		catch(...)
		{
			return false;
		} */
    }

	return true;
}


bool OCCVolute::makeExitPipe(bool isClockWise,double exit_area, bool exit_enable_transition, double exit_pipe_aspect2, bool exit_straight, bool exit_enable_extension,
							 double exit_extension_length, double exit_extension_reduction)
{
	bool success = false;
	makeExitEndSection(isClockWise,exit_enable_transition,exit_pipe_aspect2, exit_area);

	if(exit_enable_extension)
		makeExitExtensionSection(exit_extension_length, exit_extension_reduction);

	success = sweepExitPipe(exit_straight,exit_enable_transition, exit_enable_extension);

	return success;
}



bool OCCVolute::makeExitEndSection(bool isClockWise, bool exit_enable_transition, double exit_pipe_aspect2, double exit_area)
{
	//if transition is requested make the ellipse shape
    TopoDS_Wire endSectionWire;
    if(exit_enable_transition)
    {
		//exit spine segment
		
		TopExp_Explorer Ex;
		TopoDS_Edge spineEdge;
		for (Ex.Init(exitTransSpineWires[1],TopAbs_EDGE); Ex.More(); Ex.Next())
		{
			spineEdge = TopoDS::Edge(Ex.Current());
		}

        gp_Pnt endSectioncentroid;
        gp_Vec V1;
        Standard_Real U1,U2;
        Handle_Geom_Curve spineCurve = BRep_Tool::Curve(spineEdge,U1,U2);
        spineCurve->D1(U2,endSectioncentroid,V1);

        gp_Dir dir(V1);
        gp_Ax2 majorAxis(endSectioncentroid,dir);
		gp_Dir Vy;
		if(isClockWise)
		{
			Vy = gp_Dir(gp_Vec(0,0,1));
			if(exit_pipe_aspect2 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}
		else
		{
			Vy = gp_Dir(gp_Vec(0,0,-1));
			if(exit_pipe_aspect2 >= 1.0) 
				majorAxis.SetYDirection(Vy);
			else
				majorAxis.SetXDirection(Vy);
		}

        
        double majR = 0.0;
        double minR = 0.0;
        bool found = getTargetEllipse(exit_area, exit_pipe_aspect2, majR, minR);

        //we have the final ellipse shape now
        Geom_Ellipse E_converged(majorAxis,majR,minR);
        gp_Elips ellipse = E_converged.Elips();

		Handle_Geom_BSplineCurve bspEllipse;
		if(isClockWise)
		{
			if(exit_pipe_aspect2 >= 1.0) 
				bspEllipse = resetEllipsetoCurve(E_converged,-0.25/*-0.25*//*-1.0*/); 
			else
				bspEllipse = resetEllipsetoCurve(E_converged,0.0/*-0.25*/); 
		}
		else
		{
			if(exit_pipe_aspect2 >= 1.0) 
				bspEllipse = resetEllipsetoCurve(E_converged,0.25/*0.0*/); 
			else
				bspEllipse = resetEllipsetoCurve(E_converged,0.0/*0.75*/);  
		}

        BRepBuilderAPI_MakeEdge mkEdge(bspEllipse/*ellipse*/);
        BRepBuilderAPI_MakeWire mkWire(mkEdge);
        exitTranSectWires[2] = mkWire.Wire();

    }


	return true;
}

bool OCCVolute::makeExitExtensionSection(double exit_extension_length, double exit_extension_reduction)
{
	TopoDS_Wire extSectionEndWire = exitTranSectWires[2]; //start with exit end X section;
	
	TopoDS_Edge spineEdge;
	for (TopExp_Explorer Ex(exitTransSpineWires[1],TopAbs_EDGE); Ex.More(); Ex.Next())
		spineEdge = TopoDS::Edge(Ex.Current());

    gp_Pnt Pt1,Pt2;
    gp_Vec V1;
    Standard_Real U1,U2;
    Handle_Geom_Curve spineCurve = BRep_Tool::Curve(spineEdge,U1,U2);
    spineCurve->D1(U2,Pt1,V1);

	V1.Normalize();
	V1.Multiply(exit_extension_length);
	Pt2 = Pt1;
	Pt2.Translate(V1); 



	gp_Trsf  t1;
    t1.SetTranslation(Pt1,Pt2); 
    BRepBuilderAPI_Transform translate(extSectionEndWire,t1); 

	exitTranSectWires[3] = TopoDS::Wire(translate.Shape());

	//scale if necessary
	gp_Trsf  t2;
	t2.SetScale(Pt2,exit_extension_reduction/100.0);
	BRepBuilderAPI_Transform translate2(exitTranSectWires[3],t2);
	exitTranSectWires[3] = TopoDS::Wire(translate2.Shape());
	
	return true;
}


bool OCCVolute::sweepExitPipe(bool exit_straight, bool exit_enable_transition, bool exit_enable_extension)
{
	//make the exit pipe sweep
    // if exit pipe is straight use the ThruSections algorithm whic is more robust
    // if not straight, use the makePipeShell algorithm which will sweep in a non linear path

    if(!exit_straight || !exit_enable_transition) 
    {
		//////////////////
 /*		TopExp_Explorer Ex;
		TopoDS_Edge tempEdge;
		for (Ex.Init(exitTranSectWires[1],TopAbs_EDGE); Ex.More(); Ex.Next())
		{ 
			tempEdge = TopoDS::Edge(Ex.Current());
		}
		Standard_Real U1,U2;
		BRep_Tool::Range(tempEdge,U1,U2);
		Handle_Geom_Curve pCrv1 = BRep_Tool::Curve(tempEdge,U1,U2);
		Handle_Geom_BSplineCurve bsp1a = Handle(Geom_BSplineCurve)::DownCast(pCrv1->Copy());
		Handle_Geom_BSplineCurve bsp1b = Handle(Geom_BSplineCurve)::DownCast(pCrv1->Copy());
		Handle_Geom_BSplineCurve bsp1c = Handle(Geom_BSplineCurve)::DownCast(pCrv1->Copy());
		bsp1a->Segment(U1,U1+0.25*(U2-U1));
		bsp1b->Segment(U1+0.25*(U2-U1),U1+0.75*(U2-U1));
		bsp1c->Segment(U1+0.75*(U2-U1),U2);


		for (Ex.Init(exitTranSectWires[2],TopAbs_EDGE); Ex.More(); Ex.Next())
		{ 
			tempEdge = TopoDS::Edge(Ex.Current());
		}
		
		BRep_Tool::Range(tempEdge,U1,U2);
		Handle_Geom_Curve pCrv2 = BRep_Tool::Curve(tempEdge,U1,U2);
		Handle_Geom_BSplineCurve bsp2a = Handle(Geom_BSplineCurve)::DownCast(pCrv2->Copy());
		Handle_Geom_BSplineCurve bsp2b = Handle(Geom_BSplineCurve)::DownCast(pCrv2->Copy());
		Handle_Geom_BSplineCurve bsp2c = Handle(Geom_BSplineCurve)::DownCast(pCrv2->Copy());
		bsp2a->Segment(U1,U1+0.25*(U2-U1));
		bsp2b->Segment(U1+0.25*(U2-U1),U1+0.75*(U2-U1));
		bsp2c->Segment(U1+0.75*(U2-U1),U2);

		//BRepBuilderAPI_MakeWire mkW1(BRepBuilderAPI_MakeEdge(bsp1b).Edge());
		TopoDS_Edge tEdge1 = BRepBuilderAPI_MakeEdge(bsp1a).Edge();
		TopoDS_Edge tEdge2 = BRepBuilderAPI_MakeEdge(bsp1b).Edge();
		TopoDS_Edge tEdge3 = BRepBuilderAPI_MakeEdge(bsp1c).Edge();

		BRepBuilderAPI_MakeWire mkW1(tEdge1,tEdge2,tEdge3);
	//	mkW1.Add(BRepBuilderAPI_MakeEdge(bsp1a).Edge());
		TopoDS_Wire wire1 = mkW1.Wire();

	//	BRepBuilderAPI_MakeWire mkW2(BRepBuilderAPI_MakeEdge(bsp2b).Edge());
		tEdge1 = BRepBuilderAPI_MakeEdge(bsp2a).Edge();
		tEdge2 = BRepBuilderAPI_MakeEdge(bsp2b).Edge();
		tEdge3 = BRepBuilderAPI_MakeEdge(bsp2c).Edge();
		BRepBuilderAPI_MakeWire mkW2(tEdge1,tEdge2,tEdge3);
	//	mkW2.Add(BRepBuilderAPI_MakeEdge(bsp2a).Edge());
		TopoDS_Wire wire2 = mkW2.Wire();
*/

		BRepOffsetAPI_MakePipeShell exitPipeShell(exitTransSpineWires[1]);
        exitPipeShell.Add(exitTranSectWires[1],Standard_False,Standard_True) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
			exitPipeShell.Add(exitTranSectWires[2],Standard_False,Standard_True) ;
	
  		try
		{
			exitPipeShell.Build();
 			exitPipeShell.MakeSolid();
			exitPipeShape1 = exitPipeShell.Shape(); 
		}
		catch(...)   
		{
			try
			{
				BRepOffsetAPI_MakePipeShell exitPipeShell2(exitTransSpineWires[1]);
				exitPipeShell2.Add(exitTranSectWires[1],Standard_False,Standard_True) ;
				
				exitPipeShell2.Build();
 				exitPipeShell2.MakeSolid();
				exitPipeShape1 = exitPipeShell2.Shape(); 
			}
			catch(...)
			{
				return false;
			} 
			
		}
         

		 //extract the exit pipe faces except the starting end cap
		TopExp_Explorer Ex;
		int numFaces = 0;
		for (Ex.Init(exitPipeShape1,TopAbs_FACE); Ex.More(); Ex.Next())
			numFaces++;

		int i=1;
		for (Ex.Init(exitPipeShape1,TopAbs_FACE); Ex.More(); Ex.Next())
		{ 
			if(i == numFaces-1)	// do not get the starting end cap
			{
			}
			else if(i == numFaces && exit_enable_extension)	//if there is an extension to the exit pipe do not get the ending endcap either
			{
			}
			else
			{
				exitPipeFaces.push_back(TopoDS::Face(Ex.Current()));	
			}
			i++;
		}

		/////////////////////
/*        BRepOffsetAPI_MakePipeShell exitPipeShell(exitTransSpineWires[2]);
        exitPipeShell.Add(exitTranSectWires[2],Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
			exitPipeShell.Add(exitTranSectWires[3],Standard_False,Standard_False) ;
	
        exitPipeShell.Build();
 		exitPipeShell.MakeSolid();
        
         exitPipeShape1 = exitPipeShell.Shape();   */
    }
    else
    {
		/////split the starting and ending X section curves in to 2 segments before making the exit pipe segment
		//// this is needed to do the boolean operation and face making of exit pipe later
		TopExp_Explorer Ex;
		TopoDS_Edge tempEdge;
		for (Ex.Init(exitTranSectWires[1],TopAbs_EDGE); Ex.More(); Ex.Next())
		{ 
			tempEdge = TopoDS::Edge(Ex.Current());
		}
		Standard_Real U1,U2;
		BRep_Tool::Range(tempEdge,U1,U2);
		Handle_Geom_Curve pCrv1 = BRep_Tool::Curve(tempEdge,U1,U2);
		Handle_Geom_BSplineCurve bsp1a = Handle(Geom_BSplineCurve)::DownCast(pCrv1->Copy());
		Handle_Geom_BSplineCurve bsp1b = Handle(Geom_BSplineCurve)::DownCast(pCrv1->Copy());
		bsp1a->Segment(U1,U1+0.5*(U2-U1));
		bsp1b->Segment(U1+0.5*(U2-U1),U2);


		for (Ex.Init(exitTranSectWires[2],TopAbs_EDGE); Ex.More(); Ex.Next())
		{ 
			tempEdge = TopoDS::Edge(Ex.Current());
		}
		
		BRep_Tool::Range(tempEdge,U1,U2);
		Handle_Geom_Curve pCrv2 = BRep_Tool::Curve(tempEdge,U1,U2);
		Handle_Geom_BSplineCurve bsp2a = Handle(Geom_BSplineCurve)::DownCast(pCrv2->Copy());
		Handle_Geom_BSplineCurve bsp2b = Handle(Geom_BSplineCurve)::DownCast(pCrv2->Copy());
		bsp2a->Segment(U1,U1+0.5*(U2-U1));
		bsp2b->Segment(U1+0.5*(U2-U1),U2);

		BRepBuilderAPI_MakeWire mkW1(BRepBuilderAPI_MakeEdge(bsp1b).Edge());
		mkW1.Add(BRepBuilderAPI_MakeEdge(bsp1a).Edge());
		TopoDS_Wire wire1 = mkW1.Wire();

		BRepBuilderAPI_MakeWire mkW2(BRepBuilderAPI_MakeEdge(bsp2b).Edge());
		mkW2.Add(BRepBuilderAPI_MakeEdge(bsp2a).Edge());
		TopoDS_Wire wire2 = mkW2.Wire();

        BRepOffsetAPI_ThruSections thruSections1(Standard_True,Standard_False,1.0e-6);
        thruSections1.AddWire(wire1/*exitTranSectWires[2]*/);
        thruSections1.AddWire(wire2/*exitTranSectWires[3]*/);
        thruSections1.Build();
        exitPipeShape1 = thruSections1.Shape();

		//extract the exit pipe faces except the starting end cap
		int numFaces = 0;
		for (Ex.Init(exitPipeShape1,TopAbs_FACE); Ex.More(); Ex.Next())
			numFaces++;

		int i=1;
		for (Ex.Init(exitPipeShape1,TopAbs_FACE); Ex.More(); Ex.Next())
		{ 
			if(i == numFaces-1)	// do not get the starting end cap
			{
			}
			else if(i == numFaces && exit_enable_extension)	//if there is an extension to the exit pipe do not get the ending endcap either
			{
			}
			else
			{
				exitPipeFaces.push_back(TopoDS::Face(Ex.Current()));
			}
			i++;
		}
    }

	/// make the exit extension if it is there
	if(exit_enable_extension)
	{
		BRepOffsetAPI_ThruSections thruSectExitExtension(Standard_True,Standard_True,1.0e-6);
        thruSectExitExtension.AddWire(exitTranSectWires[2]);
        thruSectExitExtension.AddWire(exitTranSectWires[3]);
        thruSectExitExtension.Build();
        exitPipeShape2 = thruSectExitExtension.Shape();

		//extract the exit pipe extension faces
		int i=0;
		for (TopExp_Explorer Ex(exitPipeShape2,TopAbs_FACE); Ex.More(); Ex.Next())
		{ 
			if(i!=1)	//do not get the starting end cap
				exitPipeFaces.push_back(TopoDS::Face(Ex.Current()));
			 i++;
		}
	}

	

	return true;
}




int OCCVolute::tongueBoolean(int filletMethod,bool exit_straight, double exit_length, double transitionLength, double firstTransPos, int firstTransPosOption, double earlyTurnAngle)
{
	int errorCode = 0;
	bool success = 0;

	///extract Scroll surface for boolean
	TopExp_Explorer Ex1;
	TopoDS_Face face1;
	int facenum = 0;
	int i=0;
    for (Ex1.Init(scrollShape1,TopAbs_FACE); Ex1.More(); Ex1.Next())
    {        
		if(i== facenum)
			face1 = TopoDS::Face(Ex1.Current());
		else
			tongueRegionScrollFaces.push_back(TopoDS::Face(Ex1.Current()));

		i++;
    }
	///extract Transition segment surface for boolean
    TopoDS_Face face2;
	facenum = 0;
	i=0;
    for (Ex1.Init(transPipeShape1,TopAbs_FACE); Ex1.More(); Ex1.Next())
    {       
		if(i==facenum)
			face2 = TopoDS::Face(Ex1.Current());
		else 
		{
			tongueRegionPipeFaces.push_back(TopoDS::Face(Ex1.Current()));
		}
		i++;
    }


	int flag = JoinTransitionPipe2B(filletMethod,face1,face2,exit_length, transitionLength, firstTransPos, firstTransPosOption, earlyTurnAngle);
	success = (flag == 0);
	if(!success)
	{
		if(flag == 1)
			errorCode = 51;
		else if(flag == 2)
			errorCode = 52;
		else
			errorCode = 5;
	}


	return errorCode;	
}


int OCCVolute::JoinTransitionPipe2()
{
	//extract geom faces
	TopExp_Explorer Ex1;
    TopoDS_Face face1;
	Standard_Real U1,U2,V1,V2;

	///Scroll surface
	int facenum = 0;
	int i=0;
    for (Ex1.Init(scrollShape1,TopAbs_FACE); Ex1.More(); Ex1.Next())
    {        
		if(i== facenum)
			face1 = TopoDS::Face(Ex1.Current());
		else
			tongueRegionScrollFaces.push_back(TopoDS::Face(Ex1.Current()));

		i++;
    }
    Handle_Geom_Surface pSurf1 = BRep_Tool::Surface(face1);

	pSurf1->Bounds(U1,U2,V1,V2);

	


	///Transition segment surface
    TopoDS_Face face2;
	Standard_Real U1a,U2a,V1a,V2a;

	facenum = 0;
	i=0;
    for (Ex1.Init(transPipeShape1,TopAbs_FACE); Ex1.More(); Ex1.Next())
    {       
		if(i==facenum)
			face2 = TopoDS::Face(Ex1.Current());
		else 
		{
			tongueRegionPipeFaces.push_back(TopoDS::Face(Ex1.Current()));
		}
		i++;
    }

    Handle_Geom_Surface pSurf2 = BRep_Tool::Surface(face2);
	pSurf2->Bounds(U1a,U2a,V1a,V2a);

 
	BRepAlgoAPI_Section makeBool1(face1,face2);
	//makeBool1.Build();	//this is not neccesary, automatically call from construction
	Standard_Boolean isDone = makeBool1.IsDone();
	TopoDS_Shape sh1;
	if(isDone)
		sh1 = makeBool1.Shape();
	else
		return false;


	//pick the intersection curve (select the longest edge if there is more than one
	TopoDS_Edge cutEdge;
	i=0;
	for (TopExp_Explorer ex(sh1,TopAbs_EDGE); ex.More(); ex.Next())
	{
		TopoDS_Edge aEdge = TopoDS::Edge(ex.Current());  
		if(i==0)
			cutEdge = aEdge;
		if(i > 0)
		{
			GProp_GProps LProps;
			BRepGProp::LinearProperties(aEdge,LProps);
			double length1 = LProps.Mass();

			BRepGProp::LinearProperties(cutEdge,LProps);
			double length2 = LProps.Mass();

			if(length1 > length2)
				cutEdge = aEdge;
		}
		i++;
	}

	Standard_Real Ut1,Ut2;
	BRep_Tool::Range(cutEdge,Ut1,Ut2);

	Handle_Geom_BSplineCurve ptempCrv = makeSingleCurve(BRepBuilderAPI_MakeWire(cutEdge).Wire());
	Ut1 = ptempCrv->FirstParameter();
	Ut2 = ptempCrv->LastParameter();


	cutEdge = BRepBuilderAPI_MakeEdge(ptempCrv).Edge();
	BRep_Tool::Range(cutEdge,Ut1,Ut2);

	

	//make the cut face on scroll
	Standard_Real Ua1,Ua2;
	Handle_Geom_Curve pCutCrv = BRep_Tool::Curve(cutEdge,Ua1,Ua2);
	gp_Pnt Pta1 = pCutCrv->Value(Ua1);
	gp_Pnt Pta2 = pCutCrv->Value(Ua2);

	ShapeAnalysis_Curve SAC;
	gp_Pnt projPnt1,projPnt2;
	Standard_Real Uparam1,Uparam2;
	double dist = SAC.Project(pSurf1->VIso(V1),Pta1,1.0e-5,projPnt1,Uparam1);
	dist = SAC.Project(pSurf1->VIso(V1),Pta2,1.0e-5,projPnt2,Uparam2);

	Handle_Geom_BSplineCurve pCrv1 = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(V1)); 
	pCrv1->Segment(U1,min(Uparam1,Uparam2));
	TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(pCrv1).Edge();

	Handle_Geom_BSplineCurve pCrv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(V1)); 
	pCrv2->Segment(std::max(Uparam1,Uparam2),U2);
	TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(pCrv2).Edge();

	BRepBuilderAPI_MakeWire mkWire( cutEdge );
	mkWire.Add(edge1);
	Standard_Boolean success = mkWire.IsDone();
	mkWire.Add(edge2);
	success = mkWire.IsDone();

	//add the balance edges to close the wire on scroll surface
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->UIso(U1)).Edge());
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->UIso(U2)).Edge());
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->VIso(V2)).Edge());

	BRepBuilderAPI_MakeFace mkFace1(pSurf1,mkWire.Wire());
	success = mkFace1.IsDone();

	if(success)
	{
		TopoDS_Face tempface = mkFace1.Face();
		ShapeFix_Face fixFace1(tempface);
		fixFace1.Perform();
		Standard_Boolean status = fixFace1.Status(ShapeExtend_OK);
		trimFace1 = fixFace1.Face();
	}
	else
		return false;


	//make the cut face on exit pipe trans segment
	dist = SAC.Project(pSurf2->VIso(V1a),Pta1,1.0e-5,projPnt1,Uparam1);
	dist = SAC.Project(pSurf2->VIso(V1a),Pta2,1.0e-5,projPnt2,Uparam2);

	Handle_Geom_BSplineCurve pCrv3 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv3->Segment(U1a,min(Uparam1,Uparam2));
	TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(pCrv3).Edge();

	Handle_Geom_BSplineCurve pCrv4 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv4->Segment(std::max(Uparam1,Uparam2),U2a);
	TopoDS_Edge edge4 = BRepBuilderAPI_MakeEdge(pCrv4).Edge();

	BRepBuilderAPI_MakeWire mkWire2(cutEdge);
	mkWire2.Add(edge3);
	success = mkWire2.IsDone();
	mkWire2.Add(edge4);
	success = mkWire2.IsDone();

	//add the balance edges to close the wire on scroll surface
	mkWire2.Add(BRepBuilderAPI_MakeEdge(pSurf2->UIso(U1a)).Edge());
	mkWire2.Add(BRepBuilderAPI_MakeEdge(pSurf2->UIso(U2a)).Edge());
	mkWire2.Add(BRepBuilderAPI_MakeEdge(pSurf2->VIso(V2a)).Edge());

	BRepBuilderAPI_MakeFace mkFace2(pSurf2,mkWire2.Wire());
	success = mkFace2.IsDone();

	if(success)
	{
		TopoDS_Face tempface = mkFace2.Face();
		ShapeFix_Face fixFace2(tempface);
		fixFace2.Perform();
		Standard_Boolean status = fixFace2.Status(ShapeExtend_OK);
		trimFace2 = fixFace2.Face();
	}
	else
		return false;
	
	return 0;
}


int OCCVolute::JoinTransitionPipe2B(int filletMethod,TopoDS_Face face1, TopoDS_Face face2,double exit_length, double transitionLength, double firstTransPos, int firstTransPosOption, double earlyTurnAngle)
{
	Standard_Real U1,U2,V1,V2;
	int i=0;
	
//	debugList.clear();

    Handle_Geom_Surface pSurf1 = BRep_Tool::Surface(face1);
	Handle_Geom_BSplineSurface pSurf1a = Handle(Geom_BSplineSurface)::DownCast(pSurf1->Copy());

	pSurf1a->Bounds(U1,U2,V1,V2);
	if(earlyTurnAngle < 1.0e-2)
		pSurf1a->Segment(U1+0.05*(U2-U1),U1+0.95*(U2-U1),V1,V2);


	Standard_Real U1a,U2a,V1a,V2a;
 /*   Handle_Geom_Surface pSurf2a = BRep_Tool::Surface(face2);
  //  Handle_Geom_BSplineSurface pSurf2 = Handle(Geom_BSplineSurface)::DownCast(pSurf2a->Copy());
	pSurf2a->Bounds(U1a,U2a,V1a,V2a);

	////////////////TEST
	
	gp_Pnt PntC1=pCrv1->Value(min(Uparam1,Uparam2));
	gp_Pnt PntC2=pCrv2->Value(max(Uparam1,Uparam2));
	double TestD1=Pta1.Distance(PntC1);
	double TestD2=Pta2.Distance(PntC2);   
	Handle_Geom_BSplineCurve Edge1Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf1a->VIso(V1)); 
	TopoDS_Edge StEdge1=BRepBuilderAPI_MakeEdge(Edge1Crv).Edge();
	
	Handle_Geom_BSplineCurve Edge2Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf2a->VIso(V1a)); 
	TopoDS_Edge StEdge2=BRepBuilderAPI_MakeEdge(Edge2Crv).Edge();
	
	BRepExtrema_ExtCC ExtCrv2(StEdge1,StEdge2);
		
		int NumofExt2=ExtCrv2.NbExt();
	//	if(NumofExt2==1)		
	//	{
		gp_Pnt Pnt1E2=ExtCrv2.PointOnE2(1);
	//	gp_Pnt Pnt2E2=ExtCrv2.PointOnE2(2);
		gp_Pnt Pnt1=ExtCrv2.PointOnE1(1);
		double E1Ux=ExtCrv2.ParameterOnE1(1);		*/

    
    Handle_Geom_Surface pSurf2 = BRep_Tool::Surface(face2);
  //  Handle_Geom_BSplineSurface pSurf2 = Handle(Geom_BSplineSurface)::DownCast(pSurf2a->Copy());
	pSurf2->Bounds(U1a,U2a,V1a,V2a);

    
 
	TopoDS_Shape sh1;
	try
	{
		//BRepAlgoAPI_Section makeBool1(pSurf1a/*face1*/,face2);
		BRepAlgoAPI_Section makeBool1(face1,face2);	//if we use surface result may odd here
		//makeBool1.Build(); ////this is not neccesary, automatically call from construction
		Standard_Boolean isDone = makeBool1.IsDone();	
		if(isDone)
			sh1 = makeBool1.Shape();
		else
			return -1;
	}
	catch(...)
	{
		return -1;
	}
	TopAbs_ShapeEnum shapeType = sh1.ShapeType();
	
 //	debugList.push_back(sh1);
 	ShapeFix_Shape fixShape(sh1);
 	fixShape.Perform();
 	sh1 = fixShape.Shape();

	//pick the intersection curve (select the longest edge if there is more than one
	TopoDS_Edge cutEdge;
	BRepBuilderAPI_MakeWire mkCutWire;
	i=0;
	for (TopExp_Explorer ex(sh1,TopAbs_EDGE); ex.More(); ex.Next())
	{
		TopoDS_Edge aEdge = TopoDS::Edge(ex.Current()); 
		mkCutWire.Add(aEdge);

		debugList.push_back(aEdge);
		
		if(i==0)
			cutEdge = aEdge;
		if(i > 0)
		{
			GProp_GProps LProps;
			BRepGProp::LinearProperties(aEdge,LProps);
			double length1 = LProps.Mass();

			BRepGProp::LinearProperties(cutEdge,LProps);
			double length2 = LProps.Mass();

			if(length1 > length2)
				cutEdge = aEdge;
		}
		i++;
	}
//	return -1;

	debugList.push_back(face1);
	debugList.push_back(face2);
	bool okEdge=false;
	if(i>0)
		{	
			double Ux1,Ux2;
			Handle_Geom_Curve pCutCrv = BRep_Tool::Curve(cutEdge,Ux1,Ux2);
			gp_Pnt Pta1 = pCutCrv->Value(Ux1);
			gp_Pnt Pta2 = pCutCrv->Value(Ux2);

			ShapeAnalysis_Curve SAC;
			gp_Pnt projPnt1,projPnt2;
			double Uparam1,Uparam2;
			
			double dist1_temp = SAC.Project(pSurf1->VIso(V1),Pta1,1.0e-5,projPnt1,Uparam1);
			double dist2_temp = SAC.Project(pSurf1->VIso(V1),Pta2,1.0e-5,projPnt2,Uparam2);
			if(dist1_temp <1.0e-6 && dist2_temp < 1.0e-6)
				okEdge=true;
		}
	if(mkCutWire.IsDone() && i > 1)
	{
		Handle_Geom_BSplineCurve tempCutCrv = makeSingleCurve(mkCutWire.Wire());
		if(!okEdge)
		cutEdge = BRepBuilderAPI_MakeEdge(tempCutCrv).Edge();
	}


//	debugList.push_back(cutEdge);
	Standard_Real Ut1,Ut2;
	BRep_Tool::Range(cutEdge,Ut1,Ut2);

	Handle_Geom_BSplineCurve ptempCrv = makeSingleCurve(BRepBuilderAPI_MakeWire(cutEdge).Wire());
	Ut1 = ptempCrv->FirstParameter();
	Ut2 = ptempCrv->LastParameter();


	cutEdge = BRepBuilderAPI_MakeEdge(ptempCrv).Edge();
	
	InterSecEdge=cutEdge;    // This is used for fillet operation
	
	BRep_Tool::Range(cutEdge,Ut1,Ut2);
//	debugList.push_back(cutEdge);

	//make the cut face on scroll
	Standard_Real Ua1,Ua2;
	Handle_Geom_Curve pCutCrv = BRep_Tool::Curve(cutEdge,Ua1,Ua2);
	gp_Pnt Pta1 = pCutCrv->Value(Ua1);
	gp_Pnt Pta2 = pCutCrv->Value(Ua2);

	ShapeAnalysis_Curve SAC;
	gp_Pnt projPnt1,projPnt2;
	Standard_Real Uparam1,Uparam2;
	double dist1_temp = SAC.Project(pSurf1->VIso(V1),Pta1,1.0e-5,projPnt1,Uparam1);
	double dist2_temp = SAC.Project(pSurf1->VIso(V1),Pta2,1.0e-5,projPnt2,Uparam2);
	if(dist1_temp > 1.0e-5 || dist2_temp > 1.0e-5)
		return false;

	

	Handle_Geom_BSplineCurve pCrv1 = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(V1)); 
	TopoDS_Edge edge1,edge2;
	if(fabs(U1-min(Uparam1,Uparam2)) > 0.01)
	{
		pCrv1->Segment(U1,min(Uparam1,Uparam2));
		edge1 = BRepBuilderAPI_MakeEdge(pCrv1).Edge();
	}

	Handle_Geom_BSplineCurve pCrv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(V1)); 
	if(fabs(U2-std::max(Uparam1,Uparam2)) > 0.01)
	{
		pCrv2->Segment(std::max(Uparam1,Uparam2),U2);
		edge2 = BRepBuilderAPI_MakeEdge(pCrv2).Edge();
	}
	
	BRepBuilderAPI_MakeWire mkWire(cutEdge);
	if(!edge1.IsNull())
		mkWire.Add(edge1);
	Standard_Boolean success = mkWire.IsDone();
	if(!edge2.IsNull())
		mkWire.Add(edge2);
	Standard_Boolean success2 = mkWire.IsDone();
	if(!success || !success2)
	{
		debugList.push_back(cutEdge);
		debugList.push_back(edge1);
		debugList.push_back(edge2);
		return -1;
	}

	//add the balance edges to close the wire on scroll surface
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->UIso(U1)).Edge());
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->UIso(U2)).Edge());
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSurf1->VIso(V2)).Edge());

	if(!mkWire.IsDone())
		return -1;

	//fix the wire shape
	TopoDS_Wire boundWire = mkWire.Wire();
	ShapeFix_Wire fixWire(boundWire,face1,1.0e-6);
	fixWire.Perform();
	boundWire = fixWire.Wire();

	BRepBuilderAPI_MakeFace mkFace1(pSurf1,boundWire/*mkWire.Wire()*/);
	success = mkFace1.IsDone();

	if(success)
	{
		TopoDS_Face tempface = mkFace1.Face();
		ShapeFix_Face fixFace1(tempface);
		fixFace1.Perform();
		Standard_Boolean status = fixFace1.Status(ShapeExtend_OK);
		trimFace1 = fixFace1.Face();
	}
	else
		return -1;

	///break in to segments (program calculated or user specified)
	double segFactor = 0.05;

	if(earlyTurnAngle > 1.0e-3)
	{
		//estimate a proper seg factor
		Handle_Geom_Curve pUIsoCrv = pSurf2->UIso(U1a+0.5*(U2a-U1a));

		GeomAPI_ExtremaCurveCurve GExt(pUIsoCrv,ptempCrv);
		int numPnts = GExt.NbExtrema();
		Standard_Real U1_x,U2_x;
		if(numPnts > 0)
		{
			GExt.Parameters(1,U1_x,U2_x);
			segFactor = 0.5*(U1_x-U1a)/(U2a-U1a);
		}
		else
		{
			segFactor = 0.05;
		}


	}
	else if(firstTransPosOption == 0) //program calculated
	{
		double pipe0_length = exit_length*transitionLength/100.0;

		gp_Pnt startSectionCentroid;

		GProp_GProps SProps;
        BRepBuilderAPI_MakeFace mkFace(exitTranSectWires[0]);
        BRepGProp::SurfaceProperties(mkFace.Face(),SProps);
        startSectionCentroid = SProps.CentreOfMass();

		//take the scroll cirumference
		double radius = sqrt(startSectionCentroid.X()*startSectionCentroid.X()+startSectionCentroid.Y()*startSectionCentroid.Y());
		double circum = 2.0*M_PI*radius;
		double firstSegLength = circum/36.0;
		segFactor = firstSegLength/pipe0_length;
		if(earlyTurnAngle > 1.0e-3)
			segFactor = segFactor*0.5;
		
		//segFactor = 0.0;
	}
	else //user specified
		segFactor = firstTransPos/100.0;




	if (filletMethod==0)  segFactor = 0.0;	// Modification to include Fillet operation using Thrusections
	

	Handle_Geom_Curve pVIsoCrv = pSurf2->VIso(V1a+segFactor*(V2a-V1a));

	//find the intersection of this curve with the cut Edge
	GeomAPI_ExtremaCurveCurve GECC(ptempCrv,pVIsoCrv);
	int numPoints = GECC.NbExtrema();
	Standard_Real U1_a,U1_b,U2_a,U2_b;
	if(numPoints < 2 )
	{
		//debuging
		debugList.push_back(BRepBuilderAPI_MakeEdge(ptempCrv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pVIsoCrv).Edge());

		//find intersection point for each segment
		Handle_Geom_BSplineCurve pCrv1 = Handle(Geom_BSplineCurve)::DownCast( pVIsoCrv->Copy()); 
		double u1,u2;
		u1 = pVIsoCrv->FirstParameter();
		u2 = pVIsoCrv->LastParameter();

		pCrv1->Segment(u1,(u1 + u2)/2.0);
		GeomAPI_ExtremaCurveCurve GECCFirst(ptempCrv,pCrv1);
		int numSegments = GECCFirst.NbExtrema();
		if(numSegments > 0)
		{
			GECCFirst.LowerDistanceParameters(U1_a,U2_a);
		}

		Handle_Geom_BSplineCurve pCrv2 = Handle(Geom_BSplineCurve)::DownCast( pVIsoCrv->Copy()); 
		pCrv2->Segment((u1 + u2)/2.0 , u2);
		GECCFirst.Init(ptempCrv,pCrv2);
		numSegments = GECCFirst.NbExtrema();
		if(numSegments > 0)
		{
			GECCFirst.LowerDistanceParameters(U1_b,U2_b);
		}
	}else
	{
		int index1,index2;
		double dist1 = 1.0E+5;
		double dist2 = 1.0E+5;

		for(int i=1; i<= numPoints; i++)
		{
			double dist = GECC.Distance(i);
			if(dist < dist1)
			{
				if(dist1 < dist2)
				{
					dist2 = dist1;
					index2 = index1;
				}
				dist1 = dist;
				index1 = i;
			}
			else if(dist < dist2)
			{
				dist2 = dist;
				index2 = i;
			}
		}

		GECC.Parameters(index1,U1_a,U2_a);
		GECC.Parameters(index2,U1_b,U2_b);
	}
/*	BRepBuilderAPI_MakeEdge cutEdge2(pVIsoCrv);
	BRepExtrema_ExtCC ExtCrv3(cutEdge,cutEdge2.Edge());
		
		int NumofExt3=ExtCrv3.NbExt();
	double Ux1=ExtCrv3.ParameterOnE2(1);
	double Ux2=ExtCrv3.ParameterOnE2(2);
	debugList.push_back(BRepBuilderAPI_MakeVertex(ExtCrv3.PointOnE2(1)).Vertex());
	debugList.push_back(BRepBuilderAPI_MakeVertex(ExtCrv3.PointOnE2(2)).Vertex());
	debugList.push_back(cutEdge2.Edge());
	debugList.push_back(cutEdge);
	*/
	//now make the line segments necessary
	Handle_Geom_BSplineCurve pSegCrv1 = Handle(Geom_BSplineCurve)::DownCast(pVIsoCrv->Copy());
	Handle_Geom_BSplineCurve pSegCrv2 = Handle(Geom_BSplineCurve)::DownCast(pVIsoCrv->Copy());
	
	pSegCrv1->Segment(U1a,min(U2_a,U2_b));
	pSegCrv2->Segment(std::max(U2_a,U2_b),U2a);

	Handle_Geom_BSplineCurve ptempCrvSeg = Handle(Geom_BSplineCurve)::DownCast(ptempCrv->Copy());
	Handle_Geom_BSplineCurve ptempCrvSeg2 = Handle(Geom_BSplineCurve)::DownCast(ptempCrv->Copy());
	Handle_Geom_BSplineCurve ptempCrvSeg3 = Handle(Geom_BSplineCurve)::DownCast(ptempCrv->Copy());
	

	ptempCrvSeg->Segment(min(U1_a,U1_b),std::max(U1_a,U1_b));
	ptempCrvSeg2->Segment(Ut1,min(U1_a,U1_b));
	ptempCrvSeg3->Segment(std::max(U1_a,U1_b),Ut2);


	BRepBuilderAPI_MakeWire mkSegWire;
	BRepBuilderAPI_MakeWire mkSegWire2;
	BRepBuilderAPI_MakeWire mkSegWire3;

	
	//add the balance edges to close the wire on exit pipe surface
	Handle_Geom_BSplineCurve pSide1Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U1a)->Copy());
	Handle_Geom_BSplineCurve pSide2Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U2a)->Copy());

	Handle_Geom_BSplineCurve pSide1Crv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U1a)->Copy());
	Handle_Geom_BSplineCurve pSide2Crv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U2a)->Copy());

	pSide1Crv->Segment(V1a+segFactor*(V2a-V1a),V2a);
	pSide2Crv->Segment(V1a+segFactor*(V2a-V1a),V2a);

	pSide1Crv2->Segment(V1a,V1a+segFactor*(V2a-V1a));
	pSide2Crv2->Segment(V1a,V1a+segFactor*(V2a-V1a));

	if(earlyTurnAngle > 1.0e-3)
	{
		Handle_Geom_Curve pCrvC = pSurf1->VIso(V1);
		Pta1 = pCrvC->Value(pCrvC->FirstParameter());
		Pta2 = pCrvC->Value(pCrvC->LastParameter());
	}


	dist1_temp = SAC.Project(pSurf2->VIso(V1a),Pta1,1.0e-5,projPnt1,Uparam1);
	dist2_temp = SAC.Project(pSurf2->VIso(V1a),Pta2,1.0e-5,projPnt2,Uparam2);
	if(dist1_temp > 1.0e-5 || dist2_temp > 1.0e-5)
		return false;

	Handle_Geom_BSplineCurve pCrv3 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv3->Segment(U1a,min(Uparam1,Uparam2));
	TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(pCrv3).Edge();

	Handle_Geom_BSplineCurve pCrv4 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv4->Segment(std::max(Uparam1,Uparam2),U2a);
	TopoDS_Edge edge4 = BRepBuilderAPI_MakeEdge(pCrv4).Edge();


	//main wire
	bool WireOK;
	
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(ptempCrvSeg /*ptempCrv*/).Edge());
	WireOK=mkSegWire.IsDone();
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSegCrv1).Edge());
//	mkSegWire.Add(BRepBuilderAPI_MakeEdge(ptempCrv->Value(ptempCrv->FirstParameter()),pSide1Crv->Value(pSide1Crv->FirstParameter())).Edge());
	WireOK=mkSegWire.IsDone();
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSegCrv2).Edge());
	WireOK=mkSegWire.IsDone();
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSide1Crv).Edge());
	WireOK=mkSegWire.IsDone();
	
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSide2Crv).Edge());
	WireOK=mkSegWire.IsDone();
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSurf2->VIso(V2a)).Edge());
	WireOK=mkSegWire.IsDone();

	
	
	if(!mkSegWire.IsDone())
	{
		debugList.push_back(BRepBuilderAPI_MakeEdge(ptempCrvSeg).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSegCrv1).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSegCrv2).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSide1Crv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSide2Crv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSurf2->VIso(V2a)).Edge());
		return -1;
	}


	//first sub wire
	if (!filletMethod==0)			// Modification to include Fillet operation using Thrusections
	{
		mkSegWire2.Add(BRepBuilderAPI_MakeEdge(pSegCrv1).Edge());
		mkSegWire2.Add(BRepBuilderAPI_MakeEdge(pSide1Crv2).Edge());
		mkSegWire2.Add(BRepBuilderAPI_MakeEdge(ptempCrvSeg2).Edge());
		if(!mkSegWire2.IsDone())
			mkSegWire2.Add(BRepBuilderAPI_MakeEdge(ptempCrvSeg3).Edge());
		if(earlyTurnAngle > 1.0e-3 && !edge1.IsNull())
			mkSegWire2.Add(edge1);
		if(!mkSegWire2.IsDone())
			return -1;

		mkSegWire2.Add(edge3);


		//second sub wire
		mkSegWire3.Add(BRepBuilderAPI_MakeEdge(pSegCrv2).Edge());
		mkSegWire3.Add(BRepBuilderAPI_MakeEdge(pSide2Crv2).Edge());
		mkSegWire3.Add(BRepBuilderAPI_MakeEdge(ptempCrvSeg3).Edge());
		if(!mkSegWire3.IsDone())
			mkSegWire3.Add(BRepBuilderAPI_MakeEdge(ptempCrvSeg2).Edge());
		if(earlyTurnAngle > 1.0e-3 && !edge2.IsNull())
			mkSegWire3.Add(edge2);
		if(!mkSegWire3.IsDone())
			return -1;

		mkSegWire3.Add(edge4);
	}


	//fix the wire shape
	boundWire = mkSegWire.Wire();	
	
	ShapeFix_Wire fixWire2(boundWire,face2,1.0e-6);
	fixWire2.Perform();
	boundWire = fixWire2.Wire();
	
	
	BRepBuilderAPI_MakeFace mkFace2(pSurf2,boundWire/*mkSegWire.Wire()*/);	
	if(mkFace2.IsDone())
	{
		TopoDS_Face tempface = mkFace2.Face();
		ShapeFix_Face fixFace2(tempface);
		fixFace2.Perform();
		Standard_Boolean status = fixFace2.Status(ShapeExtend_OK);
		trimFace2 = fixFace2.Face();
	}
	else
		return -1;
	
	if (!filletMethod==0)		// Modification to include Fillet operation using Thrusections
	{	
		BRepBuilderAPI_MakeFace mkFace3(pSurf2,mkSegWire2.Wire());	
		if(mkFace3.IsDone())
		{
			TopoDS_Face tempface = mkFace3.Face();
			ShapeFix_Face fixFace3(tempface);
			fixFace3.Perform();
			Standard_Boolean status = fixFace3.Status(ShapeExtend_OK);
			tongueRegionPipeFaces.push_back(fixFace3.Face());
		}
		else
			return -1;

	 
		BRepBuilderAPI_MakeFace mkFace4(pSurf2,mkSegWire3.Wire());	
		if(mkFace4.IsDone())
		{
			TopoDS_Face tempface = mkFace4.Face();
			ShapeFix_Face fixFace4(tempface);
			fixFace4.Perform();
			Standard_Boolean status = fixFace4.Status(ShapeExtend_OK);
			tongueRegionPipeFaces.push_back(fixFace4.Face());
		}
		else
			return -1;
	}
	
	
	
	return 0;


}


bool OCCVolute::sewAllFaces(bool isClockWise, bool doFillets,int filletMethod,double tongueCurvControlParam,double tongueSizeControlParam, int filletOption, double filletRmin, double filletRmax, bool exit_enable_extension)
{
	//do the sewing
	Standard_Boolean isFilleted = Standard_False;
	TopoDS_Shape filletedShape;
	filletedShapeSaved.Nullify();

	if(!doFillets)
	{
		BRepBuilderAPI_Sewing sewing1(0.005/*1.0e-6*/);

		for(int i=0; i<scrollFaces.size(); i++)
   	 		sewing1.Add(scrollFaces[i]);

		sewing1.Add(trimFace1);
		sewing1.Add(trimFace2);

		for(int i=0; i<tongueRegionScrollFaces.size(); i++)
			sewing1.Add(tongueRegionScrollFaces[i]);

		for(int i=0; i<tongueRegionPipeFaces.size(); i++)
			sewing1.Add(tongueRegionPipeFaces[i]);

		for(int i=0; i<exitPipeFaces.size(); i++)
			sewing1.Add(exitPipeFaces[i]);

		sewing1.Perform();


		TopoDS_Shape sh1 = sewing1.SewedShape();
		ShapeFix_Shape shf1(sh1);
		shf1.Perform();
		sh1 = shf1.Shape();

		totalVolute = sh1;   

		m_exitPlane = exitPipeFaces[exitPipeFaces.size() - 1];	//used for stl file export
	}
	
	else  // add fillets to tongue region
	{
		Standard_Real Radmin,Radmax;
		if(filletOption == 0)
		{
			Radmin = 0.001;
			Radmax = 0.002;
		}
		else
		{
			Radmin = filletRmin;
			Radmax = filletRmax;
		}

		BRepBuilderAPI_Sewing sewing(0.005/*1.0e-6*/);


		//do the fillet operation
		
		int edgeNum = 0;
		
		if(filletMethod==0)	
		{
		
			isFilleted=makeTongueFilletUsingThrughSections(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2),tongueCurvControlParam,tongueSizeControlParam, filletedShape);
			
			
		} 
		
				
		else if(filletMethod==1)
		{
			try
			{
			/*	isFilleted = makeTongueFillet(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2), 1,Radmin, Radmax, Radmin, filletedShape, edgeNum);
				if(filletOption == 0 && !isFilleted) //if the default fillet values did not work try another set of values
				{
					Radmin = 0.0001;
					Radmax = 0.001;
					isFilleted = makeTongueFillet(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2), 1,Radmin, Radmax, Radmin, filletedShape, edgeNum);
				}
				if(filletOption == 0 && !isFilleted) //if the default fillet values did not work try another set of values
				{
					Radmin = 0.0001;
					Radmax = 0.0005;
					isFilleted = makeTongueFillet(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2), 1,Radmin, Radmax, Radmin, filletedShape, edgeNum);
				}
				if(filletOption == 0 && !isFilleted) //if the default fillet values did not work try another set of values
				{
					Radmin = 0.0001;
					Radmax = 0.0001;
					isFilleted = makeTongueFillet(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2), 1,Radmin, Radmax, Radmin, filletedShape, edgeNum);
				} */
					
				if(filletOption == 0 && !isFilleted)
				{
					Radmin = 0.0075;
					Radmax = 0.0075;
					for(int i=0; i<10; i++)
					{
						isFilleted = makeTongueFillet(TopoDS::Face(trimFace1), TopoDS::Face(trimFace2), 1,Radmin, Radmax, Radmin, filletedShape, edgeNum);
						
						if(isFilleted)
							break;
						else
						{
							Radmin = Radmin*0.75;
							Radmax = Radmax*0.75;
						}
					}
				}

						
			}
			catch(...)
			{
				isFilleted = Standard_False;
			}
		}

		if(isFilleted)
		{
			sewing.Add(filletedShape);
			filletedShapeSaved = filletedShape;		
		}
		else
		{
 			sewing.Add(trimFace1);		
 			sewing.Add(trimFace2);
		}

 
  		//add the balance faces	

		for(int i=0; i<tongueRegionScrollFaces.size(); i++)
 			sewing.Add(tongueRegionScrollFaces[i]);  
  		for(int i=0; i<tongueRegionPipeFaces.size(); i++)
 			sewing.Add(tongueRegionPipeFaces[i]);
	
   		for(int i=0; i<scrollFaces.size(); i++)
			sewing.Add(scrollFaces[i]);

		m_exitPlane = exitPipeFaces[exitPipeFaces.size() - 1];

		for(int i=0; i<exitPipeFaces.size() ; i++)
			sewing.Add(exitPipeFaces[i]);

 
		sewing.Perform();
		TopoDS_Shape sh1 = sewing.SewedShape();
 		ShapeFix_Shape shf1(sh1);
		shf1.Perform();
		sh1 = shf1.Shape();
 
		totalVolute = sh1;


	}// end fillet process
		
//	debugList.push_back(totalVolute);
		

	
	//do a cleanup of the sewed shape and make a solid
  	if (totalVolute.ShapeType() == TopAbs_SHELL)
	{ 
		TopoDS_Shell shell = TopoDS::Shell(totalVolute); 

 		//start test
		try
		{
			ShapeUpgrade_ShellSewing SUSS;
			Standard_Real tol = 0.0;
			TopoDS_Shape shape = SUSS.ApplySewing(shell,tol);
			Standard_Boolean clsd = shape.Closed();

	 		shell = TopoDS::Shell(shape);
		}
		catch(...)
		{
		}


		ShapeAnalysis_Shell sas;
		sas.LoadShells(shell);
		//analysis of the shell , second parameter is set to True for //getting free edges,(default False)
		sas.CheckOrientedShells(shell,Standard_True);
		//getting the result of analysis
		Standard_Boolean hasbadEdges = sas.HasBadEdges();
		Standard_Boolean hasfreeEdges = sas.HasFreeEdges();

		

/*		ShapeFix_Wireframe fixWireFrame(shell);
 		fixWireFrame.SetPrecision ( 1.0e-8 );
 		fixWireFrame.SetMaxTolerance ( 1.0e-3 );
		fixWireFrame.ModeDropSmallEdges() = Standard_True;
		Standard_Boolean isdone = fixWireFrame.FixSmallEdges();
		isdone = fixWireFrame.FixWireGaps();
		
		shell = TopoDS::Shell(fixWireFrame.Shape());

		Standard_Boolean isClosed = shell.Closed();		

		ShapeFix_Shell fixShell(shell);
 		fixShell.SetPrecision ( 1.0e-8 );
 		fixShell.SetMaxTolerance ( 1.0e-3 );

		Standard_Boolean isDone = fixShell.Perform();
*/		
		BRepBuilderAPI_MakeSolid mkSolid(shell/*fixShell.Shell()*/);
 		
 		totalVolute = mkSolid.Solid();
		Standard_Boolean isClosed = totalVolute.Closed();		
		
		//end test







/*
		TopoDS_Shell shell = TopoDS::Shell(totalVolute); 
		ShapeFix_Shell fixShell(shell);
		Standard_Boolean isDone = fixShell.Perform();
		BRepBuilderAPI_MakeSolid mkSolid(fixShell.Shell());
       // TopoDS_Solid solid = ShapeFix_Solid().SolidFromShell(shell);
		totalVolute = mkSolid.Solid();//solid;
		*/
	}


//    writeStepFile();	
	
	return true;
}

bool OCCVolute::makeTongueFillet(TopoDS_Face face1, TopoDS_Face face2, int filletProfile, double R1, double R2, double R3, TopoDS_Shape& filleted, int edgeNum)
{	
	
	
	BRepBuilderAPI_Sewing sewing(1.0e-4);
	
	sewing.Add(face1);
	sewing.Add(face2);
	sewing.Perform();
	

	TopoDS_Shape sh1 = sewing.SewedShape();
 	if (sh1.ShapeType() == TopAbs_SHELL)
	{ 
		TopoDS_Shell shell = TopoDS::Shell(sh1); 
		ShapeFix_Shell fixShell(shell);
		Standard_Boolean isDone = fixShell.Perform();
		sh1 = fixShell.Shape();
	}
	
	
//	debugList.push_back(sh1);
	BRepFilletAPI_MakeFillet mkFillet(sh1);

	Standard_Boolean success = 0;
	TopoDS_Edge edgeToFillet;
	int j=0;

	
	//Now figure out which edge to fillet. Ideally this should be straightforward with sewing.NbContiguousEdges() funciton. But it doesn't work properly for now.
	// So take the long route of DIY.
	

	TopoDS_Face subFace1, subFace2;
	TopExp_Explorer Facexplorer(sh1 , TopAbs_FACE);
	j=0;
	while( Facexplorer.More() )
	{
		TopoDS_Face aFace = TopoDS::Face(Facexplorer.Current());
		if(j==0 )
			subFace1 = aFace;
		else
			subFace2 = aFace;
		j++;
		Facexplorer.Next();
	}


	bool foundEdge = false;
	TopExp_Explorer EdgeExplorer(subFace1 , TopAbs_EDGE);
	while( EdgeExplorer.More() )
	{
		TopoDS_Edge aEdge = TopoDS::Edge(EdgeExplorer.Current());

		Standard_Real Ut1,Ut2;
		BRep_Tool::Range(aEdge,Ut1,Ut2);

		TopExp_Explorer EdgeExplorer2(subFace2 , TopAbs_EDGE);
		while( EdgeExplorer2.More() )
		{
			TopoDS_Edge aEdge2 = TopoDS::Edge(EdgeExplorer2.Current());
			if(aEdge.IsSame(aEdge2))
			{
				edgeToFillet = aEdge;
				foundEdge = true;
				break;
			}
			EdgeExplorer2.Next();		
		}

		EdgeExplorer.Next();
		if(foundEdge)
			break;
	}

	if(!foundEdge)
		return false;


	if(filletProfile == 0) //linear profile
	{
		mkFillet.Add(R1,R1,edgeToFillet);
	}
	else if(filletProfile == 1) //3 radius values
	{
		Standard_Real U1edge,U2edge;
		BRep_Tool::Range(edgeToFillet,U1edge,U2edge);

		//the fillet profile will start and end with 0.0 R values for pump volute
	 	R1 = 0.0;
	 	R3 = 0.0;

		TColgp_Array1OfPnt2d UandR(1,3);
		UandR.SetValue(1,gp_Pnt2d(U1edge,R1));
		UandR.SetValue(2,gp_Pnt2d(0.5*(U1edge+U2edge),R2));
		UandR.SetValue(3,gp_Pnt2d(U2edge,R3));

		mkFillet.Add(UandR,edgeToFillet);
	}
	else if(filletProfile == 2) //5 radius values
	{
		Standard_Real U1edge,U2edge;
		BRep_Tool::Range(edgeToFillet,U1edge,U2edge);

		//the fillet profile will start and end with 0.0 R values for pump volute
	 	R1 = 0.0001;
	 	R3 = 0.0001;

		TColgp_Array1OfPnt2d UandR(1,7);
		UandR.SetValue(1,gp_Pnt2d(U1edge,0.0));
		UandR.SetValue(2,gp_Pnt2d(U1edge+0.1*(U2edge-U1edge),0.0005));
		UandR.SetValue(3,gp_Pnt2d(U1edge+0.25*(U2edge-U1edge),0.001));
		UandR.SetValue(4,gp_Pnt2d(0.5*(U1edge+U2edge),R2));
		UandR.SetValue(5,gp_Pnt2d(U1edge+0.75*(U2edge-U1edge),0.001));
		UandR.SetValue(6,gp_Pnt2d(U1edge+0.9*(U2edge-U1edge),0.0005));
		UandR.SetValue(7,gp_Pnt2d(U2edge,0.0));

		mkFillet.Add(UandR,edgeToFillet);
	}

	try
	{
		mkFillet.Build();
	}
	catch(...)
	{
		return false;
	}
	success = mkFillet.IsDone();

 	if(success)
	{
		Standard_Boolean isdone;
		filleted = mkFillet.Shape();
	
		//this is a test to make sure the fillet was okay. this will not catch all problems but some
		int ctr = 0;
		for(TopExp_Explorer Expl(filleted,TopAbs_FACE);Expl.More(); Expl.Next())
		{
			TopoDS_Face tempFace = TopoDS::Face(Expl.Current());
			ctr++;
		}

		if(ctr < 3)	//something wrong. there should be at least 3 faces
			return 0;

		try
		{
			ShapeFix_Shape fixShape(filleted);
			isdone = fixShape.Perform();
			filleted = fixShape.Shape();
	

 			ShapeFix_Wireframe fixWireFrame(filleted);
			fixWireFrame.ModeDropSmallEdges() = Standard_True;
			isdone = fixWireFrame.FixSmallEdges();
			isdone = fixWireFrame.FixWireGaps();
			
			filleted = fixWireFrame.Shape();

		}
		catch(...)
		{
			return 0;
		}
	}  
 
	
	return success;
}



int OCCVolute::makeOuterShell(double wallThickness, bool isTurbine, bool exit_straight, bool isClockWise)
{
	int errorCode = 0;
//	debugList.clear();

	//offset testing
	Standard_Real Offset1 = wallThickness; 
	Standard_Real Offset2 = -1.0*wallThickness;
	if((isTurbine && !isClockWise) || (!isTurbine && isClockWise))
		Offset2 = wallThickness;
	Standard_Real Tol = 1.0e-6;
	BRepOffset_Mode Mode = BRepOffset_Skin;
	Standard_Boolean Intersection = Standard_False;
	Standard_Boolean SelfInter = Standard_False;
	GeomAbs_JoinType Join = GeomAbs_Arc;

	///make a sewed shell
	Standard_Real sewTol = 5.0e-3;
	BRepBuilderAPI_Sewing sewing3(/*0.005*/sewTol);
	BRepBuilderAPI_Sewing sewing4(/*0.005*/1.0e-3);
	BRepBuilderAPI_Sewing innerShellSewing( 0.005 /*1.0e-4*/);

	TopoDS_Shape resultShape;
	TopoDS_Face scrollFace1Exp;
	TopoDS_Shape exitPipeExp;
	TopoDS_Shape innerScrollShell;


	//remake the 1st scroll segment expanded
	int i_mid = 0;
	double angleBreakDegrees = 60.0;
	double angleBreak = angleBreakDegrees*M_PI/180.0; //M_PI/3.0;
	for(int i=0; i<scrollSectWires.size();i++)
	{
		if(thetaSects[i] > angleBreak)
		{
			i_mid = i-1;
			break;
		}
    }

 	BRepOffsetAPI_ThruSections thruSectsScrollSeg1(Standard_False,Standard_True,1.0e-8);

	for(int i=0; i< scrollSectWires.size(); i++)
		thruSectsScrollSeg1.AddWire(scrollSectWires[i]);

	thruSectsScrollSeg1.Build();
 	
	Standard_Boolean isRuled = Standard_False;

	BRepOffsetAPI_ThruSections BTS1(Standard_False,isRuled,1.0e-6);
	BRepOffsetAPI_ThruSections BTS2(Standard_False,isRuled,1.0e-6);
	int j = 0;
	int numFaces = 0;
	Handle_Geom_Curve pScrollSectLastExp;
	Handle_Geom_Curve pScrollSectGoodExp;
	for(TopExp_Explorer explr(thruSectsScrollSeg1.Shape(),TopAbs_FACE); explr.More(); explr.Next())
	{
		numFaces++;
	}

	//volute start radius
	double baseRadius = 0.0;
	for(TopExp_Explorer explr(seg0Wires[seg0Wires.size()-1],TopAbs_EDGE); explr.More(); explr.Next())
	{
		TopoDS_Edge aEdge = TopoDS::Edge(explr.Current());
		Standard_Real Ut1,Ut2;
		Handle_Geom_Curve ptCrv = BRep_Tool::Curve(aEdge,Ut1,Ut2);
		gp_Pnt tempPt = ptCrv->Value(ptCrv->FirstParameter());
		baseRadius = sqrt(tempPt.X()*tempPt.X() + tempPt.Y()*tempPt.Y());
	}


	bool startExpSecOK = false; //indicates if start expanded section is ready
	 for(TopExp_Explorer explr(thruSectsScrollSeg1.Shape(),TopAbs_FACE); explr.More(); explr.Next())
	{
		TopoDS_Face tFace = TopoDS::Face(explr.Current());

		if(j%2 == 0)
		{
			TopoDS_Face tFaceExp;
			try
			{
				BRepOffsetAPI_MakeOffsetShape mkOffset3(tFace,Offset2,Tol,Mode,Intersection,SelfInter,Join);
				for(TopExp_Explorer explr(mkOffset3.Shape(),TopAbs_FACE);explr.More();explr.Next())
				{
					tFaceExp = TopoDS::Face(explr.Current());
				}
			}
			catch(...)
			{
				debugList.push_back(tFace);
				if(j == 0)
				{
				}
				else if(j == 2*i_mid)
				{
					i_mid++;
					tFaceExp.Nullify();
				}
				else if((j == numFaces-2) || (j == numFaces-1)) //replace the failed curve with the last good one
				{
					Handle_Geom_Surface pFailedSurf = BRep_Tool::Surface(tFace);
					Standard_Real U1f,U2f,V1f,V2f;
					pFailedSurf->Bounds(U1f,U2f,V1f,V2f);
					Handle_Geom_Curve pCrvB = pFailedSurf->VIso(V2f);

					gp_Pnt Pt1 = pCrvB->Value(U1f);
					gp_Vec2d VecF1 = gp_Vec2d(gp_Pnt2d(0,0),gp_Pnt2d(Pt1.X(),Pt1.Y()));

					gp_Pnt Pt2 = pScrollSectGoodExp->Value(pScrollSectGoodExp->FirstParameter());
					gp_Vec2d VecG1 = gp_Vec2d(gp_Pnt2d(0,0),gp_Pnt2d(Pt2.X(),Pt2.Y()));

					Standard_Real angRot = VecG1.Angle(VecF1);
					gp_Ax1 Ax1(gp_Pnt(0,0,0),gp_Vec(0,0,1));
					pScrollSectLastExp = Handle(Geom_Curve)::DownCast(pScrollSectGoodExp->Copy());
					pScrollSectLastExp->Rotate(Ax1,angRot);

					BTS2.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pScrollSectLastExp).Edge()).Wire());
				}
				else
					tFaceExp.Nullify();

			}
			
			

	//		if(j >= 2*i_mid)
	//			innerShellSewing.Add(tFace);

			Standard_Real U1b,U2b,V1b,V2b;
			Handle_Geom_Curve pCrvA;
			Handle_Geom_Surface presultSurf;
			if(!tFaceExp.IsNull())
			{
				presultSurf = BRep_Tool::Surface(tFaceExp);
				presultSurf->Bounds(U1b,U2b,V1b,V2b);
			}
			

			if(presultSurf.IsNull())
			{
			}
			else if(/*j == 0*/ !startExpSecOK && (j == 4 || j == 6 || j == 8)) //j=0 section fails on many expansion operations
			{
				//get any expanded section close to start section at rotate it to align with the start angle
				startExpSecOK = true;
				pCrvA = presultSurf->VIso(0.5*(V1b+V2b));
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				
				//rotate this curve to the scroll start position
				gp_Pnt Pt1(0,0,pCrvA->Value(pCrvA->FirstParameter()).Z());
				gp_Pnt Pt2(pCrvA->Value(pCrvA->FirstParameter()).X(),0,pCrvA->Value(pCrvA->FirstParameter()).Z());
				gp_Vec v1(Pt1,pCrvA->Value(pCrvA->FirstParameter()));
				gp_Vec v2(Pt1,Pt2);
				Standard_Real rotAng = v1.Angle(v2);
				if((isTurbine && isClockWise) || (!isTurbine && isClockWise))
					rotAng = -1.0*rotAng;
				pCrvA->Rotate(gp::OZ(),rotAng);
				BTS1.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());

			 /*	pCrvA = presultSurf->VIso(V1b);
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				BTS1.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());  */
			}
			else if((j == i_mid) || (j == i_mid+1))
			{
				pCrvA = presultSurf->VIso(0.5*(V1b+V2b));
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				BTS1.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());
			}
			else if(j == 2*i_mid)
			{
				pCrvA = presultSurf->VIso(0.5*(V1b+V2b));
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				BTS1.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());
				BTS2.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());
			}
			else if((j == numFaces-2) || (j == numFaces-1))
			{
				pCrvA = presultSurf->VIso(V2b);
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				BTS2.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());
				pScrollSectLastExp = pCrvA; // save for later use
			}
			else if(j > 2*i_mid)
			{
				pCrvA = presultSurf->VIso(0.5*(V1b+V2b));
				pCrvA = makeProjCurve(pCrvA,baseRadius);
				pScrollSectGoodExp = pCrvA;
				if(j%4 == 0)
					BTS2.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(pCrvA).Edge()).Wire());
			}
		}
		j++;
	}

	try
	{
		BTS1.Build();
		BTS2.Build();
	}
	catch(...)
	{
		return 4;
	}
	

	for(TopExp_Explorer explr(BTS1.Shape(),TopAbs_FACE);explr.More();explr.Next())
	{
		scrollFace1Exp = TopoDS::Face(explr.Current());
	}

	TopoDS_Shape innerShellShape;	
	if(!hasSplitter)
		innerShellShape = makeInnerShellNoSplitter(pScrollSectLastExp);
	else innerShellShape = makeInnerShellWithSplitter(pScrollSectLastExp);
	

	//offset exit pipe
	for(int i=0; i<exitPipeFaces.size()-1; i++){
 		sewing4.Add(exitPipeFaces[i]);
	}

	sewing4.Perform();

	BRepOffsetAPI_MakeOffsetShape mkOffset4(sewing4.SewedShape(),Offset1,Tol,Mode,Intersection,SelfInter,Join);
	exitPipeExp = mkOffset4.Shape();


	//make the new transition segment(expanded)
	TopoDS_Shape transPipeShapeExp = makeTransitionExp(Offset1,exit_straight,scrollFace1Exp,exitPipeExp, pScrollSectLastExp);
	if(transPipeShapeExp.IsNull())
		return 2;


	//make the exit flange
	TopoDS_Shape flangedExitPipe = makeFlangePipe(exitPipeExp,Offset1);
	if(flangedExitPipe.IsNull())
		errorCode = 3;

	//do the boolean operation
	i = 0;
	TopoDS_Face trans1Face;
	for(TopExp_Explorer explr(transPipeShapeExp,TopAbs_FACE);explr.More();explr.Next())
	{
		if(i == 0)
			trans1Face = TopoDS::Face(explr.Current());
		else
			sewing3.Add(TopoDS::Face(explr.Current()));

		i++;
	}
	
 	TopoDS_Shape filletedShapeExp = makeBooleanExp(scrollFace1Exp, trans1Face, transStartWireExp);
	if(filletedShapeExp.IsNull())
		return 1;

 	sewing3.Add(filletedShapeExp);

	//add the rest
 	sewing3.Add(BTS2.Shape() );
	sewing3.Add(flangedExitPipe);
 	sewing3.Add(innerShellShape);


	//get the sewed shape
	sewing3.Perform();
	totalVoluteThick = sewing3.SewedShape();
	

	//make a solid
 	if(totalVoluteThick.ShapeType() == TopAbs_SHELL)
	{
		ShapeFix_Wireframe fixWireFrame(totalVoluteThick);
 		fixWireFrame.SetPrecision ( 1.0e-8 );
 		fixWireFrame.SetMaxTolerance ( 1.0e-3 );
		fixWireFrame.ModeDropSmallEdges() = Standard_True;
		Standard_Boolean isdone = fixWireFrame.FixSmallEdges();
		isdone = fixWireFrame.FixWireGaps();
		
		TopoDS_Shell shell = TopoDS::Shell(fixWireFrame.Shape());		

		Standard_Boolean isClosed = shell.Closed();		

		ShapeFix_Shell fixShell(shell);
 		fixShell.SetPrecision ( 1.0e-8 );
 		fixShell.SetMaxTolerance ( 1.0e-3 );

		Standard_Boolean isDone = fixShell.Perform();
		
		BRepBuilderAPI_MakeSolid mkSolid(fixShell.Shell());
 		
 		totalVoluteThick = mkSolid.Solid();
	}
 
//	debugList.push_back(totalVoluteThick);
	
//	writeStepFile();

	return errorCode;
}

TopoDS_Shape OCCVolute::makeTransitionExp(double offSet1,bool exit_straight,TopoDS_Face scrollFace1Exp, TopoDS_Shape exitPipeExp, Handle_Geom_Curve pScrollSectLastExp)
{
	TopoDS_Shape transPipeShapeExp;

	//get the starting X section
	TopAbs_ShapeEnum shapeType = scrollFace1Exp.ShapeType();
	Handle_Geom_Surface pGeomSurf = BRep_Tool::Surface(scrollFace1Exp);
	Standard_Real U1,U2,V1,V2;
	pGeomSurf->Bounds(U1,U2,V1,V2);
	Handle_Geom_Curve pIsoCrv1 = pGeomSurf->VIso(V1);

	Handle_Geom_Curve pIsoCrv2 = pScrollSectLastExp;

	transStartWireExp = makeExitStartSectionExp(pIsoCrv1,pIsoCrv2);
	TopoDS_Wire transEndWireExp = makeExitEndSectionExp(exitPipeExp);
	 
	Handle_Geom_BSplineCurve ptrEndBsp = makeSingleCurve(transEndWireExp);
 	transEndWireExp = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(ptrEndBsp).Edge()).Wire();
	 


	if(transStartWireExp.IsNull() || transEndWireExp.IsNull())
		return transPipeShapeExp;


	

	//do some tests
	TopoDS_Wire wire2 = TopoDS::Wire(exitTranSectWires[1]);
	BRepBuilderAPI_MakeFace mkFace1(wire2);
	mkFace1.Build();
	TopoDS_Face XsecFace;
	if(mkFace1.IsDone())
		XsecFace = mkFace1.Face();

	GProp_GProps SProps;
    BRepGProp::SurfaceProperties(XsecFace,SProps);
    double area1 = SProps.Mass();
	gp_Pnt Pt1 = SProps.CentreOfMass();

	BRepBuilderAPI_MakeFace mkFace2(transEndWireExp);
	mkFace2.Build();
	TopoDS_Face XsecFaceOld;
	double area2;
	double scaleRatio;
	if(mkFace2.IsDone())
	{
		XsecFaceOld = mkFace2.Face();
		BRepGProp::SurfaceProperties(XsecFaceOld,SProps);
		area2 = SProps.Mass();
		scaleRatio = area2/area1;
	}
	else
	{
		double R1 = sqrt(area1/M_PI);
		scaleRatio = (R1+offSet1)*(R1+offSet1)/(R1*R1);
	}
	
	scaleRatio = sqrt(scaleRatio); //(1.0+scaleRatio)/2.0;


	gp_Trsf  t1;
	if(mkFace2.IsDone())
		t1.SetScale(Pt1,scaleRatio);
	else
		t1.SetScale(Pt1,1.05);
    BRepBuilderAPI_Transform translate(wire2,t1); 
	wire2 = TopoDS::Wire(translate.Shape());


	if(exit_straight)
	{
		BRepOffsetAPI_ThruSections thruSections1(Standard_False,Standard_True,1.0e-6);
		thruSections1.AddWire(transStartWireExp );
		thruSections1.AddWire(wire2);

		
		try
		{
			thruSections1.Build();
			transPipeShapeExp = thruSections1.Shape();
		}
		catch(...)
		{
			debugList.push_back(transStartWireExp);
			debugList.push_back(transEndWireExp);
			return TopoDS_Shape();
		}
	}
	else
	{
		BRepOffsetAPI_MakePipeShell transPipeShellExp(exitTransSpineWires[0]);
        transPipeShellExp.Add(transStartWireExp,Standard_False,Standard_False) ;
        transPipeShellExp.Add(transEndWireExp/*wire2*/,Standard_False,Standard_True) ;
        
		try
		{
			transPipeShellExp.Build();
			transPipeShapeExp = transPipeShellExp.Shape();
		}
		catch(...)
		{
			debugList.push_back(transStartWireExp);
			debugList.push_back(transEndWireExp);
			debugList.push_back(exitTransSpineWires[0]);
			return TopoDS_Shape();
		}
	}


	return transPipeShapeExp;

}



TopoDS_Wire OCCVolute::makeExitStartSectionExp(Handle_Geom_Curve pscrollStartCrv, Handle_Geom_Curve pscrollEndCrv)
{
	Handle_Geom_BSplineCurve pSec1Curve = Handle_Geom_BSplineCurve::DownCast(pscrollStartCrv->Copy());
	Handle_Geom_BSplineCurve pSecLastCurve = Handle_Geom_BSplineCurve::DownCast(pscrollEndCrv->Copy());
	
	int i=0;

/*	Handle_TColgp_HArray1OfPnt pSec1Points = new TColgp_HArray1OfPnt(1,scrollStartSect.npt);

	int i=0;
	int index = 1;
	for(i=0; i< scrollStartSect.npt; i++)
	{
		gp_Pnt tempPt(scrollStartSect.x[i],scrollStartSect.y[i],scrollStartSect.z[i]);
		pSec1Points->SetValue(index,tempPt);
		index++;
	}
			
					
	GeomAPI_Interpolate interpSec1(pSec1Points,Standard_False,1.0e-6);
	interpSec1.Perform();
	Handle_Geom_BSplineCurve pSec1Curve;
	if(interpSec1.IsDone())
		pSec1Curve = interpSec1.Curve();
*/
	Standard_Real U1,U2;
	U1 = pSec1Curve->FirstParameter();
	U2 = pSec1Curve->LastParameter();

	gp_Pnt Pnt1;
	//check for symmetry
 	gp_Vec Vst,Vend;
	gp_Vec VecX = gp_Vec(gp_Pnt(0,0,0),gp_Pnt(1,0,0));
 	pSec1Curve->D1(U1,Pnt1,Vst);
	pSec1Curve->D1(U2,Pnt1,Vend);
 	double Ang1 = fabs(Vst.Angle(VecX))*180/M_PI;
	if(Ang1 > 90.0)
		Ang1 = 180.0-Ang1;
	double Ang2 = fabs(Vend.Angle(VecX))*180/M_PI;
	if(Ang2 > 90.0)
		Ang2 = 180.0-Ang2;


	Standard_Boolean isSymmetrical = fabs(Ang1-Ang2) < 1.0;

	//start test
	Pnt1 = pSec1Curve->Value(U1);
	gp_Pnt Pnt2 = pSec1Curve->Value(U1+0.03*(U2-U1));
	gp_Pnt2d Pnt1_2d(Pnt1.Y(),Pnt1.Z());
	gp_Pnt2d Pnt2_2d(Pnt2.Y(),Pnt2.Z());

	gp_Vec2d Vec1_2d(Pnt1_2d,Pnt2_2d);

	Standard_Real Ubreak1;
	for(i=5; i<100; i++)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/100.0;
		gp_Pnt Ptiter = pSec1Curve->Value(Uiter);
		gp_Pnt2d Pntiter_2d(Ptiter.Y(),Ptiter.Z());
		gp_Vec2d Veciter_2d(Pnt1_2d,Pntiter_2d);
		double Angle = Vec1_2d.Angle(Veciter_2d);

		if(fabs(Angle) > 1.0e-4)
		{
			Ubreak1 = U1+(i-1)*(U2-U1)/100.0;
			break;
		}
	}

	Pnt1 = pSec1Curve->Value(U2);
	Pnt2 = pSec1Curve->Value(U2-0.03*(U2-U1));
	Pnt1_2d = gp_Pnt2d(Pnt1.Y(),Pnt1.Z());
	Pnt2_2d = gp_Pnt2d(Pnt2.Y(),Pnt2.Z());

	Vec1_2d = gp_Vec2d(Pnt1_2d,Pnt2_2d);

	Standard_Real Ubreak2;
	for(i=95; i>0; i--)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/100.0;
		gp_Pnt Ptiter = pSec1Curve->Value(Uiter);
		gp_Pnt2d Pntiter_2d(Ptiter.Y(),Ptiter.Z());
		gp_Vec2d Veciter_2d(Pnt1_2d,Pntiter_2d);
		double Angle = Vec1_2d.Angle(Veciter_2d);

		if(fabs(Angle) > 1.0e-4)
		{
			Ubreak2 = U1+(i+1)*(U2-U1)/100.0;
			break;
		}
	}



	


	// end test
	
	////start from middle and go towards each end to find the U param of straight line segment
/*	Standard_Real Ubreak1;

	double angPrev = 0.0;
	gp_Vec VecPrev(gp_Pnt(0.0,0.0,0.0),gp_Pnt(0.0,0.0,1.0));
	for(i=100; i> 0; i--)
	{
		Standard_Real Uint = (U1+0.4*(U2-U1)) -(100-i)*(U2-U1)/100.0;
		gp_Vec V2;
		pSec1Curve->D1(Uint,Pnt1,V2);
		
		double ang = fabs(VecPrev.Angle(V2));
		if(i < 95 && ang < angPrev)
		{
			Ubreak1 = Uint;
			break;
		}
		else
		{
			angPrev = ang;
			VecPrev = V2;
		}
	}
*/
//	Standard_Real Ubreak2 = U2-(Ubreak1-U1);

/*	if(isSymmetrical)
		Ubreak2 = U2-(Ubreak1-U1);
	else
	{
		for(i=0; i<100; i++)
		{
			Standard_Real Uint = (U1+0.6*(U2-U1)) + i*(U2-U1)/100.0;
			gp_Vec V2;
			pSec1Curve->D1(Uint,Pnt1,V2);
			
			double ang = fabs(VecPrev.Angle(V2));
			if(i > 5 && ang < angPrev)
			{
				Ubreak2 = Uint;
				break;
			}
			else
			{
				angPrev = ang;
				VecPrev = V2;
			}
		}
	}
*/

	//sometimes the Ubreak1 or Ubreak2 seems to locate points too early. therefore get the max of both
	pSec1Curve->Segment(Ubreak1,Ubreak2);

/*	double gap = max((Ubreak1-U1),(U2-Ubreak2));
	pSec1Curve->Segment(U1+gap,U2-gap);
*/	

/*	//make the curve from end section
	Handle_TColgp_HArray1OfPnt pSecLastPoints = new TColgp_HArray1OfPnt(1,scrollEndSect.npt);
	i=0;
	index = 1;
	for(i=0; i< scrollEndSect.npt; i++)
	{
		gp_Pnt tempPt(scrollEndSect.x[i],scrollEndSect.y[i],scrollEndSect.z[i]);
		pSecLastPoints->SetValue(index,tempPt);
		index++;
	}
			
					
	GeomAPI_Interpolate interpSecLast(pSecLastPoints,Standard_False,1.0e-6);
	interpSecLast.Perform();
	Handle_Geom_BSplineCurve pSecLastCurve;
	if(interpSecLast.IsDone())
		pSecLastCurve = interpSecLast.Curve();
*/
	ShapeAnalysis_Curve SAC;
	gp_Pnt projPnt1,projPnt2;
	Standard_Real Uparam1,Uparam2;
	double dist = SAC.Project(pSecLastCurve,pSec1Curve->StartPoint(),1.0e-8,projPnt1,Uparam1);
	dist = SAC.Project(pSecLastCurve,pSec1Curve->EndPoint(),1.0e-8,projPnt2,Uparam2);

	Standard_Real Ua = min(Uparam1,Uparam2);
	Standard_Real Ub = std::max(Uparam1,Uparam2);
	Standard_Real Ugap = (Ub-Ua)/10.0;

	pSecLastCurve->Segment(Ua,Ub);
	//re assign the exact start and end point on pSecLastCurve
	projPnt1 = pSecLastCurve->StartPoint();
	projPnt2 = pSecLastCurve->EndPoint();

	//now make the bottom segment (we have to pick one of the two choices below)

	//choice 1 (circular arc 1 which goes tangent to the top segment)
	gp_Vec dirVec;
	gp_Pnt tempPt;
	pSecLastCurve->D1(Ua,tempPt,dirVec);
	dirVec.Reverse();
	GC_MakeArcOfCircle mkCirc(projPnt1,dirVec,projPnt2);
	
	//get the mid point of the bottom arc
	Handle_Geom_Curve pbottomArc = mkCirc.Value();
	gp_Pnt PtmidArc = pbottomArc->Value(0.5*(pbottomArc->FirstParameter()+pbottomArc->LastParameter()));

 	Handle_Geom_Curve pSelectCrv = mkCirc.Value();

/*
 	//choice 2 (circular arc 2 which goes tangent to the bottom level)
	gp_Pnt Pt1bottom = projPnt1;
	gp_Pnt Pt2bottom = projPnt2;
	BRepBuilderAPI_MakeEdge mkEdge(Pt1bottom,Pt2bottom);
	TopoDS_Edge bottomEdge = mkEdge.Edge();
	Standard_Real U1bottom,U2bottom;
	Handle_Geom_Curve pbottomCrv = BRep_Tool::Curve(bottomEdge,U1bottom,U2bottom);
	gp_Pnt Ptmidbottom = pbottomCrv->Value(0.5*(U1bottom+U2bottom));
	
	Handle_Geom_Curve pSelectCrv;
 	
	

 
 	if((Ptmidbottom.X()*Ptmidbottom.X()+Ptmidbottom.Y()*Ptmidbottom.Y()) > PtmidArc.X()*PtmidArc.X()+PtmidArc.Y()*PtmidArc.Y() ) //arc is too low
	{
		GC_MakeArcOfCircle mkCirc2(projPnt1,Ptmidbottom,projPnt2);
		pSelectCrv = mkCirc2.Value();
	}
	else
	{
		pSelectCrv = mkCirc.Value();
	}
 */
	//take parts off of the top curve to add to the bottom curve and make one single bottom curve
	Handle_Geom_BSplineCurve pPartACrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy());  
	Handle_Geom_BSplineCurve pPartBCrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy());     
	Handle_Geom_BSplineCurve pPartCCrv = Handle(Geom_BSplineCurve)::DownCast(pSecLastCurve->Copy()); 

	pPartACrv->Segment(Ua,Ua+Ugap);
	pPartBCrv->Segment(Ub-Ugap,Ub);
	pPartCCrv->Segment(Ua+Ugap,Ub-Ugap); 

	BRepBuilderAPI_MakeWire mkBottomWire(BRepBuilderAPI_MakeEdge(pSelectCrv).Edge());
	mkBottomWire.Add(BRepBuilderAPI_MakeEdge(pPartACrv).Edge());
	Standard_Boolean isdone = mkBottomWire.IsDone();
	mkBottomWire.Add(BRepBuilderAPI_MakeEdge(pPartBCrv).Edge());
	isdone = mkBottomWire.IsDone();
	if(!isdone)
	{
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSelectCrv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pPartACrv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pPartBCrv).Edge());
		return TopoDS_Wire();
	}

	Handle_Geom_BSplineCurve pBottomCrvFinal = makeSingleCurve(mkBottomWire.Wire());

	BRepBuilderAPI_MakeWire mkStartSectWire(BRepBuilderAPI_MakeEdge(pBottomCrvFinal).Edge());
	mkStartSectWire.Add(BRepBuilderAPI_MakeEdge(pPartCCrv).Edge());

	
	return mkStartSectWire.Wire();
}

TopoDS_Wire OCCVolute::makeExitEndSectionExp(TopoDS_Shape exitPipeExp)
{
	//extract the faces
	BRepBuilderAPI_MakeWire mkWire1;
	for(TopExp_Explorer explr(exitPipeExp,TopAbs_FACE); explr.More();explr.Next())
	{
		TopoDS_Face tempFace = TopoDS::Face(explr.Current());
		Handle_Geom_Surface pTempSurf = BRep_Tool::Surface(tempFace);
		Standard_Real U1,U2,V1,V2;
		pTempSurf->Bounds(U1,U2,V1,V2);

		Handle_Geom_Curve ptempCrv = pTempSurf->VIso(V1);
		mkWire1.Add(BRepBuilderAPI_MakeEdge(ptempCrv).Edge());
	}

	if(mkWire1.IsDone())
	{
		TopoDS_Wire endSectWire = mkWire1.Wire();
		return endSectWire;
	}

	return TopoDS_Wire();




}

TopoDS_Shape OCCVolute::makeBooleanExp(TopoDS_Face scroll1Face, TopoDS_Face trans1Face, TopoDS_Wire planeWire)
{
	//extract geom faces
	TopExp_Explorer Ex1;
    TopoDS_Face face1 = scroll1Face;
	Standard_Real U1,U2,V1,V2;
	int i=0;

	

    Handle_Geom_Surface pSurf1 = BRep_Tool::Surface(face1);
	pSurf1->Bounds(U1,U2,V1,V2);


	///Transition segment surface
    TopoDS_Face face2 = trans1Face;
	Standard_Real U1a,U2a,V1a,V2a;

    Handle_Geom_Surface pSurf2 = BRep_Tool::Surface(face2);
	pSurf2->Bounds(U1a,U2a,V1a,V2a);
	
	ShapeFix_Face Fixface1(face1);
	Fixface1.Perform();
	ShapeFix_Face Fixface2(face2);
	Fixface2.Perform();
	face1 = Fixface1.Face();
	face2 = Fixface2.Face();
//	debugList.push_back(face1);
//	debugList.push_back(face2);
	
	TopoDS_Shape sh1;
	try
	{
		BRepAlgoAPI_Section makeBool1(face1,face2 );
		makeBool1.Build();
		Standard_Boolean isDone = makeBool1.IsDone();
		if(isDone)
			sh1 = makeBool1.Shape();
		else
			return TopoDS_Shape();
	}
	catch(...)
	{
		debugList.push_back(face1);
		debugList.push_back(face2);
		return TopoDS_Shape();
	}


	//pick the intersection curve (select the longest edge if there is more than one
	TopoDS_Edge cutEdge;
	i=0;
	for (TopExp_Explorer ex(sh1,TopAbs_EDGE); ex.More(); ex.Next())
	{
		TopoDS_Edge aEdge = TopoDS::Edge(ex.Current()); 
		if(i==0)
			cutEdge = aEdge;
		if(i > 0)
		{
			GProp_GProps LProps;
			BRepGProp::LinearProperties(aEdge,LProps);
			double length1 = LProps.Mass();

			BRepGProp::LinearProperties(cutEdge,LProps);
			double length2 = LProps.Mass();

			if(length1 > length2)
				cutEdge = aEdge;
		}
		i++;
	}

	if(cutEdge.IsNull())
	{
		debugList.push_back(face1);
		debugList.push_back(face2);
		debugList.push_back(sh1);
		return TopoDS_Shape();
	}

	Standard_Real Ut1,Ut2;
	BRep_Tool::Range(cutEdge,Ut1,Ut2);

	Handle_Geom_BSplineCurve ptempCrv = makeSingleCurve(BRepBuilderAPI_MakeWire(cutEdge).Wire());
	Ut1 = ptempCrv->FirstParameter();
	Ut2 = ptempCrv->LastParameter();


	cutEdge = BRepBuilderAPI_MakeEdge(ptempCrv).Edge();
	BRep_Tool::Range(cutEdge,Ut1,Ut2);



	//make the cut face on scroll
	Standard_Real Ua1,Ua2;
	Handle_Geom_Curve pCutCrv = BRep_Tool::Curve(cutEdge,Ua1,Ua2);
	gp_Pnt Pta1 = pCutCrv->Value(Ua1);
	gp_Pnt Pta2 = pCutCrv->Value(Ua2);

	//get the first curve on pSurf1
	Standard_Real Usf1,Usf2,Vsf1,Vsf2;
	Usf1 = Usf2 = Vsf1 = Vsf2 = 0.0;
	BRepAdaptor_Surface BRAface1(face1);
	Extrema_ExtPS EEPS(Pta1,BRAface1,1.0e-8,1.0e-8);
	Standard_Integer numXPts = 0;
	if(EEPS.IsDone())
		numXPts = EEPS.NbExt();

	double sqDist = 100.0;
	int solLocator = 0;
	for(i=1; i<= numXPts; i++)
	{
		double sqDistNow = EEPS.SquareDistance(i);
		if(sqDist > sqDistNow)
		{
			sqDist = sqDistNow;
			solLocator = i;
		}
	}

	Extrema_POnSurf Ps = EEPS.Point(solLocator);
	Ps.Parameter(Usf1,Vsf1);



	EEPS.Perform(Pta2);
	if(EEPS.IsDone())
		numXPts = EEPS.NbExt();
	
	sqDist = 100.0;
	for(i=1; i<= numXPts; i++)
	{
		double sqDistNow = EEPS.SquareDistance(i);
		if(sqDist > sqDistNow)
		{
			sqDist = sqDistNow;
			solLocator = i;
		}
	}

	Ps = EEPS.Point(solLocator);
	Ps.Parameter(Usf2,Vsf2);

 	if((Usf1-U1) < 1.0e-3 || (U2-Usf2) < 1.0e-3)
 		return TopoDS_Shape();


	//get the curve segments on Psurf1
	Handle_Geom_BSplineCurve pSideCrvA = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(Vsf1));
	Handle_Geom_BSplineCurve pSideCrvB = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(Vsf2));

	Handle_Geom_BSplineCurve pSideCrvC = Handle(Geom_BSplineCurve)::DownCast(pSurf1->UIso(U1));
	Handle_Geom_BSplineCurve pSideCrvD = Handle(Geom_BSplineCurve)::DownCast(pSurf1->UIso(U2));

	if(Usf1 < Usf2)
	{
		pSideCrvA->Segment(U1,Usf1);
		pSideCrvB->Segment(Usf2,U2);
		pSideCrvC->Segment(Vsf1,V2);
		pSideCrvD->Segment(Vsf2,V2);
	}
	else
	{
		pSideCrvA->Segment(Usf1,U2);
		pSideCrvB->Segment(U1,Usf2);
		pSideCrvC->Segment(Vsf2,V2);
		pSideCrvD->Segment(Vsf1,V2);
	//	pSideCrvA->Segment(U1,Usf2);
	//	pSideCrvB->Segment(Usf1,U2);
	}
/*	Handle_Geom_BSplineCurve pSideCrvC = Handle(Geom_BSplineCurve)::DownCast(pSurf1->UIso(U1));
	pSideCrvC->Segment(Vsf1,V2);
	Handle_Geom_BSplineCurve pSideCrvD = Handle(Geom_BSplineCurve)::DownCast(pSurf1->UIso(U2));
	pSideCrvD->Segment(Vsf2,V2); */
	Handle_Geom_BSplineCurve pSideCrvE = Handle(Geom_BSplineCurve)::DownCast(pSurf1->VIso(V2));

	BRepBuilderAPI_MakeWire mkWire( cutEdge );
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSideCrvA).Edge());
	Standard_Boolean success = mkWire.IsDone();
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSideCrvB).Edge());
	success = mkWire.IsDone();
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSideCrvC).Edge());
	success = mkWire.IsDone();
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSideCrvD).Edge());
	success = mkWire.IsDone();
	mkWire.Add(BRepBuilderAPI_MakeEdge(pSideCrvE).Edge());
	success = mkWire.IsDone();

	if(!mkWire.IsDone())
	{
		debugList.push_back(cutEdge);
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSideCrvA).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSideCrvB).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSideCrvC).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSideCrvD).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSideCrvE).Edge());
		debugList.push_back(face1);
		debugList.push_back(face2);
		return TopoDS_Shape();
	}

	BRepBuilderAPI_MakeFace mkFace1(pSurf1,mkWire.Wire());
	success = mkFace1.IsDone();

	if(success)
	{
		TopoDS_Face tempface = mkFace1.Face();
		ShapeFix_Face fixFace1(tempface);
		fixFace1.Perform();
		Standard_Boolean status = fixFace1.Status(ShapeExtend_OK);
		trimFace1Exp = fixFace1.Face();
	}
	else
		return TopoDS_Shape();


 
	//make the face on transition segment	
	Handle_Geom_Curve pVIsoCrv = pSurf2->VIso(V1a);

	/////////////////
	BRepAdaptor_Surface BRAface2(face2);
	Extrema_ExtPS EEPS2(Pta1,BRAface2,1.0e-8,1.0e-8);
	numXPts = 0;
	if(EEPS2.IsDone())
		numXPts = EEPS2.NbExt();

	sqDist = 100.0;
	solLocator = 0;
	for(i=1; i<= numXPts; i++)
	{
		double sqDistNow = EEPS2.SquareDistance(i);
		if(sqDist > sqDistNow)
		{
			sqDist = sqDistNow;
			solLocator = i;
		}
	}

	Ps = EEPS2.Point(solLocator);
	Ps.Parameter(Usf1,Vsf1);



	EEPS2.Perform(Pta2);
	if(EEPS2.IsDone())
		numXPts = EEPS2.NbExt();
	
	sqDist = 100.0;
	for(i=1; i<= numXPts; i++)
	{
		double sqDistNow = EEPS2.SquareDistance(i);
		if(sqDist > sqDistNow)
		{
			sqDist = sqDistNow;
			solLocator = i;
		}
	}

	Ps = EEPS2.Point(solLocator);
	Ps.Parameter(Usf2,Vsf2);
	Handle_Geom_BSplineCurve pCrv1 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(Vsf1)); 
	Handle_Geom_BSplineCurve pCrv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(Vsf2)); 
	if(Usf1 < Usf2)
	{
		pCrv1->Segment(U1a,Usf1);
		pCrv2->Segment(Usf2,U2a);
	}
	else
	{
		pCrv1->Segment(Usf1,U2a);
		pCrv2->Segment(U1a,Usf2);
//		pCrv1->Segment(U1a,Usf2);
//		pCrv2->Segment(Usf1,U2a);
	}



	/////////////////

/*	ShapeAnalysis_Curve SAC;
	gp_Pnt projPnt1,projPnt2;
	Standard_Real Uparam1,Uparam2;
	double dist = SAC.Project(pSurf2->VIso(V1a),Pta1,1.0e-5,projPnt1,Uparam1);
	dist = SAC.Project(pSurf2->VIso(V1a),Pta2,1.0e-5,projPnt2,Uparam2);

	Handle_Geom_BSplineCurve pCrv1 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv1->Segment(U1a,min(Uparam1,Uparam2)); */
	TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(pCrv1).Edge();

/*	Handle_Geom_BSplineCurve pCrv2 = Handle(Geom_BSplineCurve)::DownCast(pSurf2->VIso(V1a)); 
	pCrv2->Segment(max(Uparam1,Uparam2),U2a); */
	TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(pCrv2).Edge();



	//add the balance edges to close the wire on exit pipe surface
	Handle_Geom_BSplineCurve pSide1Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U1a)->Copy());
	Handle_Geom_BSplineCurve pSide2Crv = Handle(Geom_BSplineCurve)::DownCast(pSurf2->UIso(U2a)->Copy());
	if(Usf1 < Usf2)
	{
		pSide1Crv->Segment(Vsf1,V2a);
		pSide2Crv->Segment(Vsf2,V2a);
	}
	else
	{
		pSide1Crv->Segment(Vsf2,V2a);
		pSide2Crv->Segment(Vsf1,V2a);
	}


	BRepBuilderAPI_MakeWire mkSegWire;
	//main wire
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(ptempCrv).Edge());
	mkSegWire.Add(edge1/*BRepBuilderAPI_MakeEdge(pSegCrv1).Edge()*/);
	Standard_Boolean isdone1 = mkSegWire.IsDone();	
	mkSegWire.Add(edge2/*BRepBuilderAPI_MakeEdge(pSegCrv2).Edge()*/);
	Standard_Boolean isdone2 = mkSegWire.IsDone();

	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSide1Crv).Edge());
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSide2Crv).Edge());
	mkSegWire.Add(BRepBuilderAPI_MakeEdge(pSurf2->VIso(V2a)).Edge());

	if(!isdone1 || !isdone2 || !mkSegWire.IsDone())
	{
		debugList.push_back(BRepBuilderAPI_MakeEdge(ptempCrv).Edge());
		debugList.push_back(edge1/*BRepBuilderAPI_MakeEdge(pSegCrv1).Edge()*/);
		debugList.push_back(edge2/*BRepBuilderAPI_MakeEdge(pSegCrv2).Edge()*/);
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSide1Crv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSide2Crv).Edge());
		debugList.push_back(BRepBuilderAPI_MakeEdge(pSurf2->VIso(V2a)).Edge());
		return TopoDS_Shape();
	}

	BRepBuilderAPI_MakeFace mkFace2(pSurf2,mkSegWire.Wire());	
	if(mkFace2.IsDone())
	{
		TopoDS_Face tempface = mkFace2.Face();
		ShapeFix_Face fixFace2(tempface);
		fixFace2.Perform();
		Standard_Boolean status = fixFace2.Status(ShapeExtend_OK);
		trimFace2Exp = fixFace2.Face();
	}
	else
		return TopoDS_Shape();


	//make the fillet
	TopoDS_Shape filletedShapeExp;
	bool isFilleted = false;
	double Radmin = 0.0075;
	double Radmax = 0.0075;
 
/*	for(int i=0; i<10; i++)
	{
		isFilleted = makeTongueFillet(TopoDS::Face(trimFace1Exp), TopoDS::Face(trimFace2Exp), 1,Radmin, Radmax, Radmin, filletedShapeExp, 0);
		
		Radmin = Radmin * 0.75;
		Radmax = Radmax * 0.75;

		if(isFilleted)
		{
		 	ShapeFix_Shape fixShape(filletedShapeExp);
			fixShape.Perform();
			filletedShapeExp = fixShape.Shape();  


			ShapeAnalysis_ShapeContents SAS;
			SAS.Perform(filletedShapeExp);
			int freeWires = SAS.NbFreeWires();
			int numFaces = SAS.NbFaces();
			int freeFaces = SAS.NbFreeFaces();
			int freeEdges = SAS.NbFreeEdges();

			if(numFaces != 3)
			{
				isFilleted = false;
			}
			else
			{
				break;
			}
		}
	}
*/

//	debugList.clear();
	
	double tongueCurv=0.8;
	double tongueSize=0.2;
	
	
	
	cutEdge=BRepBuilderAPI_MakeEdge(ptempCrv).Edge();
	
	InterSecEdge=cutEdge;

	
		
	isFilleted =makeTongueFilletUsingThrughSections(TopoDS::Face(trimFace1Exp), TopoDS::Face(trimFace2Exp),tongueCurv,tongueSize,filletedShapeExp);
	
	if (isFilleted)return filletedShapeExp;
//	isFilleted = false;
	if(!isFilleted)
	{
		BRepBuilderAPI_Sewing sewing(0.005);
		sewing.Add(trimFace1Exp);
		sewing.Add(trimFace2Exp);
		sewing.Perform();
		filletedShapeExp = sewing.SewedShape();
	}

	
	return filletedShapeExp;
}


Handle_Geom_Curve OCCVolute::makeProjCurve(Handle_Geom_Curve pCrv, double R)
{
	double radSq = R*R;
	Handle_Geom_BSplineCurve pBSPCrv;

	Standard_Real U1,U2;
	U1 = pCrv->FirstParameter();
	U2 = pCrv->LastParameter();

	gp_Pnt Pt1;
	pCrv->D0(U1,Pt1);
//	gp_Vec Vec1 = gp_Vec(Pt1.X(),Pt1.Y(),0.0);
	gp_Ax1 Ax1 = gp::OZ();

	double baseAngle;// = atan(Pt1.Y()/Pt1.X());
	gp_Vec2d VecOX(gp_Pnt2d(0,0),gp_Pnt2d(0,1));
	gp_Vec2d VecOR(gp_Pnt2d(0,0),gp_Pnt2d(Pt1.X(),Pt1.Y()));
	baseAngle = VecOR.Angle(VecOX);
	double baseX = R*sin(baseAngle);
	double baseY = R*cos(baseAngle);

	//locate baseZ1 and baseZ2 points
	double baseZ1 = 0.0;
	double baseZ2 = 0.0;
	for(int i=0; i<1000; i++)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/1000.0;
		gp_Pnt tempPt = pCrv->Value(Uiter);
		double tempRad = sqrt(tempPt.X()*tempPt.X()+tempPt.Y()*tempPt.Y());
		if(fabs(R-tempRad) < 1.0e-5 || tempRad > R)
		{
			baseZ1 = tempPt.Z();
			break;
		}
	}
	for(int i=1000; i>0; i--)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/1000.0;
		gp_Pnt tempPt = pCrv->Value(Uiter);
		double tempRad = sqrt(tempPt.X()*tempPt.X()+tempPt.Y()*tempPt.Y());
		if(fabs(R-tempRad) < 1.0e-5 || tempRad > R)
		{
			baseZ2 = tempPt.Z();
			break;
		}
	}

	gp_Vec2d Vec1(gp_Pnt2d(0.0,0.0),gp_Pnt2d(Pt1.X(),Pt1.Y()));

	int index = 50;
	Handle_TColgp_HArray1OfPnt pPoints = new TColgp_HArray1OfPnt(1,index+1);

	int counter = 0;
	bool lowSet = false;
	bool highSet = false;
	for(int i=0; i<=50; i++)
	{
		Standard_Real Uiter = U1+i*(U2-U1)/(double(index));
		gp_Pnt Ptiter = pCrv->Value(Uiter);

		//gp_Vec Veciter = gp_Vec(Ptiter.X(),Ptiter.Y(),0.0);
		gp_Vec2d Veciter(gp_Pnt2d(0.0,0.0),gp_Pnt2d(Ptiter.X(),Ptiter.Y()));

		double rotAng = Veciter.Angle(Vec1);

		Ptiter.Rotate(Ax1,rotAng) ;
		//start test
		double tempRadSq = Ptiter.X()*Ptiter.X()+Ptiter.Y()*Ptiter.Y();
		if(tempRadSq < radSq)
		{
			if(i < 25 && !lowSet)
			{
				pPoints->SetValue(counter+1,gp_Pnt(baseX,baseY,baseZ1));
				lowSet = true;
				counter++;
			}
			else if(i > 25 && !highSet)
			{
				pPoints->SetValue(counter+1,gp_Pnt(baseX,baseY,baseZ2));
				highSet = true;
				counter++;
			}
		}
		else if(tempRadSq > radSq)
		{
			pPoints->SetValue(counter+1,Ptiter);
			counter++;
		}

		//end test
//		pPoints->SetValue(i+1,Ptiter);

	}

	//check the number of separate points
	Handle_TColgp_HArray1OfPnt pPoints2 = new TColgp_HArray1OfPnt(1,counter);
	for(int i=1; i<= counter; i++)
	{
		pPoints2->SetValue(i,pPoints->Value(i));
	}
	

	GeomAPI_Interpolate interp(pPoints2,Standard_False,1.0e-8);
    interp.Perform();
    if(interp.IsDone())
            pBSPCrv = interp.Curve();

	return pBSPCrv;


}


TopoDS_Shape OCCVolute::makeFlangePipe(TopoDS_Shape exitPipeExp, double Offset1)
{
	TopoDS_Shape flangePipe;

	//flange outer radius (taken as twice the exit pipe radius)
	TopoDS_Face endCap = TopoDS::Face(exitPipeFaces[exitPipeFaces.size()-1]);
	GProp_GProps SProps;
    BRepGProp::SurfaceProperties(endCap,SProps);
    gp_Pnt SectionCentroid = SProps.CentreOfMass();
	double exitPipeArea = SProps.Mass();
	double exitPipeRadius = sqrt(exitPipeArea/M_PI);

	double scale = 2.0;

	TopoDS_Wire endCapWire;
	for(TopExp_Explorer explr(endCap,TopAbs_WIRE); explr.More(); explr.Next())
	{
		endCapWire = TopoDS::Wire(explr.Current());
	}


	Handle_Geom_BSplineCurve  pcrv1 = makeSingleCurve(endCapWire/*exitTranSectWires[2]*/);
	Handle_Geom_BSplineCurve  pcrv2 = Handle(Geom_BSplineCurve)::DownCast( pcrv1->Copy());
	pcrv2->Scale(SectionCentroid,scale);
	Handle_Geom_Surface pFillSurf = GeomFill::Surface(pcrv1,pcrv2) ;

	TopoDS_Shape flangeFace1 = BRepBuilderAPI_MakeFace(pFillSurf,1.0e-6).Face();

//Correction for Spliter finish at the pipe end	
	if ((hasSplitter)&&(!Spliterendtrue))
	flangeFace1=makePipeFlangeouterFacePVMV();

//Correction for Spliter finish at the pipe end

	//trim a length equal to flange thickness off the expanded pipe
	//flange thickness is set to twice the wall thickness
	BRepBuilderAPI_Sewing exitPipeExpTrimSewing(1.0e-06,Standard_True,Standard_True,Standard_True,Standard_True);
	BRepBuilderAPI_MakeWire mkPlnWire1; //mkExitPlnWire;
	BRepBuilderAPI_MakeWire mkPlnWire2; //mkInletPlnWire;
	BRepBuilderAPI_MakeWire mkPlnWire3;
	BRepBuilderAPI_MakeWire mkPlnWire4;

	Handle_Geom_Curve pCrv4;
	int i = 0;
	BRepOffsetAPI_ThruSections thru1;
	for(TopExp_Explorer explr(exitPipeExp,TopAbs_FACE);explr.More();explr.Next())
	{
		TopoDS_Face exitFace = TopoDS::Face(explr.Current());
		Handle_Geom_Surface pexitSurfExp = BRep_Tool::Surface(exitFace);
		Standard_Real U1b,U2b,V1b,V2b;
		pexitSurfExp->Bounds(U1b,U2b,V1b,V2b);

		Handle_Geom_Curve pIsoCrv1 = pexitSurfExp->UIso(U1b);
		GeomAdaptor_Curve     GAC(pIsoCrv1 );
		double len = GCPnts_AbscissaPoint::Length( GAC);
		double Vtrim = V2b-(2.0*Offset1/len)*(V2b-V1b);
		if(Vtrim < V1b || Vtrim > V2b)
			Vtrim = V2b-0.95*(V2b-V1b);

		//make trimmed surface of exit pipe exp
		pCrv4 = pexitSurfExp->VIso(Vtrim); // this is declared outside of the for loop for later use

		
		Handle_Geom_Curve pCrva = pexitSurfExp->VIso(V1b);
		Handle_Geom_Curve pCrvb = pexitSurfExp->VIso(V1b+0.3*(V2b-V1b));
		Handle_Geom_Curve pCrvc = pexitSurfExp->VIso(V1b+0.6*(V2b-V1b));
		Handle_Geom_Curve pCrvd = pexitSurfExp->VIso(Vtrim);
		
		mkPlnWire1.Add(BRepBuilderAPI_MakeEdge(pCrva).Edge());
		mkPlnWire2.Add(BRepBuilderAPI_MakeEdge(pCrvb).Edge());
		mkPlnWire3.Add(BRepBuilderAPI_MakeEdge(pCrvc).Edge());
		mkPlnWire4.Add(BRepBuilderAPI_MakeEdge(pCrvd).Edge());
	}

	TopoDS_Wire exitPlnWire = mkPlnWire4.Wire();// mkExitPlnWire.Wire();
	BRepBuilderAPI_MakeFace mkexitPlnFace(exitPlnWire);
	mkexitPlnFace.Build();
	TopoDS_Face exitPlnFace;
	if(mkexitPlnFace.IsDone())
	{
		exitPlnFace = mkexitPlnFace.Face();
	}
	else
	{
		Handle_Geom_BSplineCurve pExitPlnCrv = makeSingleCurve(exitPlnWire);
		Standard_Real Uc1 = pExitPlnCrv->FirstParameter();
		Standard_Real Uc2 = pExitPlnCrv->LastParameter();
		gp_Pnt Ptc1 = pExitPlnCrv->Value(Uc1);
		gp_Pnt Ptc2 = pExitPlnCrv->Value(Uc1+0.3*(Uc2-Uc1));
		gp_Pnt Ptc3 = pExitPlnCrv->Value(Uc1+0.6*(Uc2-Uc1));

		gp_Circ circ = gce_MakeCirc(Ptc1,Ptc2,Ptc3).Value();
		BRepBuilderAPI_MakeFace mkcircFace(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire());
		if(mkcircFace.IsDone())
			exitPlnFace = mkcircFace.Face();
		else
			return TopoDS_Shape();
	}

	///////////////
//	TopoDS_Wire inletPlnWire = mkInletPlnWire.Wire();
	BRepOffsetAPI_ThruSections tempThruSect(Standard_False,Standard_False,1.0e-06);
	tempThruSect.AddWire(mkPlnWire1.Wire());
	tempThruSect.AddWire(mkPlnWire2.Wire());
	tempThruSect.AddWire(mkPlnWire3.Wire());
	tempThruSect.AddWire(mkPlnWire4.Wire());

	tempThruSect.Build();
	TopoDS_Shape tempFaceXX = tempThruSect.Shape();

	///////////////

	
    BRepGProp::SurfaceProperties(exitPlnFace,SProps);
    SectionCentroid = SProps.CentreOfMass();
	gp_Pnt Pnt1 = pCrv4->Value(pCrv4->FirstParameter());

	double outerPipeRadius = SectionCentroid.Distance(Pnt1);
	scale = 2.0*exitPipeRadius/outerPipeRadius;
	

	//make one flange surface
	pcrv1 = makeSingleCurve(mkPlnWire4.Wire()/*mkExitPlnWire.Wire()*/);
	Handle_Geom_BSplineCurve pcrv3 = Handle(Geom_BSplineCurve)::DownCast( pcrv1->Copy()); 
	pcrv3->Scale(SectionCentroid,scale);
	pFillSurf = GeomFill::Surface(pcrv1,pcrv3) ;
	TopoDS_Shape flangeFace2 = BRepBuilderAPI_MakeFace(pFillSurf,1.0e-6).Face();
	//split in to two
	Standard_Real U1f,U2f,V1f,V2f;
	pFillSurf->Bounds(U1f,U2f,V1f,V2f);
	Handle_Geom_BSplineSurface pFillSurfA = Handle(Geom_BSplineSurface)::DownCast(pFillSurf->Copy());
	Handle_Geom_BSplineSurface pFillSurfB = Handle(Geom_BSplineSurface)::DownCast(pFillSurf->Copy());
	pFillSurfA->Segment(U1f,U1f+0.5*(U2f-U1f),V1f,V2f);
	pFillSurfB->Segment(U1f+0.5*(U2f-U1f),U2f,V1f,V2f);
	TopoDS_Shape flangeFace2A = BRepBuilderAPI_MakeFace(pFillSurfA,1.0e-6).Face();
	TopoDS_Shape flangeFace2B = BRepBuilderAPI_MakeFace(pFillSurfB,1.0e-6).Face();


	BRepBuilderAPI_Sewing sewing1(1.0e-06,Standard_True,Standard_True,Standard_True,Standard_True);
	sewing1.Add(tempFaceXX/*exitPipeExpTrim*/);
	//sewing1.Add(flangeFace2);
	sewing1.Add(flangeFace2A);
	sewing1.Add(flangeFace2B);
	sewing1.Perform();
	TopoDS_Shape flangePart1 = sewing1.SewedShape();

	ShapeFix_Shape fixShape;

	


	

	//do the fillet
	BRepFilletAPI_MakeFillet mkFlangeFillet1(flangePart1);
	Standard_Real filletR1 = outerPipeRadius/2.0;
	
	i=0;
	int edNum1 = 2;
	TopoDS_Edge edgeToFillet;
	for(TopExp_Explorer explr(flangePart1,TopAbs_EDGE);explr.More();explr.Next())
	{
		edgeToFillet = TopoDS::Edge(explr.Current());
		if(i == edNum1 )
		{
			mkFlangeFillet1.Add(filletR1,edgeToFillet);
		}
		 	
		i++;
	}

	try
	{
		mkFlangeFillet1.Build();
	}
	catch(...)
	{
		return TopoDS_Shape();
	}

	TopoDS_Shape outerFlange;
	if(mkFlangeFillet1.IsDone())
	{
		outerFlange = mkFlangeFillet1.Shape();
	}
	else
	{
		outerFlange = flangePart1; //use the unfilleted part if fillet failed
	}
	
	
 

	//make the flange thickness face	
	//crv tests
	gp_Pnt Pt1 = pcrv2->Value(pcrv2->FirstParameter());
	gp_Pnt Pt2 = pcrv3->Value(pcrv3->FirstParameter());
	double dist = Pt1.Distance(Pt2);
	if(dist > 4.0*Offset1)
	{
		Standard_Real Ust = pcrv2->FirstParameter();
		Standard_Real Uend = pcrv2->LastParameter();
		
		Handle_Geom_BSplineCurve ptempCrv1 = Handle(Geom_BSplineCurve)::DownCast(pcrv2->Copy());
		Handle_Geom_BSplineCurve ptempCrv2 = Handle(Geom_BSplineCurve)::DownCast(pcrv2->Copy());
		ptempCrv1->Segment(Ust,Ust+0.5*(Uend-Ust));
		ptempCrv2->Segment(Ust+0.5*(Uend-Ust),Uend);

		ptempCrv2 = makeSingleCurve(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(ptempCrv2).Edge()).Wire());
		BRepBuilderAPI_MakeWire mktempWire;
		mktempWire.Add(BRepBuilderAPI_MakeEdge(ptempCrv2).Edge());
		mktempWire.Add(BRepBuilderAPI_MakeEdge(ptempCrv1).Edge());
		pcrv2 = makeSingleCurve(mktempWire.Wire());	
	}
	//check whether the curves are in same orientation
	Pt1 = pcrv2->Value(pcrv2->FirstParameter()+0.25*(pcrv2->LastParameter()- pcrv2->FirstParameter()));
	Pt2 = pcrv3->Value(pcrv3->FirstParameter()+0.25*(pcrv3->LastParameter()- pcrv3->FirstParameter()));
	dist = Pt1.Distance(Pt2);
	if(dist > 4.0*Offset1)
		pcrv2->Reverse();


	pFillSurf = GeomFill::Surface(pcrv2,pcrv3) ;

	TopoDS_Shape flangeThickFace = BRepBuilderAPI_MakeFace(pFillSurf,1.0e-6).Face();
	fixShape.Init(flangeThickFace);
	fixShape.Perform();
	flangeThickFace = fixShape.Shape();

	//sew all the faced together
	BRepBuilderAPI_Sewing flangeSewing(1.0e-6);
	flangeSewing.Add(flangeFace1);
	flangeSewing.Add(flangeThickFace);
	flangeSewing.Add(/* flangePart1*/outerFlange );

	flangeSewing.Perform();
	flangePipe = flangeSewing.SewedShape();
	fixShape.Init(flangePipe);
	fixShape.Perform();
	flangePipe = fixShape.Shape();

	return flangePipe;
}



void OCCVolute::writeStepFile()
{ 
	//test work for ctaads
	Standard_CString aFileName = (Standard_CString) "rotorBlade.stp";
	STEPControl_Reader aReader;
	IFSelect_ReturnStatus status1 = aReader.ReadFile(aFileName);
	if ( status1 == IFSelect_RetDone )
    {
	    TopoDS_Compound comp;
		TopoDS_Shape tempShape;

		aReader.TransferRoots();
		int n = aReader.NbShapes();
		for (int i = 1; i <= 1/*n*/; ++i) 
		{
			tempShape = aReader.Shape( i ); 
			BRepTools::Write(tempShape,"tempShape.brep");
		}

		//extract what we want
		TopoDS_Face mainFace;
		Handle_Geom_Surface pmainFaceSurf;
		for(TopExp_Explorer explr(tempShape,TopAbs_FACE); explr.More(); explr.Next())
		{
			mainFace = TopoDS::Face(explr.Current());
			pmainFaceSurf = BRep_Tool::Surface(mainFace);
			BRepTools::Write(mainFace,"mainFace.brep");
		}

		Standard_Real U1,U2,V1,V2;
		pmainFaceSurf->Bounds(U1,U2,V1,V2);
		gp_Pnt startPt = pmainFaceSurf->Value(U1,V1);
		gp_Pnt endPt = pmainFaceSurf->Value(U2,V1);
		double radii_st = sqrt(startPt.Y()*startPt.Y()+startPt.Z()*startPt.Z());
		double radii_end = sqrt(endPt.Y()*endPt.Y()+endPt.Z()*endPt.Z());
		double diff = radii_st-radii_end;

		gp_Pnt midPt1 = pmainFaceSurf->Value(0.5*(U1+U2),V1);
		gp_Pnt midPt2 = pmainFaceSurf->Value(0.5*(U1+U2),V2);

		//translate the file to reqd radial height (255mm)
		gp_Trsf  t1;

		t1.SetTranslation(gp_Pnt(0,0,0),gp_Pnt(0,0,255));
		BRepBuilderAPI_Transform translate(mainFace,t1); 
		TopoDS_Face faceatTDC = TopoDS::Face(translate.Shape());
		BRepTools::Write(faceatTDC,"faceatTDC.brep");

		//write an iges file
		IGESControl_Controller::Init();
		IGESControl_Writer writer( Interface_Static::CVal( "XSTEP.iges.unit" ),
                                Interface_Static::SetIVal ("write.iges.brep.mode", 1)
                             /*  Interface_Static::IVal( "XSTEP.iges.writebrep.mode" ) */);

		 Standard_Integer byvalue = Interface_Static::IVal("write.iges.brep.mode");
		 Standard_Boolean ok = writer.AddShape ( faceatTDC );

		 writer.ComputeModel();
		 writer.Write( (Standard_CString)"rotorAtTDC.igs" );

	}




	return;
	//end test work for ctaads

	//writing a step file
    STEPControl_StepModelType type = STEPControl_AsIs;
    IFSelect_ReturnStatus status;
    STEPControl_Writer writer;

	for(int i=0; i<debugList.size(); i++)
	{
		if(!debugList[i].IsNull())
			status = writer.Transfer(debugList[i] , type );
	}




  	if(!totalVoluteThick.IsNull())
	{
		status = writer.Transfer(totalVoluteThick, type);
	}
	else if( !totalVolute.IsNull())
	{
 		status = writer.Transfer(totalVolute, type);
	}
	else
	{
 		if(!scrollShape1.IsNull())
 			status = writer.Transfer( scrollShape1 , type ); 
		for(int i=0; i<scrollFaces.size();i++)
			status = writer.Transfer(scrollFaces[i],type);
	  	if(!scrollShape2.IsNull())
  			status = writer.Transfer( scrollShape2 , type );
 		
		if(!transPipeShape1.IsNull())
			status = writer.Transfer( transPipeShape1, type);

 		if(!exitPipeShape1.IsNull())
			status = writer.Transfer( exitPipeShape1, type);    
	}

 

    Standard_CString filename = "test.step";
    status = writer.Write( filename );
	//end writing a step file


/*	//write a stl file
	if(!totalVolute.IsNull())
	{
		filename = "testSTL.stl";
 		StlAPI_Writer writerSTL;
		writerSTL.RelativeMode() = Standard_False;

		Standard_Real aDeflection = 0.001;
		Standard_Real aCoefficient = 0.001;
		writerSTL.SetDeflection(aDeflection);
		writerSTL.SetCoefficient(aCoefficient);
 		writerSTL.Write(totalVolute, filename);
	} */

	//write a brep file
	if(!totalVolute.IsNull())
		BRepTools::Write(totalVolute,"test.brep");


	return;
}



int OCCVolute::makeThickVolute(double wallThickness, bool isTurbine, bool exit_straight, bool isClockWise)
{
	int errorCode = 0;

	//make the outer shell for thick solid
	errorCode = makeOuterShell(wallThickness, isTurbine, exit_straight, isClockWise);

 	 //writeStepFile();

	return errorCode;
}

TopoDS_Shape OCCVolute::makeInnerShellNoSplitter(Handle_Geom_Curve pScrollSectLastExp)
{
	BRepBuilderAPI_Sewing innerShellSewing( 0.005 /*1.0e-4*/);
	for(int i=0; i<scrollFaces.size(); i++)
	{
		for(TopExp_Explorer explr(scrollFaces[i],TopAbs_FACE); explr.More(); explr.Next())
		{
			TopoDS_Shape innerShellFace = TopoDS::Face(explr.Current());
			innerShellSewing.Add(innerShellFace);
			break;
		}
	}
	if(filletedShapeSaved.IsNull())
	{
		innerShellSewing.Add(trimFace1);
		innerShellSewing.Add(trimFace2);
	}
	else
		innerShellSewing.Add(filletedShapeSaved);

	for(int i=0; i< tongueRegionPipeFaces.size();i++)
		innerShellSewing.Add(tongueRegionPipeFaces[i]);
	for(int i=0; i<exitPipeFaces.size()-1; i++)
		innerShellSewing.Add(exitPipeFaces[i]);

	
	//make the base surfaces
	gp_Pnt PtOut1 = pScrollSectLastExp->Value(pScrollSectLastExp->FirstParameter());
	gp_Pnt PtOut2 = pScrollSectLastExp->Value(pScrollSectLastExp->LastParameter());

	gp_Pnt PtIn1, PtIn2;

	TopoDS_Wire lastScrollSect = seg0Wires[seg0Wires.size()-1];
	for(TopExp_Explorer explr(lastScrollSect,TopAbs_EDGE); explr.More(); explr.Next())
	{
		TopoDS_Edge tEdge = TopoDS::Edge(explr.Current());
		Standard_Real Uin1,Uin2;
		Handle_Geom_Curve tCrv = BRep_Tool::Curve(tEdge,Uin1,Uin2);

		PtIn1 = tCrv->Value(Uin1);
		PtIn2 = tCrv->Value(Uin2);
	}

	//make edge1
	TopoDS_Edge sweepEdge1 = BRepBuilderAPI_MakeEdge(PtOut1,PtIn1).Edge();
	TopoDS_Edge sweepEdge2 = BRepBuilderAPI_MakeEdge(PtIn2,PtOut2).Edge();



	gp_Ax1 axe = gp::OZ();
	TopoDS_Shape base1 = BRepPrimAPI_MakeRevol(sweepEdge1,axe); 
	TopoDS_Shape base2 = BRepPrimAPI_MakeRevol(sweepEdge2,axe); 


 	innerShellSewing.Add(base1);
	innerShellSewing.Add(base2);

	innerShellSewing.Perform();
	
	return innerShellSewing.SewedShape();




}

bool OCCVolute::separateFaces()
{
	//create separate parts
	TopoDS_Compound comp1,comp2,comp3;
	BRep_Builder builder1,builder2,builder3;
	builder1.MakeCompound( comp1 );//Part between in and exit
	builder2.MakeCompound( comp2 );//inlet
	builder3.MakeCompound( comp3 );//exit

	std::vector <TopoDS_Face>inputFace,exitFace;

	if(!m_inputPlane.IsNull())
	{
		for(TopExp_Explorer explr(m_inputPlane,TopAbs_FACE);explr.More(); explr.Next())
			{
				inputFace.push_back(TopoDS::Face(explr.Current()));
			}
	}

	if(!m_exitPlane.IsNull())
	{
		for(TopExp_Explorer explr(m_exitPlane,TopAbs_FACE);explr.More(); explr.Next())
		{
			exitFace.push_back(TopoDS::Face(explr.Current()));
		}
	}

	double area;
	gp_Pnt center;

	for(TopExp_Explorer explr(totalVolute,TopAbs_FACE);explr.More(); explr.Next())
	{
		TopoDS_Face aShape = TopoDS::Face(explr.Current());

		bool isFound = false;
		//First find area in inlet planes
		for(int i = 0; i < inputFace.size();i++)
		{
			if(isSameFace(inputFace[i], aShape))
			{
				builder2.Add(comp2,aShape);
				isFound = true;
				break;
			}
		}
		if(isFound) continue;

		//Then find area in exit planes
		for(int i = 0; i < exitFace.size();i++)
		{
			if(isSameFace(exitFace[i], aShape))
			{
				builder3.Add(comp3,aShape);
				isFound = true;
				break;
			}
		}
		if(isFound) continue;

		builder1.Add(comp1,aShape);
	}

	m_totalVoluteWithoutInOut = comp1;
	m_inputPlane = comp2;
	m_exitPlane = comp3;

	return true;
}

bool OCCVolute::makeTransitionPatch(bool firstRun, bool overRideTransPos,bool isClockWise,std::vector<asymSection> asymSects, double tongue_le_radius, double tongue_R_aspect_ratio, double tongue_z_aspect_ratio, int & TongueTransEndIndex, double TongueStAngle,
		double throat_area,double exit_area, double exit_length,
		bool exit_straight, double exit_pipe_angle_v, double exit_pipe_angle_h, int exit_curve_order, double* pX_exitCrv, double* pY_exitCrv, bool exit_enable_transition,
		double transitionLength, double exit_pipe_aspect1, double firstTransPos, int transSectOption, int firstTransPosOption, int& runFlag, double& calcedSecondTransPos,
		double earlyTurnAngle)
{
	//this method creates a transition patch between the scroll and the exit pipe
	double TongueR= tongue_le_radius;
	double AspectRatio = abs(tongue_R_aspect_ratio);
	const bool Tonguless = tongue_le_radius <=1.0e-6 ? true : false;
	bool IsPipeShell =!exit_straight || !exit_enable_transition;
	int NPnt = asymSects[0].x.size();   //number of points per sector
	int NumOfSect = asymSects.size(); // number of sectors
	double TongueStAng = TongueStAngle;  // in Degrees
	gp_Pnt PtSect_St = gp_Pnt(asymSects[0].x[0],asymSects[0].y[0],asymSects[0].z[0]); // Start point of Start Sector of the volute
	gp_Pnt PtSect_End = gp_Pnt(asymSects[0].x[NPnt-1],asymSects[0].y[NPnt-1],asymSects[0].z[NPnt-1]); // End point of Start Sector of the volute
	const double SectWidth = PtSect_St.Distance(PtSect_End);
	double TranslateFactor = tongue_z_aspect_ratio *SectWidth/3.0;  // according to +/- sign Tongue fillet will go forward or backward 
	
	if(tongue_le_radius <=1.0e-6)
	{
		TongueR = SectWidth/10.0;
	}
	
	if(IsPipeShell)
	{
		TongueStAng = 0.0;
		TranslateFactor = 0.0;
	}
	double tempStartAngle = TongueStAng*M_PI/180.0;
	int tempStindex;
	for(int i=0; i<NumOfSect;i++)
	{
		if(asymSects[i].thetaSect >= tempStartAngle)
		{	
			
			tempStindex = i;
			break;
		}
	}

	
	double VoluteR = Sqrt((asymSects[tempStindex].x[NPnt/2]) *( asymSects[tempStindex].x[NPnt/2]) +(asymSects[tempStindex].y[NPnt/2]) *( asymSects[tempStindex].y[NPnt/2]) );
	double AnglDelta = TongueR*AspectRatio/VoluteR;
	double AngTranslate = -TranslateFactor/VoluteR;
	//find out the stating section correspond to main volute
	
	double StartAngleSelected = (TongueStAng*M_PI/180.0 +AngTranslate)> 0 ? TongueStAng*M_PI/180.0 +AngTranslate +AnglDelta  : AnglDelta;
	int TongueStartindex = 0;
	for(int i=0; i<NumOfSect;i++)
	{
		if(asymSects[i].thetaSect >= StartAngleSelected)
		{	
			
			TongueStartindex = i;
//			i = NumOfSect;
			/*gp_Pnt tempPt = gp_Pnt(asymSects[i-1].x[0],asymSects[i-1].y[0],asymSects[i-1].z[0]);
			BRepTools::Write(BRepBuilderAPI_MakeVertex(tempPt), "tempPt.brep");
			BRepTools::Write(seg0Wires[i], "seg0WiresAt8Deg.brep");
			double AngleDeg = asymSects[i].thetaSect *180.0/M_PI;*/
			break;
		}
	}

	int TongueMidStartindex = 0;
	for(int i=0; i<NumOfSect;i++)
	{
		if(asymSects[i].thetaSect >= TongueStAng*M_PI/180.0 +AnglDelta)
		{	
			
			TongueMidStartindex = i;
//			i = NumOfSect;
			/*gp_Pnt tempPt = gp_Pnt(asymSects[i-1].x[0],asymSects[i-1].y[0],asymSects[i-1].z[0]);
			BRepTools::Write(BRepBuilderAPI_MakeVertex(tempPt), "tempPt.brep");
			BRepTools::Write(seg0Wires[i], "seg0WiresAt8Deg.brep");
			double AngleDeg = asymSects[i].thetaSect *180.0/M_PI;*/
			break;
		}
	}

	TopoDS_Wire TongueStWire = seg0Wires[TongueStartindex];
	TopoDS_Wire VoluteEndWire = seg0Wires[seg0Wires.size()-1];
	//BRepTools::Write(TongueStWire, "TongueStWire.brep");

	TopoDS_Edge TongueStEdge;
	for(TopExp_Explorer PVexplr(TongueStWire,TopAbs_EDGE); PVexplr.More(); PVexplr.Next())
		{
			TongueStEdge = TopoDS::Edge(PVexplr.Current());
		}

	double newAnglei=(asymSects[TongueStartindex+1].thetaSect)*180.0/M_PI;
	double newAngleim1=(asymSects[TongueStartindex].thetaSect)*180.0/M_PI;


	BRepOffsetAPI_ThruSections BTS1(Standard_False,Standard_False,1.0e-06);
	
	int TongueStEndIndex = TongueStartindex +16;

	TongueTransEndIndex = TongueStEndIndex;

	for(int i=TongueStartindex; i<=TongueStEndIndex; i++)
	{
	if(i== TongueStartindex || i == TongueStEndIndex || i%3 ==0 )
	BTS1.AddWire(  seg0Wires[i]);
	/*if( i scrollSectWires.size()-5)
	BTS2.AddWire(  scrollSectWires[i]);*/
	}

	BTS1.Build();
	

	//BRepTools::Write(exitTranSectWires[0], "exitTranSectWires0.brep");
	//BRepTools::Write(exitTranSectWires[1], "exitTranSectWires1.brep");
	//BRepTools::Write(exitTranSectWires[2], "exitTranSectWires2.brep");
	//BRepTools::Write(exitTranSectWires[3], "exitTranSectWires3.brep");
	TopoDS_Shape VoluteStart = BTS1.Shape();
	//BRepTools::Write(BTS1.Shape(), "VoluteStart.brep");
	
	TopoDS_Wire VoluteEndWireClosed;

	{// Cut seg0Wires[last] into two at the top middle
		TopoDS_Edge LastSeg0WiresEdge ;
		for(TopExp_Explorer PVexplr(seg0Wires[seg0Wires.size()-1],TopAbs_EDGE); PVexplr.More(); PVexplr.Next())
		{
			LastSeg0WiresEdge = TopoDS::Edge(PVexplr.Current());
		}
		double u1, u2;
		Handle_Geom_Curve LastSeg0WiresEdgeCRV=BRep_Tool::Curve(LastSeg0WiresEdge,u1, u2);
		gp_Pnt PntSt = LastSeg0WiresEdgeCRV->Value(u1);
		gp_Pnt PntEnd = LastSeg0WiresEdgeCRV->Value(u2);
		TopoDS_Edge Eg1a= BRepBuilderAPI_MakeEdge(LastSeg0WiresEdgeCRV, u1, 0.5*(u1+ u2)).Edge();
		TopoDS_Edge Eg1b= BRepBuilderAPI_MakeEdge(LastSeg0WiresEdgeCRV, 0.5*(u1+ u2) , u2).Edge();
		TopoDS_Edge BotEdge = BRepBuilderAPI_MakeEdge(PntSt, PntEnd).Edge();
		BRepBuilderAPI_MakeWire WireNew;
		WireNew.Add(Eg1a);
		WireNew.Add(Eg1b);
		WireNew.Add(BotEdge);
		WireNew.Build();
		VoluteEndWireClosed = WireNew.Wire();


	}


	TopoDS_Shape transPipeTemp;
	  if(IsPipeShell) 
    {
        BRepOffsetAPI_MakePipeShell exitPipeShell(exitTransSpineWires[0]);
        exitPipeShell.Add(scrollSectWires[scrollSectWires.size()-1],Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
            exitPipeShell.Add(exitTranSectWires[1],Standard_False,Standard_False) ;
   
		try
		{
			exitPipeShell.Build();
		//	Standard_Boolean madeSolid = exitPipeShell.MakeSolid();
		//	exitPipeShell.MakeSolid();
		}
		catch(...)
		{
			BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
			BTS2.SetContinuity(GeomAbs_C1);
			BTS2.AddWire(VoluteEndWireClosed/*exitTranSectWires[0]*/);
			BTS2.AddWire(exitTranSectWires[1]);
			BTS2.Build();

			transPipeTemp =BTS2.Shape();
		//	BRepTools::Write( transPipeTemp, "VoluteEnd.brep");
			IsPipeShell = false;
		}
        
        transPipeTemp = exitPipeShell.Shape();		
	//	BRepTools::Write( transPipeTemp, "VoluteEnd.brep");
    }

	else

	{
		BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
		BTS2.SetContinuity(GeomAbs_C1);
		BTS2.AddWire( VoluteEndWireClosed /*scrollSectWires[scrollSectWires.size()-1] exitTranSectWires[0]*/);
		BTS2.AddWire(exitTranSectWires[1]);
		BTS2.Build();

		transPipeTemp =BTS2.Shape();
		//BRepTools::Write( transPipeTemp, "VoluteEnd.brep");
	}
	gp_Vec VoluteExitSt_Vec;
	TopoDS_Edge TranslatedVolEdge = TranslateVoluteEndToVotStart(transPipeTemp, TongueStWire, VoluteEndWire, IsPipeShell, VoluteExitSt_Vec);
	//BRepTools::Write(TranslatedVolEdge, "TranslatedVolEdge.brep");

	TopoDS_Wire TranslatedVolWireMidleFull ; // this will be evaluated only if angle is greater than 10 Degres
	if(TongueStAng > 10.0)
	{
		 TopoDS_Edge TranslatedVolEdgeMidle = TranslateVoluteEndToVotStart(transPipeTemp,  seg0Wires[TongueStartindex/2], VoluteEndWire, IsPipeShell, VoluteExitSt_Vec);
		
		TranslatedVolWireMidleFull = GetVolwireAtMidDistanceToStart(TranslatedVolEdgeMidle, seg0Wires[TongueStartindex/2]);
	}

	double U1, U2;
	Handle_Geom_Curve TransVolEdgeGeomCrv=BRep_Tool::Curve(TranslatedVolEdge,U1, U2);

	double V1, V2;
	Handle_Geom_Curve TongueStGeomCrv=BRep_Tool::Curve(TongueStEdge,V1, V2);

	TopoDS_Edge StartSegA = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, V1, 0.5*(V1+ V2));

//	BRepTools::Write(StartSegA, "StartSegA.brep");

	TopoDS_Edge StartSegB = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, 0.5*(V1+ V2), V2);

//	BRepTools::Write(StartSegB, "StartSegB.brep");

	BRepExtrema_ExtCC ExtCrvA(TranslatedVolEdge,StartSegA);
	
	int NumofExt=ExtCrvA.NbExt();
	double P1=ExtCrvA.ParameterOnE1(1);
	double Q1=ExtCrvA.ParameterOnE2(1);
	gp_Pnt SegABot=ExtCrvA.PointOnE2(1);

	BRepExtrema_ExtCC ExtCrvB(TranslatedVolEdge,StartSegB);
	
	int NumofExtB=ExtCrvB.NbExt();
	double P2=ExtCrvB.ParameterOnE1(1);
	double Q2=ExtCrvB.ParameterOnE2(1);
	gp_Pnt SegBBot=ExtCrvB.PointOnE2(1);

	if(true || Q1==0.0 || Q2 ==V2)
	{
		double Vdif = V2-V1;
		gp_Pnt V1pnt, V2pnt, tempPnt1, tempPnt2;
		gp_Vec tempVec1, tempVec2;
		TongueStGeomCrv->D1(V1,V1pnt,tempVec1);
		TongueStGeomCrv->D1(V2,V2pnt,tempVec2);
		for(int k=1; k< 50; k ++)
		{
			gp_Vec tempVec;
			gp_Pnt Temp1;
			TongueStGeomCrv->D1(V1+Vdif*0.01*k,Temp1,tempVec);
			if (tempVec.IsParallel(tempVec1,0.5))
			{
				k++;
			}
			else
			{
			Q1 = V1+Vdif*0.01*k;
			break;
			}
		}
			SegABot = TongueStGeomCrv->Value(Q1);
		for(int j=1; j<50; j ++)
		{
			gp_Vec tempVec;
			gp_Pnt Temp1;
			TongueStGeomCrv->D1(V2 -Vdif*0.01*j,Temp1,tempVec);
			if (tempVec.IsParallel(tempVec2,0.5))
			{
				j++;
			}
			else
			{
			Q2 = V2-Vdif*0.01*j;
			break;
			}
		}

		SegBBot = TongueStGeomCrv->Value(Q2);

	}

	TopoDS_Edge StartSegACut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, V1, Q1);

	//BRepTools::Write(StartSegACut, "StartSegACut.brep");

	TopoDS_Edge StartSegBCut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, Q2, V2);

	//BRepTools::Write(StartSegBCut, "StartSegBCut.brep");

	TopoDS_Edge StartSegTopCut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, Q1, Q2);
	//BRepTools::Write(StartSegTopCut, "StartSegTopCut.brep");

	gp_Pnt TongueBotMidPnt = TongueStGeomCrv->Value(0.5*(Q1+Q2));
	//BRepTools::Write(BRepBuilderAPI_MakeVertex(TongueBotMidPnt).Vertex(), "TongueBotMidPnt.brep");

	double P1new, P2new;
	double DeltaP = abs(P1-P2);

	if(P1<P2)
	{
		P1new = P1+DeltaP*0.1;
		P2new = P2-DeltaP*0.1;
	}
	else
	{
		P1new = P1-DeltaP*0.1;
		P2new = P2+DeltaP*0.1;
	}

	gp_Pnt Pnt_P1new = TransVolEdgeGeomCrv->Value(P1new);

	gp_Pnt Pnt_P2new = TransVolEdgeGeomCrv->Value(P2new);

	TopoDS_Edge ConnectSegA = BRepBuilderAPI_MakeEdge(SegABot, Pnt_P1new).Edge();

//	BRepTools::Write(ConnectSegA, "ConnectSegA.brep");

	TopoDS_Edge ConnectSegB = BRepBuilderAPI_MakeEdge(SegBBot, Pnt_P2new).Edge();

//	BRepTools::Write(ConnectSegB, "ConnectSegB.brep");

	TopoDS_Edge CutVolSeg = BRepBuilderAPI_MakeEdge(TransVolEdgeGeomCrv,P1new, P2new);

//	BRepTools::Write(CutVolSeg, "CutVolSeg.brep");

	BRepBuilderAPI_MakeWire TransFullWire( CutVolSeg);
	TransFullWire.Add(ConnectSegA);
	TransFullWire.Add(ConnectSegB);
	TransFullWire.Add(StartSegACut);
	TransFullWire.Add(StartSegBCut);
	TransFullWire.Build();

	//BRepTools::Write(VoluteEndWire, "VoluteEndWire.brep");
	//BRepTools::Write(TransFullWire.Wire(), "TransFullWire.brep");
	BRepOffsetAPI_ThruSections ThruTransShape(Standard_False,Standard_False,1.0e-06);
	ThruTransShape.AddWire(VoluteEndWire);
	if( TongueStAng>10.0)
	{
		if(!TranslatedVolWireMidleFull.IsNull())
			ThruTransShape.AddWire(TranslatedVolWireMidleFull);	
	}
	ThruTransShape.AddWire(TransFullWire.Wire());
	ThruTransShape.Build();

	//BRepTools::Write(ThruTransShape.Shape(), "ThruTransShape.brep");

	double FilletTopU = GetTopPointParameterOfTongueFillet(TongueR, SegABot, Pnt_P1new, SegBBot, Pnt_P2new);

//	double V1, V2;
	Handle_Geom_Curve SlantSegAGeomCrv=BRep_Tool::Curve(ConnectSegA,V1, V2);

	TopoDS_Edge SlantSegATop = BRepBuilderAPI_MakeEdge(SlantSegAGeomCrv,FilletTopU, V2);
	gp_Pnt SegATopSt = SlantSegAGeomCrv->Value(FilletTopU);
	gp_Pnt SegATopMid= SlantSegAGeomCrv->Value(FilletTopU/2.0);

	Handle_Geom_Curve SlantSegBGeomCrv=BRep_Tool::Curve(ConnectSegB,V1, V2);

	TopoDS_Edge SlantSegBTop = BRepBuilderAPI_MakeEdge(SlantSegBGeomCrv,FilletTopU, V2);
	gp_Pnt SegBTopSt = SlantSegBGeomCrv->Value(FilletTopU);
	gp_Pnt SegBTopMid = SlantSegBGeomCrv->Value(FilletTopU/2.0);

	TopoDS_Edge TongueTop = BRepBuilderAPI_MakeEdge(SegATopSt, SegBTopSt).Edge();
	gp_Pnt TongueTopMidPnt = gp_Pnt(0.5*(SegATopSt.X()+SegBTopSt.X()), 0.5*(SegATopSt.Y()+SegBTopSt.Y()), 0.5*(SegATopSt.Z()+SegBTopSt.Z()));
	{// to make more accurate TongueTopMidPnt
		TopoDS_Edge tempTongueTop = BRepBuilderAPI_MakeEdge(Pnt_P1new, Pnt_P2new).Edge();// temp parallel edge to Tongue Bottom
		/*double c1,c2;
		Handle_Geom_Curve tempTongueTopCRV=BRep_Tool::Curve(tempTongueTop,c1,c2);*/
		BRepExtrema_ExtPC ExtPC( BRepBuilderAPI_MakeVertex(TongueBotMidPnt).Vertex() ,tempTongueTop);
		//ExtPC.Perform();
			
			 gp_Pnt  tempTopMidPnt = ExtPC.Point(1);
			 gp_Vec VectoTop = gp_Vec(TongueBotMidPnt, tempTopMidPnt);
			 VectoTop.Normalize();
			 VectoTop.Multiply(2.0*TongueR);
			 TongueTopMidPnt = TongueBotMidPnt.Translated(VectoTop);

		}
	//BRepTools::Write(BRepBuilderAPI_MakeVertex(TongueTopMidPnt).Vertex(), "TongueTopMidPnt.brep");

	BRepBuilderAPI_MakeWire TongueTopVolWire( CutVolSeg);
	if(Tonguless)
		{
		TongueTopVolWire.Add(ConnectSegA);
		TongueTopVolWire.Add(ConnectSegB);
		TongueTopVolWire.Add(StartSegTopCut);
		}
	else
		{
		TongueTopVolWire.Add(SlantSegATop);
		TongueTopVolWire.Add(SlantSegBTop);
		TongueTopVolWire.Add(TongueTop);
		}
	TongueTopVolWire.Build();

	TopoDS_Wire TongueVolWire = TongueTopVolWire.Wire();

	//BRepTools::Write(TongueVolWire, "TongueVolWire.brep");

	BRepOffsetAPI_ThruSections ThruVolTopShape(Standard_False,Standard_False,1.0e-06);
	ThruVolTopShape.AddWire(TongueVolWire);
	ThruVolTopShape.AddWire(exitTranSectWires[1]);
	ThruVolTopShape.Build();
	TopoDS_Shape  VolutefromFilletEnd = ThruVolTopShape.Shape();

	//BRepTools::Write(exitTransSpineWires[0], "exitTransSpineWires[0].brep");
	if(IsPipeShell)
	{
		GProp_GProps SProps1;
        BRepBuilderAPI_MakeFace mkFace(TongueVolWire);
		if(mkFace.IsDone())
		{
	//	BRepTools::Write(mkFace.Face() ,  "mkFaceXX.brep");
		BRepGProp::SurfaceProperties(mkFace.Face(),SProps1);
		}
		else{
					BRepFill_Filling faceFiller(3,15,2,Standard_False,0.00001,0.0001,0.01,0.1,5,9);
				
				for(TopExp_Explorer Expl(TongueVolWire,TopAbs_EDGE); Expl.More(); Expl.Next())
				{
					TopoDS_Edge tempEdge = TopoDS::Edge(Expl.Current());
					faceFiller.Add(tempEdge,GeomAbs_C0);
				}
				faceFiller.Build();

				if(faceFiller.IsDone())
				{
				//	BRepTools::Write(faceFiller.Face() ,  "mkFaceXX2.brep");
					BRepGProp::SurfaceProperties(faceFiller.Face(),SProps1);
				}
		}
        
        gp_Pnt startSectionCentroid = SProps1.CentreOfMass();
	//	BRepTools::Write(BRepBuilderAPI_MakeVertex(startSectionCentroid).Vertex(), "startSectionCentroid.brep");
		TopoDS_Edge SpineEdge;
		{
				TopExp_Explorer Ex1;
				for (Ex1.Init(exitTransSpineWires[0],TopAbs_EDGE); Ex1.More(); Ex1.Next())
					{ 
						SpineEdge = TopoDS::Edge(Ex1.Current());
					}
		}
		BRepExtrema_ExtPC extPC1(BRepBuilderAPI_MakeVertex(startSectionCentroid).Vertex(),SpineEdge);
		
		double NearestU;
		if(extPC1.IsDone())
		{
		int NBpnts = extPC1.NbExt();
		NearestU =extPC1.Parameter(1);

		}
		double v1, v2;
		Handle_Geom_Curve SpineEdgeCrv=BRep_Tool::Curve(SpineEdge,v1, v2);
		TopoDS_Edge NewSpineEdge  = BRepBuilderAPI_MakeEdge(SpineEdgeCrv, NearestU, v2);

		BRepOffsetAPI_MakePipeShell exitPipeShell1(BRepBuilderAPI_MakeWire(NewSpineEdge).Wire());
        exitPipeShell1.Add(TongueVolWire,Standard_False,Standard_False) ;
        if(exit_enable_transition) //if transition in reqd add the second shape
            exitPipeShell1.Add(exitTranSectWires[1],Standard_False,Standard_False) ;
   
		try
		{
			exitPipeShell1.Build();

			VolutefromFilletEnd = exitPipeShell1.Shape();
		//	Standard_Boolean madeSolid = exitPipeShell.MakeSolid();
		//	exitPipeShell.MakeSolid();
		}
		catch(...)
		{
			VolutefromFilletEnd = ThruVolTopShape.Shape();
		}
	}
	//BRepTools::Write(VolutefromFilletEnd, "VolutefromFilletEnd.brep");

	ShapeFix_Shape FixShape1(VolutefromFilletEnd);
	FixShape1.Perform();
	//BRepTools::Write(FixShape1.Shape(), "FixShape1.brep");

	TopoDS_Face  TongueFilletShape;
	TopoDS_Shape NewThroughTransShape= GetTongueFilletAndNewTransShape(ThruTransShape.Shape(), TongueFilletShape, TongueR, AspectRatio, TranslateFactor, 
																					TongueStartindex, TongueMidStartindex,  Tonguless,
																					SegABot, SegATopMid, SegATopSt, SegBBot, SegBTopMid , SegBTopSt,
																					 TongueBotMidPnt,  TongueTopMidPnt, VoluteExitSt_Vec,
																					 SlantSegATop,  SlantSegBTop);

	//BRepTools::Write(NewThroughTransShape, "NewTrhoughTransShape.brep");
	//BRepTools::Write(TongueFilletShape, "TongueFilletShape.brep");
	TopoDS_Edge FilletTopEg, FilletBotEg;
	{// Get TongueFillet Edges
				 int j=0;
			TopExp_Explorer Ex;
			TopoDS_Edge Eg1, Eg2, Eg4;
		TopExp_Explorer Ex1;
		
			for (Ex1.Init(TongueFilletShape,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					FilletTopEg = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					FilletBotEg = TopoDS::Edge(Ex1.Current());
			}

	//		BRepTools::Write(Eg1, "Eg1.brep");
	//		BRepTools::Write(FilletTopEg, "FilletTopEg.brep");
	//		BRepTools::Write(FilletBotEg, "FilletBotEg.brep");
	
	}

	TopoDS_Shape NewVolutefromFilletEnd;
	{
			TopoDS_Face FaceA1,FaceA2, FaceA3, FaceB2, FaceB1;
    
			int j=0;
			TopExp_Explorer Ex;
			for (Ex.Init(VolutefromFilletEnd,TopAbs_FACE); Ex.More(); Ex.Next())
			{ 
				j++;
				if(j==1)
					{FaceA1 = TopoDS::Face(Ex.Current());
			//		BRepTools::Write(FaceA1, "FaceA1.brep");	
				}
				if(j==2)
				  {  FaceA2 = TopoDS::Face(Ex.Current());
			//		BRepTools::Write(FaceA2, "FaceA2.brep");
				}
				 if(j==3)
				  {  FaceA3 = TopoDS::Face(Ex.Current());
			//		BRepTools::Write(FaceA3, "FaceA3.brep");
				}
				 if(j==4)
				  {  FaceB1 = TopoDS::Face(Ex.Current());
			//		BRepTools::Write(FaceB1, "FaceB1.brep");
				}
				 if(j==5)
				  {  FaceB2 = TopoDS::Face(Ex.Current());
			//		BRepTools::Write(FaceB2, "FaceB2.brep");
				}
			}

			TopoDS_Edge Eg1, Eg2, Eg3, Eg4;
			TopExp_Explorer Ex1;
			j=0;
			for (Ex1.Init(FaceB1,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				if(j==3)
					Eg3 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

	//		BRepTools::Write(Eg1, "Eg1.brep");
	//		BRepTools::Write(Eg2, "Eg2.brep");
	//		BRepTools::Write(Eg3, "Eg3.brep");
	//		BRepTools::Write(Eg4, "Eg4.brep");


			BRepBuilderAPI_MakeWire BotWire;
			if(Tonguless)
			{BotWire.Add(FilletBotEg);}
			else
			{BotWire.Add(FilletTopEg);}
			{
			if(IsPipeShell)
			BotWire.Add(Eg1);
			else
			BotWire.Add(Eg2);
			}
			BotWire.Add(Eg3);
			BotWire.Add(Eg4);
			BotWire.Build();

			BRepBuilderAPI_MakeFace BotFace1(BRep_Tool::Surface(FaceB1), BotWire.Wire(), true);
			BotFace1.Build();

			ShapeFix_Face fixFace(BotFace1.Face());
				fixFace.Perform();
				TopoDS_Face  NewBotFace = fixFace.Face();

		//	BRepTools::Write(NewBotFace, "NewBotFace.brep");

			TopoDS_Shape  NewFaceA1;
			{
				TopoDS_Edge Eg1, Eg2, Eg3, Eg4;
			TopExp_Explorer Ex1;
			j=0;
			for (Ex1.Init(FaceA1,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				if(j==3)
					Eg3 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

	//		BRepTools::Write(Eg1, "Eg1.brep");
	//		BRepTools::Write(Eg2, "Eg2.brep");
	//		BRepTools::Write(Eg3, "Eg3.brep");
			if(IsPipeShell)
			{NewFaceA1 = FaceA1;}
			else
			{
					double x1, x2;
					Handle_Geom_Curve Eg1Crv=BRep_Tool::Curve(Eg1,x1, x2);
					TopoDS_Edge Eg1a= BRepBuilderAPI_MakeEdge(Eg1Crv, x1, 0.5*(x1+ x2)).Edge();
					TopoDS_Edge Eg1b= BRepBuilderAPI_MakeEdge(Eg1Crv, 0.5*(x1+ x2) , x2).Edge();
					BRepBuilderAPI_MakeWire Eg1Wire(Eg1a);
					Eg1Wire.Add(Eg1b);
					Eg1Wire.Build();

					BRepOffsetAPI_ThruSections ThruNewShape(Standard_False,Standard_False,1.0e-06);
					ThruNewShape.AddWire(Eg1Wire.Wire() );
					ThruNewShape.AddWire(BRepBuilderAPI_MakeWire(Eg3).Wire() );
					ThruNewShape.Build();
			
					 NewFaceA1 = ThruNewShape.Shape();

		//			BRepTools::Write(NewFaceA1, "NewFace1.brep");
			}
		}

		BRepBuilderAPI_Sewing sewing1(0.0001); 
		sewing1.Add(NewBotFace);
		sewing1.Add(FaceB2);
		sewing1.Add(FaceA3);
		sewing1.Add(NewFaceA1);
		sewing1.Add(FaceA2);
		sewing1.Perform();

		NewVolutefromFilletEnd = sewing1.SewedShape();
	//	BRepTools::Write(NewVolutefromFilletEnd, "NewVolutefromFilletEnd.brep");

	}

	TopoDS_Face NewVoluteStart;
	{// Get the new VoluteStart after adjesting for Tongue Fillet 
			TopoDS_Face FaceA1;
    
			int j=0;
			TopExp_Explorer Ex;
			for (Ex.Init(VoluteStart,TopAbs_FACE); Ex.More(); Ex.Next())
			{ 
				j++;
				if(j==1)
					{FaceA1 = TopoDS::Face(Ex.Current());
		//			BRepTools::Write(FaceA1, "FaceA1.brep");	
					}
				/*if(j==2)
				  {  FaceA2 = TopoDS::Face(Ex.Current());
					BRepTools::Write(FaceA2, "FaceA2.brep");}
				 if(j==3)
				  {  FaceTop1 = TopoDS::Face(Ex.Current());
					BRepTools::Write(FaceTop1, "FaceTop1.brep");}
				 if(j==4)
				  {  FaceB1 = TopoDS::Face(Ex.Current());
					BRepTools::Write(FaceB1, "FaceB1.brep");}
				 if(j==5)
				  {  FaceB2 = TopoDS::Face(Ex.Current());
					BRepTools::Write(FaceB2, "FaceB2.brep");}*/
			}

			TopoDS_Edge Eg1, Eg2, Eg3, Eg4;
			TopExp_Explorer Ex1;
			j=0;
			for (Ex1.Init(FaceA1,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				if(j==3)
					Eg3 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

		//	BRepTools::Write(Eg1, "Eg1.brep");
		//	BRepTools::Write(Eg2, "Eg2.brep");
	//		BRepTools::Write(Eg4, "Eg4.brep");

			BRepBuilderAPI_MakeWire BotWire;
			BotWire.Add(FilletBotEg);
			BotWire.Add(StartSegACut);
			BotWire.Add(StartSegBCut);
			BotWire.Add(Eg2);
			BotWire.Add(Eg3);
			BotWire.Add(Eg4);
			BotWire.Build();

			BRepBuilderAPI_MakeFace BotFace1(BRep_Tool::Surface(FaceA1), BotWire.Wire(), true);
			BotFace1.Build();

			ShapeFix_Face fixFace(BotFace1.Face());
				fixFace.Perform();
				NewVoluteStart = fixFace.Face();

	//		BRepTools::Write(NewVoluteStart, "NewVoluteStart.brep");

			

	}

	TopoDS_Shape TransitionBottom;
	{// make TransitionBottom
		
	//	BRepTools::Write( seg6Wires[0], "seg6Wires0.brep");
		BRepOffsetAPI_ThruSections BTSBot(Standard_False,Standard_False,1.0e-06);
	    BTSBot.SetContinuity(GeomAbs_C1);
	
			for(int i=0; i<=TongueStEndIndex; i++)
			{
			if(i== 0 || i == TongueStEndIndex || i%3 ==0 )
			BTSBot.AddWire(  seg6Wires[i]);
			/*if( i scrollSectWires.size()-5)
			BTS2.AddWire(  scrollSectWires[i]);*/
			}

			BTSBot.Build();

			TransitionBottom = BTSBot.Shape();
	}

	{
		BRepBuilderAPI_Sewing sewingWithTongue(1.0e-06); 
		sewingWithTongue.Add(NewThroughTransShape);
		sewingWithTongue.Add(NewVoluteStart);
		if(!Tonguless)
		{sewingWithTongue.Add(TongueFilletShape);}
		sewingWithTongue.Add(NewVolutefromFilletEnd);	
		sewingWithTongue.Add(TransitionBottom);		
		sewingWithTongue.Perform();
//		BRepTools::Write(sewingWithTongue.SewedShape(), "sewingWithTongue.brep");
		
		ShapeFix_Shape Fixshape(sewingWithTongue.SewedShape());
		Fixshape.Perform();
 	   TopoDS_Shape  sh1 = Fixshape.Shape();

//	   BRepTools::Write(sh1, "sh1.brep");

	   TongueTransitionPatch= sh1;
		
	   //add first part of input plane
	   //this is used to export stl file as segmented
	   m_inputPlane = TransitionBottom;
		
	//	BRepTools::Write(TongueTransitionPatch, "TongueTransitionPatch.brep");

	}

	

	/*
	-locate where to break the first scroll section so that the tongue radius tip point will be at the user specified angle
	-extract the outer curve of the scroll section at this point
	-create arcs for the tongue
	-create arcs or face for the sides of the tongue (starting from tongue arcs going to the scroll end section)
	-make the tongue surface using the above arcs
	-make the strating curve for the exit pipe on top of the tongue curve
	-do rest of the surface filling for the transition patch
	*/



	return true;
}

TopoDS_Shape OCCVolute::GetTongueFilletAndNewTransShape(TopoDS_Shape ThroughTransShape, TopoDS_Face & TongueFilletShape, double TongueR, double AspectRatio, double TranslateFactor, 
																					int TongueStartindex, int TongueMidStartindex, bool Tonguless,
																					gp_Pnt SegABot, gp_Pnt SegATopMid, gp_Pnt SegATopSt, gp_Pnt SegBBot, gp_Pnt SegBTopMid , gp_Pnt SegBTopSt,
																					gp_Pnt TongueBotMidPnt, gp_Pnt TongueTopMidPnt, gp_Vec VoluteExitSt_Vec,
																					TopoDS_Edge SlantSegATop, TopoDS_Edge SlantSegBTop)
{
	    // Fillet Curves on FaceA and FaceB
	double V1, V2;
	TopoDS_Face FaceA0,FaceTop, FaceA, FaceB, FaceB0;
    
    int j=0;
	TopExp_Explorer Ex;
    for (Ex.Init(ThroughTransShape,TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        j++;
        if(j==1)
            {FaceA0 = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(FaceA0, "FaceA0.brep");
		}
        if(j==2)
          {  FaceA = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(FaceA, "FaceA.brep");
	       }
		 if(j==3)
          {  FaceTop = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(FaceTop, "FaceTop.brep");
			}
		 if(j==4)
          {  FaceB = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(FaceB, "FaceB.brep");
			}
		 if(j==5)
          {  FaceB0 = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(FaceB0, "FaceB0.brep");
			}
    }

	Handle_Geom_BSplineCurve PVBSplineA, PVBSplineB;
	TopoDS_Face NewFaceA, NewFaceB;
	TopoDS_Edge TongueSideEgA, TongueSideEgB;
	gp_Vec TongueFilletDir; // Vector to Tongue Leading Direction , Vec size = TongueR

	//GetFaceAandFilletEdge(FaceA, NewFaceA, SegABot, SegATopSt, TongueSideEgA);
	{		  // Fillet Curves on FaceA and NewFaceA
		TopoDS_Edge Eg1, Eg2, Eg4;
		TopExp_Explorer Ex1;
		j=0;
			for (Ex1.Init(FaceA,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

	//		BRepTools::Write(Eg1, "Eg1.brep");
	//		BRepTools::Write(Eg2, "Eg2.brep");
	//		BRepTools::Write(Eg4, "Eg4.brep");

			Handle_Geom_Curve Eg4Crv=BRep_Tool::Curve(Eg4,V1, V2);
			gp_Pnt Eg4End ;
			gp_Vec Eg4EndVec;
			Eg4Crv->D1(V2, Eg4End, Eg4EndVec);
	//		BRepTools::Write(BRepBuilderAPI_MakeVertex(Eg4End).Vertex(), "Eg4End.brep");

			Handle_Geom_Surface pSurf2A = BRep_Tool::Surface(FaceA);
		 //  Handle_Geom_BSplineSurface pSurf2 = Handle(Geom_BSplineSurface)::DownCast(pSurf2a->Copy());
			//pSurf2A->Bounds(U1a,U2a,V1a,V2a);

			double U_A,V_A;
			GeomLib_Tool::Parameters(pSurf2A, SegATopMid, 0.001,U_A,V_A);

			Handle_Geom_Curve TempCurve=pSurf2A->UIso(U_A);
			
			TopoDS_Edge TempEdgeU1=BRepBuilderAPI_MakeEdge(TempCurve).Edge();

	//		BRepTools::Write(TempEdgeU1, "TempEdgeU1.brep");

	

			Handle_Geom_Curve TempEdgeU1Crv=BRep_Tool::Curve(TempEdgeU1,V1, V2);

			gp_Pnt TempMidPnt;
			gp_Vec PntA_Vec;
			TempEdgeU1Crv->D1(V2,TempMidPnt, PntA_Vec);

			PntA_Vec.Normalize();
			PntA_Vec.Multiply(-TongueR*AspectRatio);
			TongueFilletDir = PntA_Vec;

			gp_Pnt MidPntA =TempMidPnt.Translated(PntA_Vec);

	//		BRepTools::Write(BRepBuilderAPI_MakeVertex(MidPntA).Vertex(), "MidPntA.brep");

	//		GC_MakeEllipse 	TongueSide(SegATopSt, MidPntA,  TempMidPnt);
	//		Handle(Geom_Ellipse)  cGeom= TongueSide.Value();
	//		gp_Elips GpElips1 = cGeom->Elips();

	//		TongueSideEgA = BRepBuilderAPI_MakeEdge(GpElips1, SegATopSt ,SegABot).Edge();

	////		BRepTools::Write(TongueSideEgA, "TongueSideEgA.brep");

	//		Handle_Geom_Curve TongueSideEgACrv=BRep_Tool::Curve(TongueSideEgA,V1, V2);
	//		gp_Pnt EgACrvTopMid = TongueSideEgACrv->Value((3*V1+V2)/4);

	//		gp_Pnt EgACrvBotMid = TongueSideEgACrv->Value((V1+3*V2)/4);

			Handle_TColgp_HArray1OfPnt PointArrayA = new TColgp_HArray1OfPnt(1,3);
			PointArrayA->SetValue(1, SegABot);			
			PointArrayA->SetValue(2, MidPntA);
		//	PointArrayA->SetValue(3, EgACrvTopMid);
			PointArrayA->SetValue(3, SegATopSt);
		//	PointArrayA->SetValue(4, SegATopSt);

			GeomAPI_Interpolate PVinterpA(PointArrayA,Standard_False,1.0e-6);
			PVinterpA.Load(-Eg4EndVec,VoluteExitSt_Vec);
		
			PVinterpA.Perform();
			
			if(PVinterpA.IsDone())
				{
			   PVBSplineA = PVinterpA.Curve();

				}
			PVBSplineA->SetWeight(1, 2);
			PVBSplineA->SetWeight(3, 2);
	//			BRepTools::Write(BRepBuilderAPI_MakeEdge(PVBSplineA).Edge(), "PVBSplineA.brep");
		//	
		//	Handle_Geom_BezierCurve bez1 = new Geom_BezierCurve(*PointArrayA);
		//	bez1->SetWeight(2,5.0);
		//	
		//	BRepTools::Write(BRepBuilderAPI_MakeEdge(bez1).Edge(), "PVBSplineABez.brep");

			

			//GC_MakeArcOfEllipse FilletArc(GpElips1, SegABot, SegATopSt, false);

			
				
				BRepBuilderAPI_MakeWire FaceAwire(BRepBuilderAPI_MakeEdge(PVBSplineA).Edge());
				FaceAwire.Add(SlantSegATop);
				FaceAwire.Add(Eg2);
				FaceAwire.Add(Eg1);
				FaceAwire.Add(Eg4);
				FaceAwire.Build();
				bool isClosed = FaceAwire.Wire().Closed();
				BRepBuilderAPI_MakeFace  NewFaceA1(BRep_Tool::Surface(FaceA),FaceAwire.Wire());
				NewFaceA1.Build();

				ShapeFix_Face fixFace(NewFaceA1.Face());
				fixFace.Perform();
				NewFaceA = fixFace.Face();
				if(Tonguless)
				NewFaceA =FaceA;
	//			BRepTools::Write(NewFaceA, "NewFaceA.brep");
		}

	{	 // Fillet Curves on FaceB	 and NewFaceB
			TopoDS_Edge Eg1, Eg2, Eg4;
		TopExp_Explorer Ex1;
		j=0;
			for (Ex1.Init(FaceB,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

		/*	BRepTools::Write(Eg1, "Eg1.brep");
			BRepTools::Write(Eg2, "Eg2.brep");
			BRepTools::Write(Eg4, "Eg3.brep");*/

			Handle_Geom_Curve Eg2Crv=BRep_Tool::Curve(Eg2,V1, V2);
			gp_Pnt Eg2End ;
			gp_Vec Eg2EndVec;
			Eg2Crv->D1(V2, Eg2End, Eg2EndVec);
	//		BRepTools::Write(BRepBuilderAPI_MakeVertex(Eg2End).Vertex(), "Eg2End.brep");

			Handle_Geom_Surface pSurf2B = BRep_Tool::Surface(FaceB);
		 //  Handle_Geom_BSplineSurface pSurf2 = Handle(Geom_BSplineSurface)::DownCast(pSurf2a->Copy());
			//pSurf2A->Bounds(U1a,U2a,V1a,V2a);

			double U_A,V_A;
			GeomLib_Tool::Parameters(pSurf2B, SegBTopMid, 0.001,U_A,V_A);

			Handle_Geom_Curve TempCurve=pSurf2B->UIso(U_A);
			
			TopoDS_Edge TempEdgeU1=BRepBuilderAPI_MakeEdge(TempCurve).Edge();

	//		BRepTools::Write(TempEdgeU1, "TempEdgeU1.brep");

	

			Handle_Geom_Curve TempEdgeU1Crv=BRep_Tool::Curve(TempEdgeU1,V1, V2);

			gp_Pnt TempMidPnt;
			gp_Vec PntA_Vec;
			TempEdgeU1Crv->D1(V2,TempMidPnt, PntA_Vec);

			PntA_Vec.Normalize();
			PntA_Vec.Multiply(-TongueR*AspectRatio);

			gp_Pnt MidPntB =TempMidPnt.Translated(PntA_Vec);

	//		BRepTools::Write(BRepBuilderAPI_MakeVertex(MidPntB).Vertex(), "MidPntB.brep");

	//		GC_MakeEllipse 	TongueSide(SegBTopSt, MidPntB,  TempMidPnt);
	//		Handle(Geom_Ellipse)  cGeom= TongueSide.Value();
	//		gp_Elips GpElips1 = cGeom->Elips();

	//		TongueSideEgB = BRepBuilderAPI_MakeEdge(GpElips1, SegBTopSt ,SegBBot).Edge();

	////		BRepTools::Write(TongueSideEgB, "TongueSideEgB.brep");

	//		Handle_Geom_Curve TongueSideEgBCrv=BRep_Tool::Curve(TongueSideEgB,V1, V2);
	//		gp_Pnt EgBCrvTopMid = TongueSideEgBCrv->Value((3*V1+V2)/4);

	//		gp_Pnt EgBCrvBotMid = TongueSideEgBCrv->Value((V1+3*V2)/4);

			Handle_TColgp_HArray1OfPnt PointArrayB = new TColgp_HArray1OfPnt(1,3);
			PointArrayB->SetValue(1, SegBBot);
		//	PointArrayB->SetValue(2, EgBCrvBotMid);
			PointArrayB->SetValue(2, MidPntB);
		//	PointArrayB->SetValue(3, EgBCrvTopMid);
			PointArrayB->SetValue(3, SegBTopSt);

			GeomAPI_Interpolate PVinterpB(PointArrayB,Standard_False,1.0e-6);
			PVinterpB.Load(-Eg2EndVec, VoluteExitSt_Vec);
			PVinterpB.Perform();
			
			if(PVinterpB.IsDone())
				{
			   PVBSplineB = PVinterpB.Curve();
		//	   PVBSplineB->Reverse();
				}
	//			BRepTools::Write(BRepBuilderAPI_MakeEdge(PVBSplineB).Edge(), "PVBSplineB.brep");

				PVBSplineB->SetWeight(1, 2);
				PVBSplineB->SetWeight(3, 2);

			BRepBuilderAPI_MakeWire FaceAwire(BRepBuilderAPI_MakeEdge(PVBSplineB).Edge());
				FaceAwire.Add(SlantSegBTop);
				FaceAwire.Add(Eg2);
				FaceAwire.Add(Eg1);
				FaceAwire.Add(Eg4);
				FaceAwire.Build();
				bool isClosed = FaceAwire.Wire().Closed();
				BRepBuilderAPI_MakeFace  NewFaceB1(BRep_Tool::Surface(FaceB),FaceAwire.Wire());
				NewFaceB1.Build();

				ShapeFix_Face fixFace(NewFaceB1.Face());
				fixFace.Perform();
				NewFaceB = fixFace.Face();
				if(Tonguless)
				NewFaceB = FaceB;
	//			BRepTools::Write(NewFaceB, "NewFaceB.brep");

		}

		TopoDS_Edge FilletArcedge;
		{// middle ofTongue Fillet Edge
			gp_Pnt MidFilletCentre = gp_Pnt(0.5*(TongueBotMidPnt.X() +TongueTopMidPnt.X()), 0.5*(TongueBotMidPnt.Y() +TongueTopMidPnt.Y()) ,0.5*(TongueBotMidPnt.Z() +TongueTopMidPnt.Z()));
			
			gp_Pnt MidFilletFrontPnt = MidFilletCentre.Translated(TongueFilletDir);

			{
				Handle_TColgp_HArray1OfPnt PointArray = new TColgp_HArray1OfPnt(1,3);
			PointArray->SetValue(1, TongueBotMidPnt);			
			PointArray->SetValue(2, MidFilletFrontPnt);
		
			PointArray->SetValue(3, TongueTopMidPnt);
		

			GeomAPI_Interpolate PVinterp(PointArray,Standard_False,1.0e-6);
			PVinterp.Load(-VoluteExitSt_Vec,VoluteExitSt_Vec);
		
			PVinterp.Perform();
			Handle_Geom_BSplineCurve MidSpline;
			if(PVinterp.IsDone())
				{
			   MidSpline = PVinterp.Curve();
			//   MidSpline->Reverse();
				}

				MidSpline->SetWeight(1, 3);
				MidSpline->SetWeight(3, 3);

				FilletArcedge = BRepBuilderAPI_MakeEdge(MidSpline).Edge();
						BRepTools::Write(FilletArcedge, "FilletArcedge0.brep");
			}
			/*GC_MakeCircle FilletCirc(  TongueBotMidPnt, MidFilletFrontPnt , TongueTopMidPnt   );

			Handle(Geom_Circle)  cGeom= FilletCirc.Value();
			gp_Circ cCirc = cGeom->Circ();

			GC_MakeArcOfCircle FilletArc(cCirc, TongueBotMidPnt, TongueTopMidPnt, true);
		
			FilletArcedge = BRepBuilderAPI_MakeEdge(FilletArc.Value());*/
			if(!TranslateFactor==0)
			{
				
			TopExp_Explorer Ex1;
			TopoDS_Edge TempSegEdge;
			for (Ex1.Init(seg0Wires[TongueMidStartindex],TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				
					 TempSegEdge = TopoDS::Edge(Ex1.Current());
			}
			double e1,e2;
			Handle_Geom_Curve TempSegEdgeCrv=BRep_Tool::Curve(TempSegEdge, e1, e2);
			gp_Pnt SecTopPnt =   TempSegEdgeCrv->Value(0.5*(e1+e2));
			gp_Trsf T;
			gp_Vec FilletTranslateVec = gp_Vec(TongueBotMidPnt, SecTopPnt) ; // TranslateFactor*TongueFilletDir.Normalized();
			T.SetTranslation(FilletTranslateVec);
			FilletArcedge.Move(  T);
			}

			BRepTools::Write(FilletArcedge, "FilletArcedge.brep");
		}

		BRepOffsetAPI_ThruSections ThruTongueShape(Standard_False,Standard_False,1.0e-06);
		ThruTongueShape.SetContinuity(GeomAbs_C1);
		ThruTongueShape.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(PVBSplineA).Edge()).Wire());
		ThruTongueShape.AddWire(BRepBuilderAPI_MakeWire(FilletArcedge).Wire());
		ThruTongueShape.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(PVBSplineB).Edge()).Wire());
		ThruTongueShape.Build();

		TopoDS_Shape FilletShape = ThruTongueShape.Shape();
		TopExp_Explorer Ex1;
		int k=0;
	  for (Ex1.Init(FilletShape,TopAbs_FACE); Ex1.More(); Ex1.Next())
	  { 
			k++;
			if(k==1){
            TongueFilletShape = TopoDS::Face(Ex1.Current());
			BRepTools::Write(TongueFilletShape, "TongueFilletShape.brep");
			}
		}
		
		//if(TongueStartindex>2 )
		//{
		////ModifyFaceA0AndFaceB0( FaceA0, FaceB0, TongueStartindex);
		//}

		//sew all faces except ThruTongueShape.Shape()

		BRepBuilderAPI_Sewing sewing1(0.0001); 
		sewing1.Add(FaceB0);
		sewing1.Add(NewFaceB);
//		sewing1.Add(TongueFilletShape);
		sewing1.Add(FaceTop);
		sewing1.Add(NewFaceA);
		sewing1.Add(FaceA0);
		sewing1.Perform();

		return sewing1.SewedShape();

}

void OCCVolute::ModifyFaceA0AndFaceB0(TopoDS_Face & FaceA0, TopoDS_Face & FaceB0, int TongueStartIndex)
{
	BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
	BTS2.SetContinuity(GeomAbs_C1);
	BTS2.AddWire(scrollSectWires[0]);
	if(TongueStartIndex>2)
	BTS2.AddWire(scrollSectWires[int(TongueStartIndex/2)]);
	BTS2.AddWire(scrollSectWires[TongueStartIndex]);
	BTS2.Build();

	TopoDS_Face VoluteST;
	TopExp_Explorer Ex;
    for (Ex.Init(BTS2.Shape(),TopAbs_FACE); Ex.More(); Ex.Next())
    { 
       
            {VoluteST = TopoDS::Face(Ex.Current());
			//BRepTools::Write(TransFace2, "TransFace1.brep");	
			}
	}

	TopoDS_Edge  A0Eg, B0Eg;
	TopExp_Explorer Ex1;
		int j=0;
			for (Ex1.Init(VoluteST,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				
				if(j==2)
					A0Eg = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					B0Eg = TopoDS::Edge(Ex1.Current());
			}
						
	/*		BRepTools::Write(A0Eg, "A0Eg.brep");
			BRepTools::Write(B0Eg, "B0Eg.brep");*/

		{// Modify FaceA0
			TopoDS_Edge Eg1, Eg2,Eg3, Eg4;
		//TopExp_Explorer Ex1;
			j=0;
			for (Ex1.Init(FaceA0,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				if(j==3)
					Eg3 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

			/*BRepTools::Write(Eg1, "Eg1.brep");
			BRepTools::Write(Eg2, "Eg2.brep");
			BRepTools::Write(Eg4, "Eg4.brep");*/

			/*BRepFill_Filling faceFiller;
			faceFiller.Add(Eg1,GeomAbs_C1);
			faceFiller.Add(Eg2,GeomAbs_C1);
			faceFiller.Add(Eg3,GeomAbs_C1);
			faceFiller.Add(A0Eg,GeomAbs_C1);
			faceFiller.Build();*/

			BRepAlgo_NormalProjection FilletProj(FaceA0/*HubFace*/);
			FilletProj.Add(BRepBuilderAPI_MakeWire(A0Eg).Wire());
			FilletProj.SetDefaultParams();
			FilletProj.Build();

			TopoDS_Edge projectedEdge;
			if(FilletProj.IsDone())
		{
		//	BRepTools::Write(projProjection.Shape(), "h.brep");
			for(TopExp_Explorer explrE(FilletProj.Projection(),TopAbs_EDGE);explrE.More(); explrE.Next())
			{
				projectedEdge = TopoDS::Edge(explrE.Current());
			}
			//RightNewFilletEdge = TopoDS::Edge(projProjection.Shape());
		}

			BRepBuilderAPI_MakeWire FaceA0wire(Eg1);				
				FaceA0wire.Add(Eg2);
				FaceA0wire.Add(Eg3);	
				FaceA0wire.Add(projectedEdge);			
				FaceA0wire.Build();
				bool isClosed = FaceA0wire.Wire().Closed();
				BRepBuilderAPI_MakeFace  NewFaceA01(BRep_Tool::Surface(FaceA0),FaceA0wire.Wire());
				NewFaceA01.Build();

				ShapeFix_Face fixFace(NewFaceA01.Face());
				fixFace.Perform();
				FaceA0 = fixFace.Face();

		//		BRepTools::Write(fixFace.Face(), "NewFaceA0.brep");

		}

		{		// Modify FaceB0
				TopoDS_Edge Eg1, Eg2,Eg3, Eg4;
					j=0;
			for (Ex1.Init(FaceB0,TopAbs_EDGE); Ex1.More(); Ex1.Next())
			{ 
				j++;
				if(j==1)
					Eg1 = TopoDS::Edge(Ex1.Current());
				if(j==2)
					Eg2 = TopoDS::Edge(Ex1.Current());
				if(j==3)
					Eg3 = TopoDS::Edge(Ex1.Current());
				 if(j==4)
					Eg4 = TopoDS::Edge(Ex1.Current());
			}

		/*	BRepTools::Write(Eg1, "Eg1.brep");
			BRepTools::Write(Eg2, "Eg2.brep");
			BRepTools::Write(Eg4, "Eg4.brep");*/

			BRepAlgo_NormalProjection FilletProj(FaceB0/*HubFace*/);
			FilletProj.Add(BRepBuilderAPI_MakeWire(B0Eg).Wire());
			FilletProj.SetDefaultParams();
			FilletProj.Build();

			TopoDS_Edge projectedEdge;
			if(FilletProj.IsDone())
		{
		//	BRepTools::Write(projProjection.Shape(), "h.brep");
			for(TopExp_Explorer explrE(FilletProj.Projection(),TopAbs_EDGE);explrE.More(); explrE.Next())
			{
				projectedEdge = TopoDS::Edge(explrE.Current());
			}
			//RightNewFilletEdge = TopoDS::Edge(projProjection.Shape());
		}

			BRepBuilderAPI_MakeWire FaceB0wire(Eg1);				
				FaceB0wire.Add(projectedEdge);
				FaceB0wire.Add(Eg3);	
				FaceB0wire.Add(Eg4);			
				FaceB0wire.Build();
				bool isClosed = FaceB0wire.Wire().Closed();
				BRepBuilderAPI_MakeFace  NewFaceB01(BRep_Tool::Surface(FaceB0),FaceB0wire.Wire());
				NewFaceB01.Build();

				ShapeFix_Face fixFace(NewFaceB01.Face());
				fixFace.Perform();
				FaceB0 = fixFace.Face();

		//		BRepTools::Write(fixFace.Face(), "NewFaceB0.brep");

		}

}

TopoDS_Wire OCCVolute::GetVolwireAtMidDistanceToStart(TopoDS_Edge TranslatedVolEdgeMidle, TopoDS_Wire StartMidWire)
{
	
	TopoDS_Edge TongueStEdge;
	for(TopExp_Explorer PVexplr(StartMidWire,TopAbs_EDGE); PVexplr.More(); PVexplr.Next())
		{
			TongueStEdge = TopoDS::Edge(PVexplr.Current());
		}

	double U1, U2;
	Handle_Geom_Curve TransVolEdgeGeomCrv=BRep_Tool::Curve(TranslatedVolEdgeMidle, U1, U2);

	double V1, V2;
	Handle_Geom_Curve TongueStGeomCrv=BRep_Tool::Curve(TongueStEdge,V1, V2);

	TopoDS_Edge StartSegA = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, V1, 0.5*(V1+ V2));

//	BRepTools::Write(StartSegA, "StartSegA.brep");

	TopoDS_Edge StartSegB = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, 0.5*(V1+ V2), V2);

//	BRepTools::Write(StartSegB, "StartSegB.brep");

	BRepExtrema_ExtCC ExtCrvA(TranslatedVolEdgeMidle,StartSegA);
	
	int NumofExt=ExtCrvA.NbExt();
	double P1=ExtCrvA.ParameterOnE1(1);
	double Q1=ExtCrvA.ParameterOnE2(1);
	gp_Pnt SegABot=ExtCrvA.PointOnE2(1);

	BRepExtrema_ExtCC ExtCrvB(TranslatedVolEdgeMidle,StartSegB);
	
	int NumofExtB=ExtCrvB.NbExt();
	double P2=ExtCrvB.ParameterOnE1(1);
	double Q2=ExtCrvB.ParameterOnE2(1);
	gp_Pnt SegBBot=ExtCrvB.PointOnE2(1);

	
	{
		double Vdif = V2-V1;
		gp_Pnt V1pnt, V2pnt, tempPnt1, tempPnt2;
		gp_Vec tempVec1, tempVec2;
		TongueStGeomCrv->D1(V1,V1pnt,tempVec1);
		TongueStGeomCrv->D1(V2,V2pnt,tempVec2);
		for(int k=1; k< 50; k ++)
		{
			gp_Vec tempVec;
			gp_Pnt Temp1;
			TongueStGeomCrv->D1(V1+Vdif*0.01*k,Temp1,tempVec);
			if (tempVec.IsParallel(tempVec1,0.5))
			{
				k++;
			}
			else
			{
			if(k>1)
			Q1 = V1+Vdif*0.01*(k-1);
			else
			Q1 = V1+Vdif*0.01*(k);
			break;
			}
		}
			SegABot = TongueStGeomCrv->Value(Q1);
		for(int j=1; j<50; j ++)
		{
			gp_Vec tempVec;
			gp_Pnt Temp1;
			TongueStGeomCrv->D1(V2 -Vdif*0.01*j,Temp1,tempVec);
			if (tempVec.IsParallel(tempVec2,0.5))
			{
				j++;
			}
			else
			{
			if(j>1)
			Q2 = V2-Vdif*0.01*(j-1);
			else
			Q2 = V2-Vdif*0.01*(j);
			break;
			}
		}

		SegBBot = TongueStGeomCrv->Value(Q2);

	}

	TopoDS_Edge StartSegACut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, V1, Q1);

	//BRepTools::Write(StartSegACut, "StartSegACut.brep");

	TopoDS_Edge StartSegBCut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, Q2, V2);

	//BRepTools::Write(StartSegBCut, "StartSegBCut.brep");

	TopoDS_Edge StartSegTopCut = BRepBuilderAPI_MakeEdge(TongueStGeomCrv, Q1, Q2);
	//BRepTools::Write(StartSegTopCut, "StartSegTopCut.brep");

	gp_Pnt TongueBotMidPnt = TongueStGeomCrv->Value(0.5*(Q1+Q2));
	//BRepTools::Write(BRepBuilderAPI_MakeVertex(TongueBotMidPnt).Vertex(), "TongueBotMidPnt.brep");

	double P1new, P2new;
	double DeltaP = abs(P1-P2);

	if(P1<P2)
	{
		P1new = P1+DeltaP*0.1;
		P2new = P2-DeltaP*0.1;
	}
	else
	{
		P1new = P1-DeltaP*0.1;
		P2new = P2+DeltaP*0.1;
	}

	gp_Pnt Pnt_P1new = TransVolEdgeGeomCrv->Value(P1new);

	gp_Pnt Pnt_P2new = TransVolEdgeGeomCrv->Value(P2new);

	TopoDS_Edge ConnectSegA = BRepBuilderAPI_MakeEdge(SegABot, Pnt_P1new).Edge();

//	BRepTools::Write(ConnectSegA, "ConnectSegA.brep");

	TopoDS_Edge ConnectSegB = BRepBuilderAPI_MakeEdge(SegBBot, Pnt_P2new).Edge();

//	BRepTools::Write(ConnectSegB, "ConnectSegB.brep");

	TopoDS_Edge CutVolSeg = BRepBuilderAPI_MakeEdge(TransVolEdgeGeomCrv,P1new, P2new);

//	BRepTools::Write(CutVolSeg, "CutVolSeg.brep");
	TopoDS_Edge SegA  =mergeTwoEdges(StartSegACut,ConnectSegA);
	TopoDS_Edge SegB  =mergeTwoEdges(StartSegBCut,ConnectSegB);
	TopoDS_Edge VolSeg  =mergeTwoEdges(CutVolSeg,SegA);
	VolSeg = mergeTwoEdges(VolSeg,SegB);
	BRepBuilderAPI_MakeWire TransFullWire( VolSeg);
	/*TransFullWire.Add(ConnectSegA);
	TransFullWire.Add(ConnectSegB);
	TransFullWire.Add(StartSegACut);
	TransFullWire.Add(StartSegBCut);*/
	TransFullWire.Build();
	return TransFullWire.Wire();
}


TopoDS_Edge OCCVolute::mergeTwoEdges(TopoDS_Edge edge1, TopoDS_Edge edge2)
{
	TopoDS_Edge mergedEdge;
	Handle_Geom_Curve pGeoCurve,pGeoCurve2;
	Standard_Real first, last, first2, last2;

	pGeoCurve = BRep_Tool::Curve(edge1, first, last);
	Handle(Geom_TrimmedCurve) curve = new Geom_TrimmedCurve(pGeoCurve, first, last/*pGeoCurve->FirstParameter(), pGeoCurve->LastParameter()*/);
	GeomConvert_CompCurveToBSplineCurve final_spline(curve);
	
	pGeoCurve2 = BRep_Tool::Curve(edge2, first2, last2);
	Handle(Geom_TrimmedCurve) curve2 = new Geom_TrimmedCurve(pGeoCurve2,first2, last2 /*pGeoCurve2->FirstParameter(), pGeoCurve2->LastParameter()*/);

	final_spline.Add(curve2, MyPrecision::Precision::Approximation()/*Intersection()*//*Confusion()*/, Standard_False);
	BRepBuilderAPI_MakeEdge theEdgeBuilder_final(final_spline.BSplineCurve());
	//printShape(theEdgeBuilder_final,"theEdgeBuilder_final");

	return theEdgeBuilder_final;
}

//void OCCVolute::GetFaceAandFilletEdge(TopoDS_Face FaceA, TopoDS_Face &NewFaceA, gp_Pnt PntBottomA, gp_Pnt PntTopA, TopoDS_Edge FilletEdgeSideA, double TongueR)
//{		  // Fillet Curves on FaceA
//		double V1, V2;
//		TopoDS_Edge Eg1, Eg2, Eg4;
//		TopExp_Explorer Ex1;
//		j=0;
//			for (Ex1.Init(FaceA,TopAbs_EDGE); Ex1.More(); Ex1.Next())
//			{ 
//				j++;
//				if(j==1)
//					Eg1 = TopoDS::Edge(Ex1.Current());
//				if(j==2)
//					Eg2 = TopoDS::Edge(Ex1.Current());
//				 if(j==4)
//					Eg4 = TopoDS::Edge(Ex1.Current());
//			}
//
//		//	BRepTools::Write(Eg1, "Eg1.brep");
//		//	BRepTools::Write(Eg2, "Eg2.brep");
//	//		BRepTools::Write(Eg4, "Eg3.brep");
//
//			Handle_Geom_Curve Eg4Crv=BRep_Tool::Curve(Eg4,V1, V2);
//			gp_Pnt Eg4End ;
//			gp_Vec Eg4EndVec;
//			Eg4Crv->D1(V2, Eg4End, Eg4EndVec);
//	//		BRepTools::Write(BRepBuilderAPI_MakeVertex(Eg4End).Vertex(), "Eg4End.brep");
//
//			Handle_Geom_Surface pSurf2A = BRep_Tool::Surface(FaceA);
//		 //  Handle_Geom_BSplineSurface pSurf2 = Handle(Geom_BSplineSurface)::DownCast(pSurf2a->Copy());
//			//pSurf2A->Bounds(U1a,U2a,V1a,V2a);
//
//			double U_A,V_A;
//			GeomLib_Tool::Parameters(pSurf2A, SegATopMid, 0.001,U_A,V_A);
//
//			Handle_Geom_Curve TempCurve=pSurf2A->UIso(U_A);
//			
//			TopoDS_Edge TempEdgeU1=BRepBuilderAPI_MakeEdge(TempCurve).Edge();
//
//			BRepTools::Write(TempEdgeU1, "TempEdgeU1.brep");
//
//	
//
//			Handle_Geom_Curve TempEdgeU1Crv=BRep_Tool::Curve(TempEdgeU1,V1, V2);
//
//			gp_Pnt TempMidPnt;
//			gp_Vec PntA_Vec;
//			TempEdgeU1Crv->D1(V2,TempMidPnt, PntA_Vec);
//
//			PntA_Vec.Normalize();
//			PntA_Vec.Multiply(-TongueR);
//
//			gp_Pnt MidPntA =TempMidPnt.Translated(PntA_Vec);
//
//			BRepTools::Write(BRepBuilderAPI_MakeVertex(MidPntA).Vertex(), "MidPntA.brep");
//
//			GC_MakeEllipse 	TongueSide(PntTopA, MidPntA,  TempMidPnt);
//			Handle(Geom_Ellipse)  cGeom= TongueSide.Value();
//			gp_Elips GpElips1 = cGeom->Elips();
//
//			FilletEdgeSideA = BRepBuilderAPI_MakeEdge(GpElips1, PntTopA ,PntBottomA).Edge();
//
//			BRepTools::Write(FilletEdgeSideA, "FilletEdgeSideA.brep");
//
//			Handle_Geom_Curve TongueSideEgACrv=BRep_Tool::Curve(FilletEdgeSideA,V1, V2);
//			gp_Pnt EgACrvTopMid = TongueSideEgACrv->Value((3*V1+V2)/4);
//
//			gp_Pnt EgACrvBotMid = TongueSideEgACrv->Value((V1+3*V2)/4);
//
//			Handle_TColgp_HArray1OfPnt PointArrayA = new TColgp_HArray1OfPnt(1,4);
//			PointArrayA->SetValue(1, PntBottomA);			
//			PointArrayA->SetValue(2, MidPntA);
//			PointArrayA->SetValue(3, EgACrvTopMid);
//			PointArrayA->SetValue(4, PntTopA);
//		//	PointArrayA->SetValue(4, SegATopSt);
//
//			GeomAPI_Interpolate PVinterpA(PointArrayA,Standard_False,1.0e-6);
//		//	PVinterpA.Load(-Eg4EndVec,Eg4EndVec);
//		
//			PVinterpA.Perform();
//			
//			if(PVinterpA.IsDone())
//				{
//			   PVBSplineA = PVinterpA.Curve();
//
//				}
//				BRepTools::Write(BRepBuilderAPI_MakeEdge(PVBSplineA).Edge(), "PVBSplineA.brep");
//		//	
//		//	Handle_Geom_BezierCurve bez1 = new Geom_BezierCurve(*PointArrayA);
//		//	bez1->SetWeight(2,5.0);
//		//	
//		//	BRepTools::Write(BRepBuilderAPI_MakeEdge(bez1).Edge(), "PVBSplineABez.brep");
//
//			
//
//			//GC_MakeArcOfEllipse FilletArc(GpElips1, SegABot, SegATopSt, false);
//
//			
//				
//				BRepBuilderAPI_MakeWire FaceAwire(BRepBuilderAPI_MakeEdge(PVBSplineA).Edge());
//				FaceAwire.Add(SlantSegATop);
//				FaceAwire.Add(Eg2);
//				FaceAwire.Add(Eg1);
//				FaceAwire.Add(Eg4);
//				FaceAwire.Build();
//				bool isClosed = FaceAwire.Wire().Closed();
//				BRepBuilderAPI_MakeFace  NewFaceA1(BRep_Tool::Surface(FaceA),FaceAwire.Wire());
//				NewFaceA1.Build();
//
//				ShapeFix_Face fixFace(NewFaceA1.Face());
//				fixFace.Perform();
//				NewFaceA = fixFace.Face();
//
//				BRepTools::Write(NewFaceA, "NewFaceA.brep");
//		
//}

TopoDS_Edge OCCVolute::TranslateVoluteEndToVotStart(TopoDS_Shape VoluteEndShape, TopoDS_Wire VoluteStartWire, TopoDS_Wire VoluteEndWire, bool IsPipeShell, gp_Vec & VoluteExitSt_Vec)
{
	
		TopoDS_Face TransFace1, TransFace2, TransFace3;
    
    int j=0;
	TopExp_Explorer Ex;
    for (Ex.Init(VoluteEndShape,TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        j++;
        if(j==1)
            {TransFace1 = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(TransFace1, "TransFace1.brep");	
			}
        if(j==2)
          {  TransFace2 = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(TransFace2, "TransFace2.brep");
		}
		 if(j==3)
          {  TransFace3 = TopoDS::Face(Ex.Current());
	//		BRepTools::Write(TransFace3, "TransFace3.brep");
		}
	 }
	
	//{// Testing New Method

	//	
	//	BRepBuilderAPI_MakeFace mkF1(scrollSectWires[4],true);
	//	BRepTools::Write(scrollSectWires[4], "scrollSectWires4.brep");
	//	mkF1.Build();
	//	if(mkF1.IsDone())
	//	{
	//	BRepFill_Filling faceFiller;
	//		//BRepFill_Filling faceFiller(2,50,2,Standard_False,0.00001,0.0001,0.01,0.1,5,9);
	//		for(TopExp_Explorer explr(scrollSectWires[4],TopAbs_EDGE);explr.More(); explr.Next())
	//		{
	//			faceFiller.Add(TopoDS::Edge(explr.Current()),GeomAbs_C0);
	//		}
	//		faceFiller.Build();
	//	BRepTools::Write(faceFiller.Face(), "faceFiller.brep");
	//	Handle_Geom_Surface baseSurf = BRep_Tool::Surface(faceFiller.Face());

	//	TopoDS_Shape E_Wire, Eshape;

	//	BRepAlgoAPI_Section aSec(baseSurf, VoluteEndShape);
	//	aSec.Build();
	//	aSec.HasAncestorFaceOn2(E_Wire, Eshape);

	//			

	//	BRepTools::Write(E_Wire, "E_Wire.brep");
	//	BRepTools::Write(Eshape, "Eshape.brep");
	//	}

	//}


		 TopoDS_Edge TransEg1, TransEg2, TransEg3;
	TopExp_Explorer Ex1;
	j=0;
    for (Ex1.Init(TransFace1,TopAbs_EDGE); Ex1.More(); Ex1.Next())
    { 
        j++;
        if(j==1)
            TransEg1 = TopoDS::Edge(Ex1.Current());
        if(j==2)
            TransEg2 = TopoDS::Edge(Ex1.Current());
		 if(j==3)
            TransEg3 = TopoDS::Edge(Ex1.Current());
    }

	//BRepTools::Write(TransEg1, "TransEg1.brep");
	//BRepTools::Write(TransEg2, "TransEg2.brep");
	//BRepTools::Write(TransEg3, "TransEg3.brep");
	 
	double Umin, Umax, Vmin, Vmax;
	Handle_Geom_Curve TransEg2GeomCrv/*=BRep_Tool::Curve(TransEg2,Umin, Umax)*/;
	gp_Pnt PntMin; gp_Vec VecMin;
	if(IsPipeShell){
	TransEg2GeomCrv=BRep_Tool::Curve(TransEg3,Umin, Umax);
	TransEg2GeomCrv->D1(Umin,PntMin,VecMin);
	VoluteExitSt_Vec = VecMin;
	}
	else{
	TransEg2GeomCrv=BRep_Tool::Curve(TransEg2,Umin, Umax);
	TransEg2GeomCrv->D1(Umin,PntMin,VecMin);
	VoluteExitSt_Vec = VecMin;
	}
//	BRepTools::Write(BRepBuilderAPI_MakeVertex(PntMin).Vertex(), "PntMin.brep");

	

	TopoDS_Edge TongueStEdge;
	for(TopExp_Explorer PVexplr(VoluteStartWire,TopAbs_EDGE); PVexplr.More(); PVexplr.Next())
		{
			TongueStEdge = TopoDS::Edge(PVexplr.Current());
		}

//	BRepTools::Write(TongueStEdge, "TongueStEdge.brep");

	TopoDS_Edge VoluteEndEdge;
	for(TopExp_Explorer PVexplr(VoluteEndWire,TopAbs_EDGE); PVexplr.More(); PVexplr.Next())
		{
			VoluteEndEdge = TopoDS::Edge(PVexplr.Current());
		}

//	BRepTools::Write(VoluteEndEdge, "VoluteEndEdge.brep");


	Handle_Geom_Curve TongueStGeomCrv=BRep_Tool::Curve(TongueStEdge,Vmin, Vmax);

	GeomAPI_ExtremaCurveCurve GECCFirst(TransEg2GeomCrv, TongueStGeomCrv);

	int npt = GECCFirst.NbExtrema();
	double Ux, Vx;
	GECCFirst.LowerDistanceParameters(Ux, Vx);
	gp_Pnt Px=TransEg2GeomCrv->Value(Ux);

	double v0,v1;
	Handle_Geom_Curve VoluteEndGeomCrv=BRep_Tool::Curve(VoluteEndEdge,v0,v1);
	gp_Pnt P0 = VoluteEndGeomCrv->Value(v0);
	gp_Pnt P1 = VoluteEndGeomCrv->Value(v1);
	gp_Pnt P_trfm ;
	if(Px.Distance(P0) < Px.Distance(P1))
		{
			P_trfm= P0;
		}
	else
		{
			P_trfm= P1;
		}

	gp_Trsf  translateStart;
	
	translateStart.SetTranslation(PntMin,Px);
	
	BRepBuilderAPI_Transform translatedStart(VoluteEndEdge,translateStart); 
	TopoDS_Edge TranslatedV_Endedge = TopoDS::Edge(translatedStart.Shape());

	//	BRepTools::Write(TranslatedV_Endedge, "TranslatedV_EndedgeStart.brep");

	double minU = Ux;
	double maxU = 1.5*Ux ;
	double tempU = Ux*1.01 ; //0.5*(minU +maxU) ;
	GeomAPI_ExtremaCurveCurve GECC0(BRep_Tool::Curve(TranslatedV_Endedge,v0,v1), TongueStGeomCrv);
	double MinDist = GECC0.LowerDistance();
	double DistError = 100.0;

	
	
	int iMin=0;
	for(int i=0; i<30; i++)
	{
		double u0, u1;
		Handle_Geom_Curve TempGeomCrv=BRep_Tool::Curve(TranslatedV_Endedge,u0,u1);

		gp_Pnt Ptemp=TransEg2GeomCrv->Value(Ux + 0.01*Ux*i);
		
		gp_Trsf  translateT1;
	
		translateT1.SetTranslation(Px,Ptemp);
		TempGeomCrv->Transform(translateT1);
		//Handle_Geom_Curve TempGeomCrv=BRep_Tool::Curve(TranslatedV_Endedge,u0,u1);
		GeomAPI_ExtremaCurveCurve GECC(TempGeomCrv, TongueStGeomCrv);

		int npt = GECCFirst.NbExtrema();
		double Ut, Vt;
		GECC.LowerDistanceParameters(Ut, Vt);
		
		gp_Pnt Pt=TempGeomCrv->Value(Ut);
		gp_Pnt P_tongue=TongueStGeomCrv->Value(Vt);
		double CalcDist = GECC.LowerDistance();
		
		if(CalcDist  < DistError)
				{
					
					DistError =CalcDist;

					iMin = i;
				}
				
				
				
	}

	gp_Pnt Ptrans=TransEg2GeomCrv->Value(Ux + 0.01*Ux*iMin);
		
	/*BRepTools::Write(BRepBuilderAPI_MakeVertex(Ptrans).Vertex(), "Ptrans.brep");
	BRepTools::Write(BRepBuilderAPI_MakeVertex(Px).Vertex(), "Px.brep");
	BRepTools::Write(VoluteEndEdge, "VoluteEndEdge_BeforeTrans.brep");*/
	gp_Trsf  translateTend;
	
	translateTend.SetTranslation(	Px, Ptrans);
	BRepBuilderAPI_Transform translateEnd(TranslatedV_Endedge,translateTend); 

	TranslatedV_Endedge =TopoDS::Edge(translateEnd.Shape());
//	BRepTools::Write(TranslatedV_Endedge, "TranslatedV_Endedge.brep");

	/*if(IsPipeShell)
	{
		

		return TranslatedV_Endedge;
	}*/

	double v00,v01, Ve, Vt;
	Handle_Geom_Curve TransltedGeomCrv=BRep_Tool::Curve(TranslatedV_Endedge,v00,v01);
	GeomAPI_ExtremaCurveCurve GECC1(TransltedGeomCrv, TongueStGeomCrv);
	GECC1.LowerDistanceParameters(Ve,Vt);
	gp_Pnt gpVe = TransltedGeomCrv->Value(Ve);

	gp_Pnt gpVt = TongueStGeomCrv->Value(Vt);

	translateTend.SetTranslation(	gpVe, gpVt);

	BRepBuilderAPI_Transform translateEnd1(TranslatedV_Endedge,translateTend); 

	TranslatedV_Endedge =TopoDS::Edge(translateEnd1.Shape());

	//BRepTools::Write(TranslatedV_Endedge, "TranslatedV_Endedgex.brep");

	/*BRepOffsetAPI_ThruSections Trans1(Standard_False,Standard_False,1.0e-06);
	Trans1.AddWire(BRepBuilderAPI_MakeWire(VoluteEndEdge));
	Trans1.AddWire(BRepBuilderAPI_MakeWire(TranslatedV_Endedge));
	Trans1.Build();*/
//	BRepTools::Write(Trans1.Shape(), "Trans1.brep");
 //       gp_Pnt P1;
 //       gp_Vec V1;
 //       pexitSpineCurve->D1(0.0,P1,V1);//StartPoint();
 //       T.SetTranslation(P1,startSectionCentroid);
	//	pexitSpineCurve->Transform(T);

	
	return TranslatedV_Endedge;
}

double OCCVolute::GetTopPointParameterOfTongueFillet(double TongueRadius, gp_Pnt PntLeft, gp_Pnt PntLeftTop, gp_Pnt PntRight, gp_Pnt PntRightTop)
{
	double Bot = PntLeft.Distance(PntRight);

	double Top = PntLeftTop.Distance(PntRightTop);

	double SlantH= PntLeft.Distance(PntLeftTop);

	double R_BotLeft = sqrt(PntLeft.X()*PntLeft.X() +PntLeft.Y()*PntLeft.Y());
	double R_TopLeft = sqrt(PntLeftTop.X()*PntLeftTop.X() +PntLeftTop.Y()*PntLeftTop.Y());
	double A = R_TopLeft-R_BotLeft;
	double B = PntLeft.Z() - PntLeftTop.Z();
	
	
	//double H_vertical = Sqrt(SlantH*SlantH - (Top-Bot)*(Top-Bot)/4.0);

	double Slant_thickness = (SlantH/A)*TongueRadius*2.0;

	TopoDS_Edge SlantSegA = BRepBuilderAPI_MakeEdge(PntLeft, PntLeftTop).Edge();

	double V1, V2;
	Handle_Geom_Curve SlantSegAGeomCrv=BRep_Tool::Curve(SlantSegA,V1, V2);

	double deltaV = V2-V1;

	double Vx = V1 + (deltaV/SlantH)*Slant_thickness;

	gp_Pnt T_top_PntA = SlantSegAGeomCrv->Value(Vx);

	TopoDS_Edge SlantSegB = BRepBuilderAPI_MakeEdge(PntRight, PntRightTop).Edge();

	double U1, U2;
	Handle_Geom_Curve SlantSegBGeomCrv=BRep_Tool::Curve(SlantSegB,U1, U2);

	gp_Pnt T_top_PntB = SlantSegBGeomCrv->Value(Vx);

	TopoDS_Edge TongueTop = BRepBuilderAPI_MakeEdge(T_top_PntA, T_top_PntB).Edge();

//	BRepTools::Write(TongueTop, "TongueTop.brep");

	return Vx;
}

bool OCCVolute::makeVoluteWithoutPatch(int TongueTransEndIndex)
{
	BRepOffsetAPI_ThruSections BTS1(Standard_False,Standard_False,1.0e-06);
	BRepOffsetAPI_ThruSections BTS2(Standard_False,Standard_False,1.0e-06);
	 int Midindex = (scrollSectWires.size()-1 - TongueTransEndIndex)/2;
	for(int i=TongueTransEndIndex; i<=Midindex; i++)
	{
	if(i== TongueTransEndIndex || i == Midindex || i%3 ==0 )
	BTS1.AddWire(  scrollSectWires[i]);
	/*if( i scrollSectWires.size()-5)
	BTS2.AddWire(  scrollSectWires[i]);*/
	}

	for(int i=Midindex; i<scrollSectWires.size(); i++)
	{
	if(i== Midindex || i == scrollSectWires.size()-1 || i%3 ==0 )
	BTS2.AddWire(  scrollSectWires[i]);
	/*if( i scrollSectWires.size()-5)
	BTS2.AddWire(  scrollSectWires[i]);*/
	}

	BTS1.Build();
	BTS2.Build();

	BRepBuilderAPI_Sewing sewing1(0.0001); 
		sewing1.Add(BTS1.Shape());
		sewing1.Add(BTS2.Shape());
		sewing1.Perform();

	VoluteWithOutTransPatch   =sewing1.SewedShape() ;//New transition patch method

//	BRepTools::Write(VoluteWithOutTransPatch, "VoluteWithOutTransPatch.brep");

	//adding input plane which is used to export stl file as segmented faces
	BRep_Builder aBuilder; 
	TopoDS_Compound aCompond;
	aBuilder.MakeCompound(aCompond);
	if(!m_inputPlane.IsNull())
	{
		aBuilder.Add(aCompond,m_inputPlane);
	}

	TopExp_Explorer Ex;
    int i=0;
    for (TopExp_Explorer Ex(VoluteWithOutTransPatch,TopAbs_FACE); Ex.More(); Ex.Next())
    { 
        i++;
		if( i == 2)
		{
			aBuilder.Add(aCompond,Ex.Current());
		}
		if( i == 4)
		{
			aBuilder.Add(aCompond,Ex.Current());
			break;
		}
    }

	m_inputPlane = aCompond;
	return true;
}

bool OCCVolute::sewAllNewPatchMethod()
{
	
	BRepBuilderAPI_Sewing sewing1(0.0001); 
		sewing1.Add(VoluteWithOutTransPatch);
		sewing1.Add(TongueTransitionPatch);

		for(int i=0; i<exitPipeFaces.size() ; i++)
		{	
			sewing1.Add(exitPipeFaces[i]);
		}

		sewing1.Perform();

		totalVolute= sewing1.SewedShape();

		m_exitPlane = exitPipeFaces[exitPipeFaces.size() - 1];

		return true;

}

// =======================================================================
// Author	: Mahesh
// function : isSameFace
// purpose  : Compate two faces
// =======================================================================
bool OCCVolute::isSameFace(const TopoDS_Face& theFace1, const TopoDS_Face& theFace2, const double tol)
{
  TopExp_Explorer E(theFace1, TopAbs_EDGE);
  TopTools_ListOfShape LS1, LS2;
  for(; E.More(); E.Next()) LS1.Append(E.Current());

  E.Init(theFace2, TopAbs_EDGE);
  for(; E.More(); E.Next()) LS2.Append(E.Current());

  //Compare the number of edges in the faces
  if(LS1.Extent() != LS2.Extent()) return false;

  double aMin = RealFirst(), aMax = RealLast();
  double xminB1=aMax, yminB1=aMax, zminB1=aMax, xminB2=aMax, yminB2=aMax, zminB2=aMax;
  double xmaxB1=aMin, ymaxB1=aMin, zmaxB1=aMin, xmaxB2=aMin, ymaxB2=aMin, zmaxB2=aMin;

  for(E.Init(theFace1, TopAbs_VERTEX); E.More(); E.Next()) {
    gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(E.Current()));
    if(P.X() < xminB1) xminB1 = P.X();
    if(P.Y() < yminB1) yminB1 = P.Y();
    if(P.Z() < zminB1) zminB1 = P.Z();
    if(P.X() > xmaxB1) xmaxB1 = P.X();
    if(P.Y() > ymaxB1) ymaxB1 = P.Y();
    if(P.Z() > zmaxB1) zmaxB1 = P.Z();
  }

  for(E.Init(theFace2, TopAbs_VERTEX); E.More(); E.Next()) {
    gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(E.Current()));
    if(P.X() < xminB2) xminB2 = P.X();
    if(P.Y() < yminB2) yminB2 = P.Y();
    if(P.Z() < zminB2) zminB2 = P.Z();
    if(P.X() > xmaxB2) xmaxB2 = P.X();
    if(P.Y() > ymaxB2) ymaxB2 = P.Y();
    if(P.Z() > zmaxB2) zmaxB2 = P.Z();
  }

  //Compare the bounding boxes of both faces
  if(gp_Pnt(xminB1, yminB1, zminB1).Distance(gp_Pnt(xminB2, yminB2, zminB2)) > tol)
    return false;

  if(gp_Pnt(xmaxB1, ymaxB1, zmaxB1).Distance(gp_Pnt(xmaxB2, ymaxB2, zmaxB2)) > tol)
    return false;

  Handle(Geom_Surface) S1 = BRep_Tool::Surface(theFace1);
  Handle(Geom_Surface) S2 = BRep_Tool::Surface(theFace2);

  //Check if there a coincidence of two surfaces at least in two points
  double U11, U12, V11, V12, U21, U22, V21, V22;
  BRepTools::UVBounds(theFace1, U11, U12, V11, V12);
  BRepTools::UVBounds(theFace2, U21, U22, V21, V22);

  double rangeU = U12-U11;
  double rangeV = V12-V11;
  double U = U11 + rangeU/3.0;
  double V = V11 + rangeV/3.0;
  gp_Pnt P1 = S1->Value(U, V);
  U = U11+rangeU*2.0/3.0;
  V = V11+rangeV*2.0/3.0;
  gp_Pnt P2 = S1->Value(U, V);

  if (!GeomLib_Tool::Parameters(S2, P1, tol, U, V) || U < U21 || U > U22 || V < V21 || V > V22)
    return false;

  if (P1.Distance(S2->Value(U,V)) > tol) return false;

  if (!GeomLib_Tool::Parameters(S2, P2, tol, U, V) || U < U21 || U > U22 || V < V21 || V > V22)
    return false;

  if (P2.Distance(S2->Value(U, V)) > tol) return false;

  //Check that each edge of the Face1 has a counterpart in the Face2
  TopTools_MapOfOrientedShape aMap;
  TopTools_ListIteratorOfListOfShape LSI1(LS1);
  for(; LSI1.More(); LSI1.Next()) {
    TopoDS_Edge E = TopoDS::Edge(LSI1.Value());
    bool isFound = false;
    TopTools_ListIteratorOfListOfShape LSI2(LS2);
    for(; LSI2.More(); LSI2.Next()) {
      TopoDS_Shape aValue = LSI2.Value();
      if(aMap.Contains(aValue)) continue; //To avoid checking already found edge several times
	  TopoDS_Edge anEdge1 = TopoDS::Edge(aValue);
      if(isSameEdge(E, anEdge1,1e-5)) {
        aMap.Add(aValue);
        isFound = true;
        break;
      }
    }
    if(!isFound) return false;
  }

  return true;
}

// =======================================================================
// Author	: Mahesh
// function : isSameEdge
// purpose  : Compate two edges
// =======================================================================
bool OCCVolute::isSameEdge(const TopoDS_Edge& theEdge1, const TopoDS_Edge& theEdge2, const double tol )
{
  TopoDS_Vertex V11, V12, V21, V22;
  TopExp::Vertices(theEdge1, V11, V12);
  TopExp::Vertices(theEdge2, V21, V22);
  gp_Pnt P11 = BRep_Tool::Pnt(V11);
  gp_Pnt P12 = BRep_Tool::Pnt(V12);
  gp_Pnt P21 = BRep_Tool::Pnt(V21);
  gp_Pnt P22 = BRep_Tool::Pnt(V22);
  bool coincide = false;

  //Check that ends of edges coincide
  if(P11.Distance(P21) <= tol) {
    if(P12.Distance(P22) <= tol) coincide =  true;
  }
  else if(P11.Distance(P22) <= tol) {
    if(P12.Distance(P21) <= tol) coincide = true;
  }

  if(!coincide) return false;

  if (BRep_Tool::Degenerated(theEdge1))
    if (BRep_Tool::Degenerated(theEdge2)) return true;
    else return false;
  else
    if (BRep_Tool::Degenerated(theEdge2)) return false;

  double U11, U12, U21, U22;
  Handle(Geom_Curve) C1 = BRep_Tool::Curve(theEdge1, U11, U12);
  Handle(Geom_Curve) C2 = BRep_Tool::Curve(theEdge2, U21, U22);

  //Check that both edges has the same geometry
  double range = U12-U11;
  double U = U11+ range/3.0;
  gp_Pnt P1 = C1->Value(U);     //Compute a point on one third of the edge's length
  U = U11+range*2.0/3.0;
  gp_Pnt P2 = C1->Value(U);     //Compute a point on two thirds of the edge's length

  C2 = new Geom_TrimmedCurve(C2, U21, U22);

  if(!GeomLib_Tool::Parameter(C2, P1, tol, U) ||  U < U21 || U > U22)
    return false;

  if(P1.Distance(C2->Value(U)) > tol) return false;

  if(!GeomLib_Tool::Parameter(C2, P2, tol, U) || U < U21 || U > U22)
    return false;

  if(P2.Distance(C2->Value(U)) > tol) return false;

  return true;
}




