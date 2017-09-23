
class OCCVoluteTwoLobe
{
	//initiating program,create compund shape
	void OnBearingVolute(double width,double exhaustFlankHeight,double bearingFlankHeight,double bearingSideAngle,double exhaustSideAngle,double wholeVoluteArea,
		double tipRadius,double dividerWallHeight,double dividerAngle,double exhaustThickness,
		double bearingThickness,double transitionPartLength,double exitPipeLength,double exitDividerAngle,double voluteRadius,double exitDividerWallWidth,double exitPipeRadius,
		double toungAreaPercentage);

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
	void makeTwoLobeVolute();

};