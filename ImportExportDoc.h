// ImportExportDoc.h : interface of the CImportExportDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_IMPORTEXPORTDOC_H__88A2147C_3B23_11D2_8E1E_0800369C8A03__INCLUDED_)
#define AFX_IMPORTEXPORTDOC_H__88A2147C_3B23_11D2_8E1E_0800369C8A03__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <ColoredShapes.h>
#include <OCC_3dDoc.h>
#include "BRepPrimAPI_MakeBox.hxx"
#include "Geom_BezierCurve.hxx"
#include <BRepGProp.hxx>
#include "GProp_GProps.hxx"
#include <BRepOffsetAPI_ThruSections.hxx>
#include "gp_Trsf.hxx"
#include "VoluteDialog.h"


//#include "BRepFeat_MakeCylindricalHole.hxx"
//#include "D:\OCCT\opencascade-7.0.0\samples\mfc\standard\05_ImportExport\adm\win\vc11\FilletDialog.h"
class CImportExportDoc : public OCC_3dDoc
{
	DECLARE_DYNCREATE(CImportExportDoc)
public:
	CImportExportDoc();
	virtual ~CImportExportDoc();
	virtual void Serialize(CArchive& ar);

	void ActivateFrame(CRuntimeClass* pViewClass, int nCmdShow = SW_RESTORE  );
    virtual void Popup (const Standard_Integer  x       ,
		    			const Standard_Integer  y       ,
                        const Handle(V3d_View)& aView   ); 


	double getTrapezuimHeight(double,double,double,double);
	 double getSurfaceArea(TopoDS_Wire);
	double numberTesting();
	TopoDS_Wire createNewShapeWithRightArea(double,double,double,double,double,double);
	TopoDS_Wire curveConstructor(double,double,double,double,double);
	TopoDS_Wire createVolute();
	TopoDS_Edge createCurveEdge(double,double,double,double,double);
	TopoDS_Wire createTrapazium(double,double,double,double);
	TopoDS_Wire createNewTrapazium(double,double,double,double,double,double);
	TopoDS_Wire createNewShapeAccordingToR1height(double,double,double,double,double);

// Implementation
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Generated message map functions
protected:
	//{{AFX_MSG(CImportExportDoc)
	afx_msg void OnFileImportIges();
	afx_msg void OnFileExportIges();
	afx_msg void OnFileImportStep();
	afx_msg void OnFileExportStep();
	afx_msg void OnFileImportBrep();
//	afx_msg void OnWindowNew3d();
	afx_msg void OnFileExportVrml();
	afx_msg void OnFileExportStl();
	afx_msg void OnBox();
	afx_msg void OnCylinder();
	afx_msg void OnObjectRemove();
	afx_msg void OnObjectErase();
	afx_msg void OnObjectDisplayall();
	afx_msg void OnBox1();
	afx_msg void OnBottle();
	afx_msg void OnTest();
	afx_msg void OnBREPFile();
	afx_msg void OnIGESFile();
	afx_msg void OnSTEPFile();
	afx_msg void OnFillet();
	afx_msg void OnFilletDialog();
	afx_msg void OnCut();
	afx_msg void OnMakeBoxDrill();
	afx_msg void OnVolute();
	
	
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

//Attributes
protected:
	CColoredShapes* m_pcoloredshapeList;
	BRepPrimAPI_MakeBox* boxPointer; //= new BRepPrimAPI_MakeBox();
	CVoluteDialog* voluteDlg;
	
	

};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_IMPORTEXPORTDOC_H__88A2147C_3B23_11D2_8E1E_0800369C8A03__INCLUDED_)
