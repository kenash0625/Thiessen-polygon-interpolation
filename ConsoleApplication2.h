//#define  WXUSINGDLL
#include <wx/wx.h>
#include <wx/grid.h>
#include <list>
#include <vector>
namespace geos
{
	namespace geom
	{
		class Geometry;
	}
}
class MyFrame;
class MyParentFrame;
class MyIdentFrame;
class MyCanvas :public wxPanel
{
	MyIdentFrame *pIdent;
	wxPoint pt1,pt2,ptident;
	wxBitmap bitmap;
	wxMemoryDC memdc;
	double zfac;
	void extractPolygon(OGRGeometry *pGeom,vector<int> &vParts,vector<OGRRawPoint> &vPts);
	void extractPolygon(geos::geom::Geometry *pGeom, vector<int> &vParts, vector<OGRRawPoint> &vPts);
public:
	MyCanvas(MyFrame *);
	void OnPaint(wxPaintEvent &event);
	void OnSize(wxSizeEvent& event);
	void OnMouseLDown(wxMouseEvent &event);
	void OnMouseLUp(wxMouseEvent &event);
	void OnMouseWheel(wxMouseEvent &event);
	void OnMouseRDown(wxMouseEvent &event);
	void Hello();


	double						m_World2DC, m_DC2World, m_Scale;
	OGREnvelope					m_rWorld, m_Extent;
	wxRect m_rDC;
	list<vector<int>> cells;
	list<vector<OGRRawPoint>> cellpts;
	list<vector<OGRRawPoint>> sitepolys;
	vector<OGRRawPoint> sitecoords;
	double xWorld2DC(double x, bool bRound = true);
	double yWorld2DC(double y, bool bRound = true);

	OGREnvelope Get_World(wxRect rClient);
	OGRRawPoint Get_World(wxRect rClient, wxPoint ptClient);
	void SetExtent();
	void xyWorld2DC(wxPoint *dst, OGRRawPoint *src);
	wxDECLARE_EVENT_TABLE();
};
class MyApp : public wxApp
{
public:
	MyParentFrame *frame;
	virtual bool OnInit();
};
class MyFrame : public wxMDIChildFrame
{
	MyCanvas *canvas;
public:
	MyFrame(wxMDIParentFrame *parent,const wxString& title, const wxPoint& pos, const wxSize& size);
	void RandRun();
private:
	wxDECLARE_EVENT_TABLE();
};
class MyParentFrame :public wxMDIParentFrame
{
public:
	MyParentFrame();
private:
	void OnHello(wxCommandEvent& event);
	void OnRandRun(wxCommandEvent& event);
	void OnExit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	MyFrame *subframe;
	wxDECLARE_EVENT_TABLE();
};
class MyIdentFrame:public wxMDIChildFrame
{
public:
	MyIdentFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size);
	void ShowRes(int wsindex);
	wxGrid *sitegrid,*cellgrid;
};