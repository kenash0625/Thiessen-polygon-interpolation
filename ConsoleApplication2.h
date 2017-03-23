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
class MyCanvas;
class MyApp : public wxApp
{
public:
	MyFrame *frame;
	virtual bool OnInit();
};
class MyFrame : public wxFrame
{
public:
	wxGrid *sitegrid, *cellgrid;
	MyCanvas *canvas;
	MyFrame();
	void OnRandRun(wxCommandEvent& event);
	void OnHello(wxCommandEvent& event);
	void ShowRes(int nwsindex);
	void ZoomIn(wxCommandEvent &);
	void ZoomOut(wxCommandEvent &);	
	void Pan(wxCommandEvent &);
private:
	wxDECLARE_EVENT_TABLE();
};

class MyCanvas :public wxScrolledWindow
{
	wxPoint pt1,pt2,ptident;
	int wsident;

	wxBitmap bitmap;
	wxMemoryDC memdc;
	double zfac;
	void extractPolygon(OGRGeometry *pGeom,vector<int> &vParts,vector<OGRRawPoint> &vPts);
	void extractPolygon(geos::geom::Geometry *pGeom, vector<int> &vParts, vector<OGRRawPoint> &vPts);
public:	
	enum Mode{Pan,ZoomIn,ZoomOut};
	Mode mode;
	MyCanvas(wxWindow *);
	void OnPaint(wxPaintEvent &event);
	void OnSize(wxSizeEvent& event);
	void OnMouseLDown(wxMouseEvent &event);
	void OnMouseLUp(wxMouseEvent &event);
	void OnMouseRDown(wxMouseEvent &event);
	void Hello();
	void Zoom(double);

	double						m_World2DC, m_DC2World, m_Scale;
	OGREnvelope					m_rWorld, m_Extent,m_layerEnv;
	wxRect m_rDC;
	list<vector<int>> cells;
	list<vector<OGRRawPoint>> cellpts;
	list<vector<OGRRawPoint>> sitepolys;
	vector<OGRRawPoint> sitecoords;
	double XWorld2DC(double x, bool bRound = true);
	double YWorld2DC(double y, bool bRound = true);

	OGREnvelope GetWorld(wxRect rClient);
	OGRRawPoint GetWorld(wxRect rClient, wxPoint ptClient);
	void SetExtent();
	void XYWorld2DC(wxPoint *dst, OGRRawPoint *src);
	double XDC2World(double x);
	double YDC2World(double y);
	wxDECLARE_EVENT_TABLE();
};
