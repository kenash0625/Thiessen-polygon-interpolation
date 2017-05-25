#define  WXUSINGDLL
#include <wx/wx.h>
#include <wx/grid.h>
#include <list>
#include <vector>
#include "st2ws.h"

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

class MyLayer
{
public:
	wxPen origpen, identpen;
	wxBrush origbrush, identbrush;
	vector<int> identfid;
	string ofilepath;
	OGRFile ofile;
	void Draw(wxMemoryDC &, OGRFeature *pFea,MyCanvas &canv);

	void Draw(wxMemoryDC &,MyCanvas &canv);
	MyLayer(const string &str);
};
class MyCanvas :public wxScrolledWindow
{
	wxPoint pt1,pt2,ptident;
	wxBitmap bitmap;
	wxMemoryDC memdc;
	double zfac;
	
	int useablest;
	double d1stval;
	map<string, st2ws_ws> *wscdweight;
public:	
	enum Mode{Pan,ZoomIn,ZoomOut};
	Mode mode;
	MyCanvas(wxWindow *);
	void OnPaint(wxPaintEvent &event);
	void OnSize(wxSizeEvent& event);
	void OnMouseLDown(wxMouseEvent &event);
	void OnMouseLUp(wxMouseEvent &event);
	void OnMouseRDown(wxMouseEvent &event);
	void Zoom(double);

	double						m_World2DC, m_DC2World, m_Scale;
	OGREnvelope					m_rWorld, m_Extent,m_layerEnv;
	wxRect m_rDC;
	vector<MyLayer *> m_layers;
	
	double XWorld2DC(double x, bool bRound = true);
	double YWorld2DC(double y, bool bRound = true);
	void XYWorld2DC(wxPoint *dst, OGRRawPoint *src);
	double XDC2World(double x);
	double YDC2World(double y);
	OGREnvelope GetWorld(wxRect rClient);
	OGRRawPoint GetWorld(wxRect rClient, wxPoint ptClient);
	void SetExtent();
	void RandRun();
	wxDECLARE_EVENT_TABLE();
};
