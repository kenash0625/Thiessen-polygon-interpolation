//#define  WXUSINGDLL
#include <wx/wx.h>

class MyFrame;
class MyCanvas :public wxPanel
{
	wxPoint pt1,pt2;
	wxBitmap bitmap;
	wxMemoryDC memdc;
	double zfac;
public:
	MyCanvas(MyFrame *);
	void OnPaint(wxPaintEvent &event);
	void OnSize(wxSizeEvent& event);
	void OnMouseLDown(wxMouseEvent &event);
	void OnMouseLUp(wxMouseEvent &event);
	void OnMouseWheel(wxMouseEvent &event);

	void Hello(int,int,double);


	double						m_World2DC, m_DC2World, m_Scale;
	OGREnvelope					m_rWorld, m_Extent;
	wxRect m_rDC;
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
	virtual bool OnInit();
};
class MyFrame : public wxMDIChildFrame
{
	MyCanvas *canvas;
public:
	MyFrame(wxMDIParentFrame *parent,const wxString& title, const wxPoint& pos, const wxSize& size);
private:
	wxDECLARE_EVENT_TABLE();
};
class MyParentFrame :public wxMDIParentFrame
{
public:
	MyParentFrame();
private:
	void OnHello(wxCommandEvent& event);
	void OnExit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	wxDECLARE_EVENT_TABLE();
};
