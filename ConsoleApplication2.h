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
	double xWorld2DC(double x, bool bRound = true)
	{
		x = (x - m_rWorld.MinX) * m_World2DC;

		return(bRound ? (int)(x < 0.0 ? x - 0.5 : x + 0.5) : x);
	}

	double yWorld2DC(double y, bool bRound = true)
	{
		y = (m_rWorld.MaxY - y) * m_World2DC - 1;

		return(bRound ? (int)(y < 0.0 ? y - 0.5 : y + 0.5) : y);
	}
	OGREnvelope Get_World(wxRect rClient)
	{
		double		d, dWorld, dClient;

		dClient = (double)rClient.GetHeight() / (double)rClient.GetWidth();
		dWorld = (m_Extent.MaxY - m_Extent.MinY) / (m_Extent.MaxX - m_Extent.MinX);

		if (dWorld > dClient)
		{
			d = (m_Extent.MaxX - m_Extent.MinX - (m_Extent.MaxY - m_Extent.MinY) / dClient) / 2.0;
			m_Extent.MinX += d;
			m_Extent.MaxX -= d;
		}
		else
		{
			d = (m_Extent.MaxY - m_Extent.MinY - (m_Extent.MaxX - m_Extent.MinX) * dClient) / 2.0;
			m_Extent.MinY += d;
			m_Extent.MaxY -= d;
		}
		return m_Extent;
	}
	OGRRawPoint Get_World(wxRect rClient, wxPoint ptClient)
	{
		double		d;
		OGREnvelope	rWorld(Get_World(rClient));

		ptClient.y = rClient.GetHeight() - ptClient.y;
		d = (rWorld.MaxX-rWorld.MinX) / (double)rClient.GetWidth();

		return(OGRRawPoint(
			rWorld.MinX + ptClient.x * d,
			rWorld.MinY + ptClient.y * d)
			);
	}

	//set m_Extent before calling.
	//rate: DC to World ,World to DC .
	//draw pts of shapefile
	void SetExtent()
	{
		m_rDC = GetClientRect();
		m_Extent = Get_World(m_rDC);
		m_rWorld = m_Extent;
		//-----------------------------------------------------
		// ensure cellsize in x-/y-direction are identical...
		double	dxdyDC = (double)m_rDC.GetWidth() / (double)m_rDC.GetHeight();
		double	dxdyWorld = (m_rWorld.MaxX - m_rWorld.MinX) / (m_rWorld.MaxY - m_rWorld.MinY);

		if (dxdyDC > dxdyWorld)
		{
			//m_rWorld.Inflate(0.5 * (m_rWorld.Get_YRange() * dxdyDC - m_rWorld.Get_XRange()), 0.0, false);
		}
		else if (dxdyDC < dxdyWorld)
		{
			//m_rWorld.Inflate(0.0, 0.5 * (m_rWorld.Get_XRange() / dxdyDC - m_rWorld.Get_YRange()), false);
		}

		//-----------------------------------------------------
		m_World2DC = (double)m_rDC.GetWidth() / (m_rWorld.MaxX - m_rWorld.MinX);
		m_DC2World = 1.0 / m_World2DC;

		Hello(0, 0, 0);
		Refresh();
	}
	void drdr(wxPoint *dst, OGRRawPoint *src) {
		dst->x = (int)xWorld2DC(src->x);
		dst->y = (int)yWorld2DC(src->y);
	}


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
