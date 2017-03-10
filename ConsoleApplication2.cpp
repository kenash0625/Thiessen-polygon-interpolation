// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".

#include "OGRFile.h"
#include "ConsoleApplication2.h"
#include "st2ws.h"
#include <fstream>
struct st2ws g_st2ws;
int main(int, wchar_t*[])
{
	return WinMain(::GetModuleHandle(NULL), NULL, NULL, SW_SHOWNORMAL);
}

enum
{
	ID_Hello = 1
};

wxBEGIN_EVENT_TABLE(MyFrame, wxMDIChildFrame)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyCanvas, wxPanel)
EVT_PAINT(MyCanvas::OnPaint)
EVT_SIZE(MyCanvas::OnSize)
EVT_LEFT_DOWN(MyCanvas::OnMouseLDown)
EVT_LEFT_UP(MyCanvas::OnMouseLUp)
EVT_MOUSEWHEEL(MyCanvas::OnMouseWheel)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyParentFrame, wxMDIParentFrame)
EVT_MENU(ID_Hello, MyParentFrame::OnHello)
EVT_MENU(wxID_EXIT, MyParentFrame::OnExit)
EVT_MENU(wxID_ABOUT, MyParentFrame::OnAbout)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(MyApp);

bool MyApp::OnInit()
{
	
	MyParentFrame *frame = new MyParentFrame();
	frame->Show(true);
	return true; 
}
/*
现象：如果流域多边形自交 并且 被泰森多边形切割 那么不会计算权重
即 函数WSRAIN返回-10 有外环的可以intersection 自交的不可以
重现：使用 柘荣 的 20160101-预报模型（一）_1_1 的流域shp（fid=2）和测站shp
函数： SelfIntersectGeom 用于调整点集 使多边形不自交
*/
void SelfIntersectPt(vector<OGRRawPoint> &vpts)
{
	for (int p0 = 0, p1 = 1, p2 = 2, p3 = 3; p3 < vpts.size(); ++p0, ++p1, ++p2, ++p3)
	{
		OGRLineString ls0, ls1;
		ls0.addPoint(vpts.at(p0).x, vpts.at(p0).y);
		ls0.addPoint(vpts.at(p1).x, vpts.at(p1).y);
		ls1.addPoint(vpts.at(p2).x, vpts.at(p2).y);
		ls1.addPoint(vpts.at(p3).x, vpts.at(p3).y);
		if (ls0.Touches(&ls1))
		{
			int i(0);
			--i;
		}
		else if (ls0.Intersects(&ls1))//自交的情况 0123->0x21x3
		{
			OGRGeometry *pInterGeom = ls0.Intersection(&ls1);
			OGRPoint *pInterPt = (OGRPoint *)pInterGeom;

			OGRRawPoint pt1 = vpts.at(p1), pinter(pInterPt->getX(), pInterPt->getY());
			vpts.at(p1) = vpts.at(p2);
			vpts.at(p2) = pt1;
			vpts.insert(vpts.begin() + p1, pinter);
			vpts.insert(vpts.begin() + p3, pinter);

		}
	}

}
void SelfIntersectGeom(OGRGeometry *pGeom)
{
	if (pGeom->getGeometryType() == wkbPolygon)
	{
		OGRPolygon *pPoly = (OGRPolygon *)pGeom;
		vector<OGRRawPoint> vpts(pPoly->getExteriorRing()->getNumPoints());
		pPoly->getExteriorRing()->getPoints(vpts.data());
		SelfIntersectPt(vpts);
		pPoly->getExteriorRing()->setPoints(vpts.size(), vpts.data());
	}
	else if (pGeom->getGeometryType() == wkbMultiPolygon)
	{
		OGRMultiPolygon *pMp = (OGRMultiPolygon *)pGeom;
		for (int i = 0; i < pMp->getNumGeometries(); ++i)
		{
			OGRPolygon *pPoly = (OGRPolygon *)pMp->getGeometryRef(i);
			vector<OGRRawPoint> vpts(pPoly->getExteriorRing()->getNumPoints());
			pPoly->getExteriorRing()->getPoints(vpts.data());
			SelfIntersectPt(vpts);
			pPoly->getExteriorRing()->setPoints(vpts.size(), vpts.data());
		}
	}
}

MyCanvas::MyCanvas(MyFrame *parent):wxPanel(parent, -1, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),zfac(1)
{
	OGRFile oStfile("");
	OGRFile oWsfile("D:/MyFirstProject/fantastic-guacamole.git/trunk/a/WATA.shp");	
	OGRFeature *pFeature;
	int nRainsz,nStsize;
	oWsfile.m_pLayer->GetExtent(&m_Extent);
	oWsfile.m_pLayer->ResetReading();
	oStfile.m_pLayer->ResetReading();
	for (; pFeature = oWsfile.m_pLayer->GetNextFeature();OGRFeature::DestroyFeature(pFeature))
	{
		SelfIntersectGeom(pFeature->GetGeometryRef());
		g_st2ws.add_ws((geos::geom::Geometry*)pFeature->GetGeometryRef()->exportToGEOS(oWsfile.geosctx()), pFeature->GetFieldAsString("WSCD"));
	}
	for (nStsize=0; pFeature = oWsfile.m_pLayer->GetNextFeature();nStsize++, OGRFeature::DestroyFeature(pFeature))
	{
		OGRPoint *pt = (OGRPoint*)pFeature->GetGeometryRef();
		g_st2ws.add_st(pt->getX(), pt->getY(), pFeature->GetFieldAsString("STCD"));
	}
	g_st2ws.add_end();
	ifstream ifs;
	double drun;
	for (int i = 0; i < nRainsz;i++)
	{
		g_st2ws.begin_st_rain();
		for (int j = 0; j < nStsize; j++)
		{
			ifs >> drun;
			g_st2ws.st_rain(drun);
		}
		g_st2ws.calc(st2ws::DYNAMIC, st2ws::THIESSEN);
	}

}
/* 
DC的边框(m_rDC) 图层的边框(m_Extent)
根据两个边框的纵横比 修改图层的边框
   算出DC与图层的比例
   根据比例 画点
*/
void MyCanvas::OnPaint(wxPaintEvent & event)
{
	if (bitmap.IsOk())
	{
		wxPaintDC dc(this);
		dc.DrawBitmap(bitmap, 0, 0);
	}
}
void MyCanvas::OnSize(wxSizeEvent & event)
{
	SetExtent();
}
void MyCanvas::OnMouseLDown(wxMouseEvent & event)
{
	pt1=event.GetPosition();
}
void MyCanvas::OnMouseLUp(wxMouseEvent & event)
{
	pt2 = event.GetPosition();
	OGREnvelope rWorld(this->m_Extent);
	OGRRawPoint p1 = Get_World(GetClientRect(), pt1),
		p2 = Get_World(GetClientRect(), pt2);
	double dx = p1.x - p2.x, dy = p1.y - p2.y;
	rWorld.MinX += dx;
	rWorld.MaxX += dx;
	rWorld.MinY += dy;
	rWorld.MaxY += dy;
	m_Extent = rWorld;
	SetExtent();
}

void MyCanvas::OnMouseWheel(wxMouseEvent & event)
{
	static double dfac = 0.5;
	OGRRawPoint opt = Get_World(GetClientRect(), event.GetPosition());
	OGREnvelope rWorld(m_Extent);
	double dx = opt.x - (rWorld.MinX + rWorld.MaxX) / 2.0,
		dy = opt.y - (rWorld.MaxY + rWorld.MinY) / 2.0;
	rWorld.MinX += dx;
	rWorld.MinY += dy;
	rWorld.MaxX += dx;
	rWorld.MaxY += dy;
	double tfac = dfac;
	if (event.GetWheelRotation() < 0)
	{
	}
	else if (event.GetWheelRotation() > 0)
	{
		tfac = -dfac;
	}
	dx = (rWorld.MaxX - rWorld.MinX)*tfac / 2.0;
	dy = (rWorld.MaxY - rWorld.MinY)*tfac / 2.0;
	rWorld.MinX -= dx;
	rWorld.MaxX += dx;
	rWorld.MinY -= dy;
	rWorld.MaxY += dy;

	m_Extent = rWorld;
	SetExtent();
}

void MyCanvas::Hello(int shiftx,int shifty,double z)
{
	OGRFile ofile("D:/MyFirstProject/fantastic-guacamole.git/trunk/a/WATA.shp");
	wxSize sz = GetClientSize();
	wxPen pen(wxColour(0, 0, 0), 1, wxPENSTYLE_SOLID);
	bitmap.Create(sz);
	memdc.SelectObject(bitmap);
	memdc.SetPen(pen);
	memdc.SetBackground(wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_SOLID));
	memdc.Clear();

	for (OGRFeature *pfea; pfea = ofile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pfea))
	{
		OGRLineString *ls = (OGRLineString*)pfea->GetGeometryRef();
		vector<OGRRawPoint> ptls(ls->getNumPoints());
		vector<wxPoint> ptls2(ptls.size());
		ls->getPoints(ptls.data());
		for (int i=0;i<ptls.size();i++)
		{
			this->xyWorld2DC(&ptls2[i], &ptls[i]);
		}
		memdc.DrawLines(ptls2.size(), ptls2.data());
	}
	memdc.SelectObject(wxNullBitmap);
}

double MyCanvas::xWorld2DC(double x, bool bRound /*= true*/)
{
	x = (x - m_rWorld.MinX) * m_World2DC;

	return(bRound ? (int)(x < 0.0 ? x - 0.5 : x + 0.5) : x);
}

double MyCanvas::yWorld2DC(double y, bool bRound /*= true*/)
{
	y = (m_rWorld.MaxY - y) * m_World2DC - 1;

	return(bRound ? (int)(y < 0.0 ? y - 0.5 : y + 0.5) : y);
}
//保持了和窗口的纵横比一样
OGREnvelope MyCanvas::Get_World(wxRect rClient)
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
OGRRawPoint MyCanvas::Get_World(wxRect rClient, wxPoint ptClient)
{
	double		d;
	OGREnvelope	rWorld(Get_World(rClient));

	ptClient.y = rClient.GetHeight() - ptClient.y;
	d = (rWorld.MaxX - rWorld.MinX) / (double)rClient.GetWidth();

	return(OGRRawPoint(
		rWorld.MinX + ptClient.x * d,
		rWorld.MinY + ptClient.y * d)
		);
}

//set m_Extent before calling.
//rate: DC to World ,World to DC .
//draw pts of shapefile
void MyCanvas::SetExtent()
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
void MyCanvas::xyWorld2DC(wxPoint *dst, OGRRawPoint *src) {
	dst->x = (int)xWorld2DC(src->x);
	dst->y = (int)yWorld2DC(src->y);
}

MyFrame::MyFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxMDIChildFrame(parent, wxID_ANY, "title test")
{

	canvas = new MyCanvas(this);
	canvas->SetExtent();
}
MyParentFrame::MyParentFrame() : wxMDIParentFrame(NULL, wxID_ANY, "wxWidgets MDI Sample",
	wxDefaultPosition, wxSize(500, 400))
{
	wxMenu *menuFile = new wxMenu;
	menuFile->Append(ID_Hello, "&Hello...\tCtrl-H",
		"Help string shown in status bar for this menu item");
	menuFile->AppendSeparator();
	menuFile->Append(wxID_EXIT);
	wxMenu *menuHelp = new wxMenu;
	menuHelp->Append(wxID_ABOUT);
	wxMenuBar *menuBar = new wxMenuBar;
	menuBar->Append(menuFile, "&File");
	menuBar->Append(menuHelp, "&Help");
	SetMenuBar(menuBar);
	CreateStatusBar();
	SetStatusText("Welcome to wxWidgets!");
}
void MyParentFrame::OnExit(wxCommandEvent& event)
{
	Close(true);
}
void MyParentFrame::OnAbout(wxCommandEvent& event)
{
	wxMessageBox("This is a wxWidgets' Hello world sample",
		"About Hello World", wxOK | wxICON_INFORMATION);
}
void MyParentFrame::OnHello(wxCommandEvent& event)
{
	MyFrame *subframe = new MyFrame(this,"",wxPoint(),wxSize());
	subframe->Show(true);
}
