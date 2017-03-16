// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".

#include "OGRFile.h"
#include "ConsoleApplication2.h"
#include "st2ws.h"
#include <fstream>
#include <thread>
#include <geos/geom/geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LineString.h>
#include <wx/headerctrl.h>
struct st2ws g_st2ws;
int main(int, wchar_t*[])
{
	return WinMain(::GetModuleHandle(NULL), NULL, NULL, SW_SHOWNORMAL);
}

enum
{
	ID_Hello = 1,
	ID_RandRun,
};

wxBEGIN_EVENT_TABLE(MyFrame, wxMDIChildFrame)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyCanvas, wxPanel)
EVT_PAINT(MyCanvas::OnPaint)
EVT_SIZE(MyCanvas::OnSize)
EVT_LEFT_DOWN(MyCanvas::OnMouseLDown)
EVT_LEFT_UP(MyCanvas::OnMouseLUp)
EVT_MOUSEWHEEL(MyCanvas::OnMouseWheel)
EVT_RIGHT_DOWN(MyCanvas::OnMouseRDown)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyParentFrame, wxMDIParentFrame)
EVT_MENU(ID_Hello, MyParentFrame::OnHello)
EVT_MENU(ID_RandRun, MyParentFrame::OnRandRun)
EVT_MENU(wxID_EXIT, MyParentFrame::OnExit)
EVT_MENU(wxID_ABOUT, MyParentFrame::OnAbout)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(MyApp);

bool MyApp::OnInit()
{
	
	frame = new MyParentFrame();
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

void MyCanvas::extractPolygon(OGRGeometry *pGeom, vector<int> &vParts, vector<OGRRawPoint> &vPts)
{
	if (wkbPolygon == pGeom->getGeometryType())
	{
		OGRPolygon *pPoly = (OGRPolygon*)pGeom;
		OGRLinearRing *pRing = pPoly->getExteriorRing();
		vector<OGRRawPoint> pttmp(pRing->getNumPoints());
		pRing->getPoints(pttmp.data());
		vParts.push_back(pttmp.size());
		vPts.insert(vPts.end(), pttmp.begin(), pttmp.end());
		for (int j = 0; j < pPoly->getNumInteriorRings(); j++)
		{
			OGRLinearRing *pRing = pPoly->getInteriorRing(j);
			vector<OGRRawPoint> pttmp(pRing->getNumPoints());
			pRing->getPoints(pttmp.data());
			vParts.push_back(pttmp.size());
			vPts.insert(vPts.end(), pttmp.begin(), pttmp.end());
		}
	}
	else if (wkbMultiPolygon == pGeom->getGeometryType())
	{
		OGRMultiPolygon *pmPoly = (OGRMultiPolygon*)pGeom;
		for (int i = 0; i < pmPoly->getNumGeometries(); i++)
		{
			OGRPolygon *pPoly = (OGRPolygon*)pmPoly->getGeometryRef(i);
			OGRLinearRing *pRing = pPoly->getExteriorRing();
			vector<OGRRawPoint> pttmp(pRing->getNumPoints());
			pRing->getPoints(pttmp.data());
			vParts.push_back(pttmp.size());
			vPts.insert(vPts.end(), pttmp.begin(), pttmp.end());
			for (int j = 0; j < pPoly->getNumInteriorRings(); j++)
			{
				OGRLinearRing *pRing = pPoly->getInteriorRing(j);
				vector<OGRRawPoint> pttmp(pRing->getNumPoints());
				pRing->getPoints(pttmp.data());
				vParts.push_back(pttmp.size());
				vPts.insert(vPts.end(), pttmp.begin(), pttmp.end());
			}
		}
	}
}

void MyCanvas::extractPolygon(geos::geom::Geometry *pGeom, vector<int> &polys, vector<OGRRawPoint> &polypts)
{
	geos::geom::GeometryTypeId geomType = pGeom->getGeometryTypeId();
	if (geomType==geos::geom::GeometryTypeId::GEOS_MULTIPOLYGON)
	{
		geos::geom::MultiPolygon *pmp = dynamic_cast<geos::geom::MultiPolygon*>(pGeom);
		for (int z = 0; z < pmp->getNumGeometries();z++)
		{
			const geos::geom::Polygon *geospoly = dynamic_cast<const geos::geom::Polygon*>(pmp->getGeometryN(z));
			const geos::geom::LineString *geosls = geospoly->getExteriorRing();
			polys.push_back(geosls->getNumPoints());
			for (int i = 0; i < geosls->getNumPoints(); i++)
			{
				geos::geom::Point *geosppp = geosls->getPointN(i);
				OGRRawPoint opt(geosppp->getX(), geosppp->getY());
				polypts.push_back(opt);
			}
			for (int i = 0; i < geospoly->getNumInteriorRing(); i++)
			{
				const geos::geom::LineString *geosls = geospoly->getInteriorRingN(i);
				polys.push_back(geosls->getNumPoints());
				for (int j = 0; j < geosls->getNumPoints(); j++)
				{
					geos::geom::Point *geosppp = geosls->getPointN(i);
					OGRRawPoint opt(geosppp->getX(), geosppp->getY());
					polypts.push_back(opt);
				}
			}
		}
	}
	else if (geomType == geos::geom::GeometryTypeId::GEOS_POLYGON)
	{
		geos::geom::Polygon *geospoly = dynamic_cast<geos::geom::Polygon*>(pGeom);
		const geos::geom::LineString *geosls = geospoly->getExteriorRing();
		polys.push_back(geosls->getNumPoints());
		for (int i = 0; i < geosls->getNumPoints(); i++)
		{
			geos::geom::Point *geosppp = geosls->getPointN(i);
			OGRRawPoint opt(geosppp->getX(), geosppp->getY());
			polypts.push_back(opt);
		}
		for (int i = 0; i < geospoly->getNumInteriorRing(); i++)
		{
			const geos::geom::LineString *geosls = geospoly->getInteriorRingN(i);
			polys.push_back(geosls->getNumPoints());
			for (int j = 0; j < geosls->getNumPoints(); j++)
			{
				geos::geom::Point *geosppp = geosls->getPointN(i);
				OGRRawPoint opt(geosppp->getX(), geosppp->getY());
				polypts.push_back(opt);
			}
		}
	}
	
}

MyCanvas::MyCanvas(MyFrame *parent):wxPanel(parent, -1, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),zfac(1)
{
	pIdent = new MyIdentFrame(wxGetApp().frame, "", wxPoint(), wxSize());
	OGRFile oStfile("D:/shuju/窑邦-江西/窑邦-江西/Maps/JCZD.shp");
	OGRFile oWsfile("D:/shuju/窑邦-江西/窑邦-江西/Maps/WATA.shp");	
	OGRFeature *pFeature;
	int nRainsz=0,nStsize;
	oWsfile.m_pLayer->GetExtent(&m_Extent);
	oWsfile.m_pLayer->ResetReading();
	oStfile.m_pLayer->ResetReading();
	for (; pFeature = oWsfile.m_pLayer->GetNextFeature();OGRFeature::DestroyFeature(pFeature))
	{
		OGRGeometry *pGeom = pFeature->GetGeometryRef();
		SelfIntersectGeom(pGeom);
		g_st2ws.add_ws((geos::geom::Geometry*)pGeom->exportToGEOS(oWsfile.geosctx()), pFeature->GetFieldAsString("WSCD"));
		
		vector<int> ints;
		vector<OGRRawPoint> pts;
		extractPolygon(pGeom, ints, pts);
		cells.push_back(ints);
		cellpts.push_back(pts);
	}
	for (nStsize=0; pFeature = oStfile.m_pLayer->GetNextFeature();nStsize++, OGRFeature::DestroyFeature(pFeature))
	{
		OGRPoint *pt = (OGRPoint*)pFeature->GetGeometryRef();
		sitecoords.push_back(OGRRawPoint(pt->getX(), pt->getY()));
		g_st2ws.add_st(pt->getX(), pt->getY(), pFeature->GetFieldAsString("STCD"));
	}
	g_st2ws.add_end();
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

void MyCanvas::OnMouseRDown(wxMouseEvent &event)
{
	ptident= event.GetPosition(); 
	SetExtent();
}

void MyCanvas::Hello()
{
	wxSize sz = GetClientSize();
	wxPen pen(wxColour(0, 0, 0), 1, wxPENSTYLE_SOLID);
	wxBrush brush(wxColour(255, 255, 255), wxBRUSHSTYLE_TRANSPARENT);
	bitmap.Create(sz);
	memdc.SelectObject(bitmap);
	memdc.SetBrush(brush);
	memdc.SetPen(pen);
	//memdc.SetBackground(wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_SOLID));
	memdc.Clear();
	OGRRawPoint optident = Get_World(GetClientRect(), ptident);
	//all cells
	list<vector<int>>::iterator itpolys = cells.begin();
	list<vector<OGRRawPoint>>::iterator itpts = cellpts.begin();
	for (; itpolys != cells.end(); itpolys++, itpts++)
	{
		vector<wxPoint> wxpts(itpts->size());
		std::transform(itpts->begin(), itpts->end(), wxpts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; xyWorld2DC(&p, &opt); return p; });
		memdc.DrawPolyPolygon(itpolys->size(), itpolys->data(), wxpts.data(), 0, 0, wxWINDING_RULE);
	}
	const geos::geom::GeometryFactory *geomfac(geos::geom::GeometryFactory::getDefaultInstance());
	geos::geom::Coordinate coord(optident.x, optident.y);
	geos::geom::Point *geospt = geomfac->createPoint(coord);
	int nrun(0);
	for (geos::geom::Geometry*p : g_st2ws.vwsgeom)
	{
		if (p->contains(geospt))
		{
			wxPen pen(wxColour(0, 0, 255), 2, wxPENSTYLE_SOLID);
			memdc.SetPen(pen);
			vector<int> polys;
			vector<OGRRawPoint> polypts;
			extractPolygon(p, polys, polypts);
			vector<wxPoint> wxpolypts(polypts.size());
			std::transform(polypts.begin(), polypts.end(), wxpolypts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; xyWorld2DC(&p, &opt); return p; });
			memdc.DrawPolyPolygon(polys.size(), polys.data(), wxpolypts.data(), 0, 0, wxPolygonFillMode::wxWINDING_RULE);
			
			
			pIdent->ShowRes(nrun);
			pIdent->Show();
			break;
		}
		++nrun;
	}
	//all sites
	wxBrush p1brush(wxColour(0, 0, 0), wxBRUSHSTYLE_SOLID);
	memdc.SetBrush(p1brush);
	memdc.SetPen(pen);
	vector<OGRRawPoint>::iterator iter = sitecoords.begin();
	for (; iter != sitecoords.end(); iter++)
	{
		wxPoint p;
		xyWorld2DC(&p, &*iter);
		memdc.DrawCircle(p, 4);
	}
	//voronoi cells
	if (g_st2ws.useablest > 1) {
		memdc.SetBrush(brush);
		vector<wxPoint> wxpts(g_st2ws.weightrun->voropts.size());
		std::transform(g_st2ws.weightrun->voropts.begin(), g_st2ws.weightrun->voropts.end(), wxpts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; xyWorld2DC(&p, &opt); return p; });
		memdc.DrawPolyPolygon(g_st2ws.weightrun->voropoly.size(), g_st2ws.weightrun->voropoly.data(), wxpts.data(), 0, 0, wxWINDING_RULE);
	}
	//voronoi sites
	wxPen rpen(wxColour(255, 0, 0), 1, wxPENSTYLE_SOLID);
	wxBrush p2brush(wxColour(255, 0, 0), wxBRUSHSTYLE_SOLID);
	memdc.SetBrush(p2brush);
	memdc.SetPen(rpen);
	for (int i = 0; i < g_st2ws.useablest;i++)
	{
		wxPoint p;
		xyWorld2DC(&p, &sitecoords[g_st2ws.stidxrun[i]]);
		memdc.DrawCircle(p, 4);
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
// 图层的纵横比 窗口的纵横比
// 图层的纵横比较小
// 图层的y扩大一些

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

// 图层和窗口的宽度比 计算了一下
// 把图层的坐标按照比例转为窗口的
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

	Hello();
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

void MyFrame::RandRun()
{
	g_st2ws.rand_run();
	canvas->SetExtent();
}

MyParentFrame::MyParentFrame() : wxMDIParentFrame(NULL, wxID_ANY, "wxWidgets MDI Sample",
	wxDefaultPosition, wxSize(500, 400))
{
	wxMenu *menuFile = new wxMenu;
	menuFile->Append(ID_Hello, "&Hello...\tCtrl-H",
		"Help string shown in status bar for this menu item");
	menuFile->Append(ID_RandRun, "&RandRun...\tCtrl-R",
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
	subframe = new MyFrame(this,"",wxPoint(),wxSize());
	subframe->Show(true);
}

void MyParentFrame::OnRandRun(wxCommandEvent& event)
{
	subframe->RandRun();
}
//显示与某个流域、相关的泰森多边形、多边形雨量、流域雨量
MyIdentFrame::MyIdentFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size ) : wxMDIChildFrame(parent, wxID_ANY, "title2 test"),sitegrid(nullptr),cellgrid(nullptr)
{	
}

void MyIdentFrame::ShowRes(int nwsindex)
{
	delete sitegrid;
	delete cellgrid;
	int r(1);
	sitegrid = new wxGrid(this, -1, wxPoint(0, 0), wxSize(400, 300));
	sitegrid->CreateGrid(1 + g_st2ws.weightrun->weights[nwsindex].size(), 3);
	sitegrid->SetCellValue(0, 0, "SITE");
	sitegrid->SetCellValue(0, 1, "WEIGHT");
	sitegrid->SetCellValue(0, 2, "VALUE");
	for (auto &a : g_st2ws.weightrun->weights[nwsindex])
	{
		sitegrid->SetCellValue(r, 0, g_st2ws.vstcdrun[a.first]);
		sitegrid->SetCellValue(r, 1, std::to_string(a.second / g_st2ws.wsarea[nwsindex]));
		sitegrid->SetCellValue(r, 2, std::to_string(g_st2ws.stdrprun[a.first]));
		r++;
	}

	cellgrid = new wxGrid(this, -1, wxPoint(0, 100), wxSize(400, 300));
	cellgrid->CreateGrid(2, 2);
	cellgrid->SetCellValue(0, 0, "CELL");
	cellgrid->SetCellValue(0, 1, "VALUE");
	cellgrid->SetCellValue(1, 0, g_st2ws.vwscd[nwsindex]);
	cellgrid->SetCellValue(1, 1, std::to_string(g_st2ws.wsdrp[nwsindex]));
}
