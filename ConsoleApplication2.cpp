// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".

#include "OGRFile.h"
#include "ConsoleApplication2.h"
#include "st2ws.h"
#include <fstream>
#include <thread>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LineString.h>
#include <geos/operation/buffer/BufferOp.h>
#include <wx/headerctrl.h>
#include <wx/sizer.h>
#include <wx/gbsizer.h>
#include <wx/splitter.h>
#include <wx/stdpaths.h>
#include <wx/filename.h>
struct st2ws g_st2ws;
/*int main(int, wchar_t*[])
{
	return WinMain(::GetModuleHandle(NULL), NULL, NULL, SW_SHOWNORMAL);
}
*/
enum
{
	ID_Hello = 1,
	ID_RandRun,
	ID_ZoomIn,
	ID_ZoomOut,
	ID_Pan,
};

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_MENU(ID_RandRun, MyFrame::OnRandRun)
EVT_MENU(ID_Hello, MyFrame::OnHello)
EVT_MENU(ID_ZoomIn, MyFrame::ZoomIn)
EVT_MENU(ID_ZoomOut, MyFrame::ZoomOut)
EVT_MENU(ID_Pan, MyFrame::Pan)
wxEND_EVENT_TABLE()

wxBEGIN_EVENT_TABLE(MyCanvas, wxScrolledWindow)
EVT_PAINT(MyCanvas::OnPaint)
EVT_SIZE(MyCanvas::OnSize)
EVT_LEFT_DOWN(MyCanvas::OnMouseLDown)
EVT_LEFT_UP(MyCanvas::OnMouseLUp)
EVT_RIGHT_DOWN(MyCanvas::OnMouseRDown)

wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(MyApp);

bool MyApp::OnInit()
{	
	frame = new MyFrame();
	frame->Show(true);
	return true; 
}
MyFrame::MyFrame() : wxFrame(NULL, wxID_ANY, "Sample", wxDefaultPosition, wxSize(1366, 768))
{
	wxMenu *menuFile = new wxMenu;
	menuFile->Append(ID_Hello, "&ReadShp...\tCtrl-H", "Help string shown in status bar for this menu item");
	menuFile->Append(ID_RandRun, "&RandomThiessenInterp...\tCtrl-R", "Help string shown in status bar for this menu item");
	menuFile->AppendSeparator();
	menuFile->Append(ID_ZoomIn, "&ZoomIn...\tCtrl-I", "Zoom In");
	menuFile->Append(ID_ZoomOut, "&ZoomOut...\tCtrl-O", "Zoom Out");
	menuFile->Append(ID_Pan, "&Pan...\tCtrl-P", "Pan");
	menuFile->Append(wxID_EXIT);
	wxMenu *menuHelp = new wxMenu;
	menuHelp->Append(wxID_ABOUT);
	wxMenuBar *menuBar = new wxMenuBar;
	menuBar->Append(menuFile, "&File");
	menuBar->Append(menuHelp, "&Help");
	SetMenuBar(menuBar);
	CreateStatusBar();
	SetStatusText("Welcome");

	wxSplitterWindow *mySplitter = new wxSplitterWindow(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_3D | wxSP_LIVE_UPDATE | wxCLIP_CHILDREN /* | wxSP_NO_XP_THEME */),
		*gridSplitter = new wxSplitterWindow(mySplitter, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_3D | wxSP_LIVE_UPDATE | wxCLIP_CHILDREN /* | wxSP_NO_XP_THEME */);
	canvas = new MyCanvas(mySplitter);
	cellgrid = new wxGrid(gridSplitter, wxID_ANY);
	sitegrid = new wxGrid(gridSplitter, wxID_ANY);


	cellgrid->CreateGrid(2, 2);
	cellgrid->SetCellValue(0, 0, "CELL");
	cellgrid->SetCellValue(0, 1, "VALUE");

	sitegrid->CreateGrid(1, 3);
	sitegrid->SetCellValue(0, 0, "SITE");
	sitegrid->SetCellValue(0, 1, "WEIGHT");
	sitegrid->SetCellValue(0, 2, "VALUE");

	//wxGridBagSizer *gbSizer = new wxGridBagSizer();
	//gbSizer->Add(canvas, wxGBPosition(0,0), wxGBSpan(2,2), wxALIGN_CENTER | wxALL, 5);
	//gbSizer->Add(cellgrid, wxGBPosition(0, 2), wxGBSpan(1, 1), wxALIGN_LEFT | wxALL, 2);
	//gbSizer->Add(sitegrid, wxGBPosition(1, 2), wxGBSpan(1, 1), wxALIGN_LEFT | wxALL, 2);
	//wxGridSizer *gridsizer = new wxGridSizer(1, 5, 5);
	//gridsizer->Add(cellgrid,wxSizerFlags().Align(wxALIGN_CENTER_VERTICAL));
	//gridsizer->Add(sitegrid,wxSizerFlags().Align(wxGROW | wxALIGN_CENTER_VERTICAL));
	//wxGridSizer *gsizer = new wxGridSizer(2, 5,5 );
	//gsizer->Add(canvas, wxSizerFlags().Align(wxALIGN_TOP));
	//gsizer->Add(gridsizer);

	//SetSizerAndFit(gbSizer);
	//Centre();

	// If you use non-zero gravity you must initialize the splitter with its
	// correct initial size, otherwise it will change the sash position by a
	// huge amount when it's resized from its initial default size to its real
	// size when the frame lays it out. This wouldn't be necessary if default
	// zero gravity were used (although it would do no harm neither).
	mySplitter->SetSize(GetClientSize());
	mySplitter->SetSashGravity(1.0);
	gridSplitter->SetSize(GetClientSize());
	gridSplitter->SetSashGravity(1.0);

	gridSplitter->Initialize(sitegrid);
	gridSplitter->SplitHorizontally(sitegrid, cellgrid, 100);
	gridSplitter->SetSashPosition(GetClientSize().GetHeight()*0.7);
	mySplitter->Initialize(canvas);
	mySplitter->SplitVertically(canvas, gridSplitter, 500);
	mySplitter->SetSashPosition(GetClientSize().GetWidth()*0.7);
}
void MyFrame::OnHello(wxCommandEvent& event)
{
	canvas->SetExtent();
}
void MyFrame::OnRandRun(wxCommandEvent& event)
{
	g_st2ws.rand_run();
	canvas->SetExtent();
}

void MyFrame::ZoomIn(wxCommandEvent &)
{
	canvas->mode = MyCanvas::Mode::ZoomIn;
}

void MyFrame::ZoomOut(wxCommandEvent &)
{
	canvas->mode = MyCanvas::Mode::ZoomOut;
}

void MyFrame::Pan(wxCommandEvent &)
{
	canvas->mode = MyCanvas::Mode::Pan;
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

MyCanvas::MyCanvas(wxWindow *parent):wxScrolledWindow(parent, -1, wxPoint(0,0), wxSize(300,400), wxFULL_REPAINT_ON_RESIZE),zfac(1),wsident(-1),mode(Mode::Pan)
{
	wxString strAppPath;
	wxStandardPathsBase& stdp = wxStandardPaths::Get();
	wxFileName exeFile(stdp.GetExecutablePath());
	string str = exeFile.GetPath(wxPATH_GET_VOLUME | wxPATH_GET_SEPARATOR).ToStdString();
	OGRFile oStfile(str+"/test/JCZD.shp");
	OGRFile oWsfile(str+"/test/WATA.shp");	
	if(!oStfile || !oWsfile)
	{
		cout<<str<<endl;
		return;
	}
	OGRFeature *pFeature;
	int nRainsz=0,nStsize;
	oWsfile.m_pLayer->GetExtent(&m_Extent);	
	oStfile.m_pLayer->GetExtent(&m_layerEnv);
	m_layerEnv.Merge(m_Extent);
	oWsfile.m_pLayer->ResetReading();
	oStfile.m_pLayer->ResetReading();
	for (; pFeature = oWsfile.m_pLayer->GetNextFeature();OGRFeature::DestroyFeature(pFeature))
	{
		OGRGeometry *pGeom = pFeature->GetGeometryRef();
		geos::geom::Geometry *pGeosPoly  = (geos::geom::Geometry*)pGeom->exportToGEOS(oWsfile.geosctx());
		g_st2ws.add_ws(pGeosPoly, std::to_string(pFeature->GetFID()));
		
		vector<int> ints;
		vector<OGRRawPoint> pts;
		extractPolygon(pGeom, ints, pts);
		cells.emplace_back(ints);
		cellpts.emplace_back(pts);
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
	if (mode==Mode::ZoomIn)
	{
		Zoom(-0.5);
	}
	else if (mode==Mode::ZoomOut)
	{
		Zoom(0.5);
	}
	else if (mode == Mode::Pan) {
		pt2 = event.GetPosition();
		OGREnvelope rWorld(this->m_Extent);
		OGRRawPoint p1 = GetWorld(GetClientRect(), pt1),
			p2 = GetWorld(GetClientRect(), pt2);
		double dx = p1.x - p2.x, dy = p1.y - p2.y;
		rWorld.MinX += dx;
		rWorld.MaxX += dx;
		rWorld.MinY += dy;
		rWorld.MaxY += dy;
		m_Extent = rWorld;
		SetExtent();
	}
}

void MyCanvas::Zoom(double tfac)
{
	OGRRawPoint opt = GetWorld(GetClientRect(), pt1);
	OGREnvelope rWorld(m_Extent);
	double dx = opt.x - (rWorld.MinX + rWorld.MaxX) / 2.0,
		dy = opt.y - (rWorld.MaxY + rWorld.MinY) / 2.0;
	rWorld.MinX += dx;
	rWorld.MinY += dy;
	rWorld.MaxX += dx;
	rWorld.MaxY += dy;
	
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
	OGRRawPoint optident = GetWorld(GetClientRect(), ptident);
	const geos::geom::GeometryFactory *geomfac(geos::geom::GeometryFactory::getDefaultInstance());
	geos::geom::Coordinate coord(optident.x, optident.y);
	geos::geom::Point *geospt = geomfac->createPoint(coord);
	wsident = 0;
	for (geos::geom::Geometry*p : g_st2ws.vwsgeom)
	{
		if (p->contains(geospt))
		{		
			break;
		}
		++wsident;
	}

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
	memdc.SetBackground(wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_SOLID));
	memdc.Clear();
	//all cells
	list<vector<int>>::iterator itpolys = cells.begin();
	list<vector<OGRRawPoint>>::iterator itpts = cellpts.begin();
	vector<wxPoint> wxpts;
	wxpts.reserve(std::max_element(cellpts.begin(), cellpts.end(), [](vector<OGRRawPoint> &a, vector<OGRRawPoint> &b)->bool {return a.size() < b.size(); })->size());
	for (; itpolys != cells.end(); itpolys++, itpts++)
	{
		wxpts.resize(itpts->size());
		std::transform(itpts->begin(), itpts->end(), wxpts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; XYWorld2DC(&p, &opt); return p; });
		if(itpolys->size()>1)
		{
			memdc.DrawPolyPolygon(itpolys->size(), itpolys->data(), wxpts.data(), 0, 0, wxWINDING_RULE);
		}
		else if(itpolys->size()==1)
		{
			memdc.DrawPolygon(wxpts.size(),wxpts.data());
		}
		
	}
	//identified cell
	if(wsident>=0 && wsident<g_st2ws.vwsgeom.size())
	{
		wxPen pen(wxColour(0, 0, 255), 2, wxPENSTYLE_SOLID);
		memdc.SetPen(pen);
		vector<int> polys;
		vector<OGRRawPoint> polypts;
		extractPolygon(g_st2ws.vwsgeom[wsident], polys, polypts);
		vector<wxPoint> wxpolypts(polypts.size());
		std::transform(polypts.begin(), polypts.end(), wxpolypts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; XYWorld2DC(&p, &opt); return p; });
		if(polys.size()>1)
		{
			memdc.DrawPolyPolygon(polys.size(), polys.data(), wxpolypts.data(), 0, 0, wxPolygonFillMode::wxWINDING_RULE);
		}
		else if(polys.size()==1)
		{
			memdc.DrawPolygon(wxpolypts.size(), wxpolypts.data(), 0, 0, wxPolygonFillMode::wxWINDING_RULE);
		}
	}
	//all sites
	wxBrush p1brush(wxColour(0, 0, 0), wxBRUSHSTYLE_SOLID);
	memdc.SetBrush(p1brush);
	memdc.SetPen(pen);
	vector<OGRRawPoint>::iterator iter = sitecoords.begin();
	for (; iter != sitecoords.end(); iter++)
	{
		wxPoint p;
		XYWorld2DC(&p, &*iter);
		memdc.DrawCircle(p, 4);
	}

	if (g_st2ws.useablest > 1) {
		//voronoi cells
		memdc.SetBrush(brush);
		wxPen pen(wxColour(255, 128, 64), 2, wxPENSTYLE_SOLID);
		memdc.SetPen(pen);
		vector<wxPoint> wxpts(g_st2ws.weightrun->voropts.size());
		std::transform(g_st2ws.weightrun->voropts.begin(), g_st2ws.weightrun->voropts.end(), wxpts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; XYWorld2DC(&p, &opt); return p; });
		if(g_st2ws.weightrun->voropoly.size()>1)
		{	
			memdc.DrawPolyPolygon(g_st2ws.weightrun->voropoly.size(), g_st2ws.weightrun->voropoly.data(), wxpts.data(), 0, 0, wxWINDING_RULE);
		}
		else if(g_st2ws.weightrun->voropoly.size()==1)
		{
			memdc.DrawPolygon(wxpts.size(), wxpts.data(), 0, 0, wxWINDING_RULE);
		}
		//voronoi sites
		wxPen rpen(wxColour(255, 0, 0), 1, wxPENSTYLE_SOLID);
		wxBrush p2brush(wxColour(255, 0, 0), wxBRUSHSTYLE_SOLID);
		memdc.SetBrush(p2brush);
		memdc.SetPen(rpen);
		for (int i = 0; i < g_st2ws.useablest; i++)
		{
			wxPoint p;
			XYWorld2DC(&p, &sitecoords[g_st2ws.stidxrun[i]]);
			memdc.DrawCircle(p, 4);
		}
		if (wsident>=0&&wsident<g_st2ws.weightrun->weights.size())
		{
			wxGrid *sitegrid = wxGetApp().frame->sitegrid,*cellgrid=wxGetApp().frame->cellgrid;
		
			sitegrid->DeleteRows(0,sitegrid->GetNumberRows());
			sitegrid->InsertRows(0,1 + g_st2ws.weightrun->weights[wsident].size());
		
			sitegrid->SetCellValue(0, 0, "SITE");
			sitegrid->SetCellValue(0, 1, "WEIGHT");
			sitegrid->SetCellValue(0, 2, "VALUE");
			int rowrun = 1;
			for (auto &a:g_st2ws.weightrun->weights[wsident])
			{
				sitegrid->SetCellValue(rowrun,0,g_st2ws.vstcdrun[a.first]);
				sitegrid->SetCellValue(rowrun, 1, std::to_string(a.second/g_st2ws.wsarea[wsident]));
				sitegrid->SetCellValue(rowrun, 2, std::to_string(g_st2ws.stdrprun[a.first]));
				rowrun++;
			}
			cellgrid->SetCellValue(1, 0, g_st2ws.vwscd[wsident]);
			cellgrid->SetCellValue(1, 1,  std::to_string(g_st2ws.wsdrp[wsident]));
		}

	}

	memdc.SelectObject(wxNullBitmap);
}


double MyCanvas::XWorld2DC(double x, bool bRound /*= true*/)
{
	x = (x - m_rWorld.MinX) * m_World2DC;

	return(bRound ? (int)(x < 0.0 ? x - 0.5 : x + 0.5) : x);
}
double MyCanvas::YWorld2DC(double y, bool bRound /*= true*/)
{
	y = (m_rWorld.MaxY - y) * m_World2DC - 1;

	return(bRound ? (int)(y < 0.0 ? y - 0.5 : y + 0.5) : y);
}
void MyCanvas::XYWorld2DC(wxPoint *dst, OGRRawPoint *src) {
	dst->x = (int)XWorld2DC(src->x);
	dst->y = (int)YWorld2DC(src->y);
}
double MyCanvas::XDC2World(double x)
{
	return(m_rWorld.MinX + m_DC2World * x);
}
double MyCanvas::YDC2World(double y)
{
	return(m_rWorld.MaxY - m_DC2World * (y + 1));
}
// 图层的纵横比 窗口的纵横比
// 图层的纵横比较小
// 图层的y扩大一些
OGREnvelope MyCanvas::GetWorld(wxRect rClient)
{
	return m_Extent;
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
OGRRawPoint MyCanvas::GetWorld(wxRect rClient, wxPoint ptClient)
{
	return OGRRawPoint(XDC2World(ptClient.x), YDC2World(ptClient.y));
	double		d;
	OGREnvelope	rWorld(GetWorld(rClient));

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

	m_Extent = GetWorld(m_rDC);
	m_rWorld = m_Extent;
	//-----------------------------------------------------
	// ensure cellsize in x-/y-direction are identical...
	//double	dxdyDC = (double)m_rDC.GetWidth() / (double)m_rDC.GetHeight();
	//double	dxdyWorld = (m_rWorld.MaxX - m_rWorld.MinX) / (m_rWorld.MaxY - m_rWorld.MinY);

	//if (dxdyDC > dxdyWorld)
	//{
	//	//m_rWorld.Inflate(0.5 * (m_rWorld.Get_YRange() * dxdyDC - m_rWorld.Get_XRange()), 0.0, false);
	//}
	//else if (dxdyDC < dxdyWorld)
	//{
	//	//m_rWorld.Inflate(0.0, 0.5 * (m_rWorld.Get_XRange() / dxdyDC - m_rWorld.Get_YRange()), false);
	//}

	//-----------------------------------------------------
	m_World2DC = (double)m_rDC.GetWidth() / (m_rWorld.MaxX - m_rWorld.MinX);
	m_DC2World = 1.0 / m_World2DC;

	Hello();
	Refresh();
}


