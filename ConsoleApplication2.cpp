// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".

#include "OGRFile.h"
#include "ConsoleApplication2.h"

#include <fstream>
#include <thread>

#include <wx/headerctrl.h>
#include <wx/sizer.h>
#include <wx/gbsizer.h>
#include <wx/splitter.h>
#include <wx/stdpaths.h>
#include <wx/filename.h>
//#include <occi.h>
//#pragma comment(lib,"oraocci11.lib")
string ExePath;
//int main(int, wchar_t*[])
//{
//	//oracle::occi::Environment::createEnvironment();
//	return WinMain(::GetModuleHandle(NULL), NULL, NULL, SW_SHOWNORMAL);
//}


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

wxIMPLEMENT_APP_CONSOLE(MyApp);

bool MyApp::OnInit()
{
	wxStandardPathsBase& stdp = wxStandardPaths::Get();
	wxFileName exeFile(stdp.GetExecutablePath());
	ExePath = exeFile.GetPath(wxPATH_GET_VOLUME | wxPATH_GET_SEPARATOR).ToStdString();

	OGRFile oStfile(ExePath + "test/JCZD.shp"), oWsfile(ExePath + "test/WATA.shp");
	OGRFeature *pFeature;
	oWsfile.m_pLayer->ResetReading();
	oStfile.m_pLayer->ResetReading();
	for (; pFeature = oWsfile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		string wscd = (pFeature->GetFieldAsString("WSCD"));
		add_ws(wscd);
	}
	for (; pFeature = oStfile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		OGRPoint *pt = (OGRPoint*)pFeature->GetGeometryRef();
		string stcd = pFeature->GetFieldAsString("STCD");
		add_st(stcd, pt->getX(), pt->getY());
	}
	set_ws_shp(ExePath + "test/WATA.shp");
	oStfile.close();
	oWsfile.close();
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
	canvas->RandRun();
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

MyCanvas::MyCanvas(wxWindow *parent):wxScrolledWindow(parent, -1, wxPoint(0,0), wxSize(300,400), wxFULL_REPAINT_ON_RESIZE),zfac(1),mode(Mode::Pan),wscdweight(nullptr)
{
	int u;
	double d;
	wscdweight = &rand_run(ExePath + "test\\a", u, d);

	MyLayer *pst = new MyLayer(ExePath + "test/ast.shp"), *pws = new MyLayer(ExePath + "test/WATA.shp"),*pvor=new MyLayer(ExePath+"test/avoronoi.shp");
	wxPen pen(wxColour(0, 0, 0), 1, wxPENSTYLE_SOLID), identpen(wxColour(255, 0, 0), 2, wxPENSTYLE_SOLID);
	wxBrush brush(wxColour(255, 0, 0), wxBRUSHSTYLE_TRANSPARENT), p1brush(wxColour(0, 0, 0), wxBRUSHSTYLE_SOLID);

	pst->origpen = pen;
	pst->identpen = identpen;
	pst->origbrush = p1brush;
	pst->identbrush = brush;

	pws->origpen = pen;
	pws->identpen = identpen;
	pws->origbrush = wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_TRANSPARENT);
	pws->identbrush = pws->origbrush;

	pvor->origpen = pen;
	pvor->identpen = identpen;
	pvor->origbrush = wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_TRANSPARENT);
	pvor->identbrush = pws->origbrush;

	m_layers.push_back(pst);
	m_layers.push_back(pws);
	m_layers.push_back(pvor);
	pws->ofile.init(pws->ofilepath);
	pws->ofile.m_pLayer->GetExtent(&m_Extent);
	m_layerEnv.Merge(m_Extent);	

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
	OGRPoint opt(optident.x, optident.y);
	//for (MyLayer *p : m_layers)
	{
		m_layers[0]->ofile.init(m_layers[0]->ofilepath);
		MyLayer *p = m_layers[1];
		p->ofile.init(p->ofilepath);
		p->ofile.m_pLayer->SetAttributeFilter(nullptr);
		p->ofile.m_pLayer->ResetReading();
		for (OGRFeature *pFea; pFea = p->ofile.m_pLayer->GetNextFeature();OGRFeature::DestroyFeature(pFea))
		{
			if (pFea->GetGeometryRef()->Contains(&opt))
			{
				p->identfid.clear();
				p->identfid.push_back(pFea->GetFID());
				string cd = pFea->GetFieldAsString("WSCD");
				if (wscdweight==nullptr || wscdweight->find(cd)==wscdweight->end())
				{
					continue;
				}
				wxGrid *sitegrid = wxGetApp().frame->sitegrid, *cellgrid = wxGetApp().frame->cellgrid;
				sitegrid->DeleteRows(0, sitegrid->GetNumberRows());
				sitegrid->InsertRows(0, 1 + wscdweight->at(cd).stcdweight.size());

				sitegrid->SetCellValue(0, 0, "SITE");
				sitegrid->SetCellValue(0, 1, "WEIGHT");
				sitegrid->SetCellValue(0, 2, "VALUE");
				int rowrun = 1;
				m_layers[0]->identfid.clear();
				for (auto &a : wscdweight->at(cd).stcdweight)
				{
					sitegrid->SetCellValue(rowrun, 0, a.first);
					string sql = "STCD='" + a.first + "'";
					m_layers[0]->ofile.m_pLayer->SetAttributeFilter(nullptr);
					m_layers[0]->ofile.m_pLayer->SetAttributeFilter(sql.c_str());
					OGRFeature *pf = m_layers[0]->ofile.m_pLayer->GetNextFeature();
					int cnt = m_layers[0]->ofile.m_pLayer->GetFeatureCount();
					m_layers[0]->identfid.push_back(pf->GetFID());
					OGRFeature::DestroyFeature(pf);
					sitegrid->SetCellValue(rowrun, 1, std::to_string(a.second / wscdweight->at(cd).area));
					sitegrid->SetCellValue(rowrun, 2, std::to_string(get_stcds2().at(a.first).z));
					rowrun++;
				}
				cellgrid->SetCellValue(1, 0, cd);
				cellgrid->SetCellValue(1, 1, std::to_string(wscdweight->at(cd).p));

				break;
			}
		}
		
		m_layers[0]->ofile.close();
		m_layers[1]->ofile.close();
	}
	SetExtent();
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

	wxSize sz = GetClientSize();

	bitmap.Create(sz);
	memdc.SelectObject(bitmap);
	memdc.SetBackground(wxBrush(wxColour(255, 255, 255), wxBRUSHSTYLE_SOLID));
	memdc.Clear();
	for (MyLayer *p:m_layers)
	{
		p->ofile.init(p->ofilepath);
		p->Draw(memdc,*this);
		p->ofile.close();
	}
	memdc.SelectObject(wxNullBitmap);
	Refresh();
}

void MyCanvas::RandRun()
{
	wscdweight = &rand_run(ExePath+"test/a",useablest, d1stval);
	SetExtent();
}


void MyLayer::Draw(wxMemoryDC &memdc,OGRFeature * pFea,MyCanvas &canv)
{
	OGRGeometry *pGeom = pFea->GetGeometryRef();
	OGRwkbGeometryType geomType = pGeom->getGeometryType();

	switch (geomType)
	{
	default:
		break;
	case wkbPoint:
	{
		OGRPoint *opt = (OGRPoint *)pGeom;
		OGRRawPoint orpt(opt->getX(), opt->getY());
		wxPoint p;
		canv.XYWorld2DC(&p, &orpt);
		memdc.DrawCircle(p, 4);
		break;
	}
	case wkbPolygon:
	{
		vector<int> vparts;
		vector<OGRRawPoint> vpts;
		vector<wxPoint> wxpts;

		OGRFile::extractPolygon(pGeom, vparts, vpts);
		wxpts.resize(vpts.size());
		std::transform(vpts.begin(), vpts.end(), wxpts.begin(), [&](OGRRawPoint &opt)->wxPoint {wxPoint p; canv.XYWorld2DC(&p, &opt); return p; });
		if (vparts.size() > 1)
		{
			memdc.DrawPolyPolygon(vparts.size(), vparts.data(), wxpts.data(), 0, 0, wxWINDING_RULE);
		}
		else if (vparts.size() == 1)
		{
			memdc.DrawPolygon(wxpts.size(), wxpts.data());
		}
		break;
	}
	}
}

void MyLayer::Draw(wxMemoryDC &memdc,MyCanvas &canv)
{
	ofile.m_pLayer->SetAttributeFilter(nullptr);
	ofile.m_pLayer->ResetReading();	
	memdc.SetPen(origpen);
	memdc.SetBrush(origbrush);
	for (OGRFeature *pFea; pFea = ofile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFea))
	{
		Draw(memdc, pFea,canv);
	}
	for (int fid : identfid)
	{
		OGRFeature *pFea = ofile.m_pLayer->GetFeature(fid);
		memdc.SetPen(identpen);
		memdc.SetBrush(identbrush);
		Draw(memdc, pFea, canv);
		OGRFeature::DestroyFeature(pFea);
	}
}

MyLayer::MyLayer(const string &str):ofilepath(str)
{
}


