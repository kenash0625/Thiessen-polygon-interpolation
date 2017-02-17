// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".

#include "OGRFile.h"
#include "ConsoleApplication2.h"
#include "ac.h"

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
ac_calibration cal;
MyFrame::MyFrame(wxMDIParentFrame *parent, const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxMDIChildFrame(parent,wxID_ANY,"title test")
{

	canvas = new MyCanvas(this);
	
	
	{
		OGRFile o("C:/Users/awtf/Desktop/a/RIVL.shp");
		o.m_pLayer->GetExtent(&canvas->m_Extent);
		for (OGRFeature *pfe;pfe=o.m_pLayer->GetNextFeature();OGRFeature::DestroyFeature(pfe))
		{
			string cd = pfe->GetFieldAsString("RVCD");
			string frvcd = pfe->GetFieldAsString("FRVCD"), trvcd = pfe->GetFieldAsString("TRVCD");
			
			cal.allnode[cd].cd = cd;
			cal.stringsplit(frvcd, cal.allnode[cd].from);
			cal.allnode[cd].to = trvcd;
		}
		string rt("AGA2100001a00000");
		
		{
			cal.trace(rt, rt);
			cal.wait();
		}

	}
	canvas->SetExtent();
	//canvas->SetScrollbars(10, 10, 100, 240);
}



int main(int ,wchar_t*[])
{
	return WinMain(::GetModuleHandle(NULL), NULL, NULL, SW_SHOWNORMAL);
}

MyCanvas::MyCanvas(MyFrame *parent):wxPanel(parent, -1, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),zfac(1)
{
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
double dfac = 0.5;
void MyCanvas::OnMouseWheel(wxMouseEvent & event)
{
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
	OGRFile ofile("C:/Users/awtf/Desktop/a/RIVL.shp");
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
			this->drdr(&ptls2[i], &ptls[i]);
		}
		memdc.DrawLines(ptls2.size(), ptls2.data());
	}
	memdc.SelectObject(wxNullBitmap);
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