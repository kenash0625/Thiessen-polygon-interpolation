
#define MINDIST 10

#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/geom/Point.h>
#include <algorithm>
#include<sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "st2ws.h"
#include "OGRFile.h"

#ifdef GDALDEBUG
#pragma comment(lib,"../debug/geos_i_361_msvc1800.lib")
#else 
#pragma comment(lib,"../release/geos_i_361_msvc1800.lib")
#endif

double great_circle_dist(double fi1/*y*/, double lam1, double fi2, double lam2);
void set_useable_station(int stindex);
//vector<double> wsx, wsy, wsdrp, stx, sty, stdrp, wsarea;
//set<string> vwscd, vstcd;
//vector<double> stxrun, styrun, stdrprun;
//vector<string> vstcdrun;
//vector<int> stidxrun;

//int useablest, setrainindex;
//useable_station _usable_station;
//station_interp _station_interp;
string filename, _wsfile = "WATA.shp";
OGRFile ofilews;


map<string, st2ws_st> stcds2;
map<string, st2ws_st> &get_stcds2() { return stcds2; }
map<string, st2ws_ws> wscds2;
//使用这些站点时 各个流域的weight
map<string, map<string, st2ws_ws>> vweight;

void add_st(string &stcd, double x, double y)
{
	stcds2[stcd].x = x;
	stcds2[stcd].y = y;
}
//放到vector
void add_ws(string &wscds)
{
	wscds2[wscds];
}
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
		else if (ls0.Intersects(&ls1))// 0123->0x21x3
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

void SelfIntersectGeom(const string &fn)
{
	OGRFile ooo(fn, OGRFile::app);
	for (OGRFeature *p; p = ooo.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(p))
	{
		OGRGeometry *pg = p->GetGeometryRef();
		SelfIntersectGeom(pg);
		p->SetGeometry(pg);
	}
}

void set_ws_shp(const string &spath)
{
	ofilews.init(spath);
	ofilews.m_pLayer->ResetReading();
	for (OGRFeature *pFeature; pFeature = ofilews.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		string cd = std::to_string(pFeature->GetFID());
		auto iter = wscds2.find(cd);
		if (iter == wscds2.end())
		{
			continue;
		}

		OGRGeometry *pGeom = pFeature->GetGeometryRef();
		OGRwkbGeometryType gt = pGeom->getGeometryType();
		if (gt == wkbPolygon)
		{
			iter->second.area = ((OGRPolygon *)pGeom)->get_Area();
		}
		else if (gt == wkbMultiPolygon)
		{
			iter->second.area = ((OGRMultiPolygon *)pGeom)->get_Area();
		}
		//SelfIntersectGeom(pGeom);
		//geos::geom::Geometry* pp = (geos::geom::Geometry*)pGeom->exportToGEOS(ofilews.geosctx());
		//vwsgeom[idx] = pp;		
	}
	string sql = string("CREATE SPATIAL INDEX ON ") + ofilews.m_pLayer->GetName();
	ofilews.m_poDS->ExecuteSQL(sql.c_str(), nullptr, nullptr);
	ofilews.open();
	vweight[""] = wscds2;
	ofstream ofs(filename + ".nodata");
	for (auto &a : wscds2)
		if (a.second.area == 0)
			ofs << a.first << '\t';
	ofs << endl;
}
//循环设置所有测站雨量
void begin_st_rain()
{
	for (auto &a : stcds2)
	{
		a.second.z = -1;
		if (a.second.pGeom)
		{
			OGRGeometryFactory::destroyGeometry(a.second.pGeom);
			a.second.pGeom = nullptr;
		}
	}
}
void st_rain(const string &cd, double d)
{
	auto iter = stcds2.find(cd);
	if (iter != stcds2.end())
	{
		iter->second.z = d;
	}
}

void ThiessenCalcZ2(map<string, st2ws_ws> *pweight, int nfind)
{
	if (!nfind)
	{
		*pweight = vweight.at("");
		geos::geom::Geometry **geomRun;
		geos::geom::CoordinateArraySequence seq;
		geos::geom::Envelope env;
		geos::triangulate::VoronoiDiagramBuilder builder;
		const geos::geom::GeometryFactory& geomFact(*geos::geom::GeometryFactory::getDefaultInstance());
		std::auto_ptr<geos::geom::GeometryCollection> polys;
		//sitePoly2 : index - site index ;val - cell index
		OGRFile opt(filename + "st.shp", OGRFile::out, "", "ESRI Shapefile", wkbPoint);
		OGRFieldDefn odfn("STCD", OFTString);
		opt.m_pLayer->CreateField(&odfn);
		for (auto &a : stcds2)
		{
			if (a.second.z < 0) continue;
			geos::geom::Coordinate coordrun(a.second.x, a.second.y);
			seq.add(coordrun);
			env.expandToInclude(coordrun);
			OGRFeature of(opt.m_pLayer->GetLayerDefn());
			OGRPoint oopt(a.second.x, a.second.y);
			of.SetGeometry(&oopt);
			of.SetField("STCD", a.first.c_str());
			opt.m_pLayer->CreateFeature(&of);
		}
		opt.close();
		OGREnvelope oenv;
		ofilews.m_pLayer->GetExtent(&oenv);
		env.expandToInclude(oenv.MinX, oenv.MinY);
		env.expandToInclude(oenv.MaxX, oenv.MaxY);

		builder.setSites(seq);
		builder.setTolerance(0);
		builder.setClipEnvelope(&env);
		polys = builder.getDiagram(geomFact);
		for (auto &a : stcds2)
		{
			if (a.second.z < 0) continue;
			for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
			{
				const geos::geom::Point *pcoord = (const geos::geom::Point *)(polys->getGeometryN(i)->getUserData());
				if (a.second.x == pcoord->getX() && a.second.y == pcoord->getY())
				{
					a.second.pGeom = OGRGeometryFactory::createFromGEOS(ofilews.geosctx(), (GEOSGeom)polys->getGeometryN(i));
					break;
				}
			}
		}
		{
			OGRFile o(filename + "voronoi.shp", OGRFile::out, "", "ESRI Shapefile", wkbPolygon);
			OGRFieldDefn odfn("stcd", OFTString);
			o.m_pLayer->CreateField(&odfn);
			for (auto &a : stcds2)
			{
				if (a.second.z < 0) continue;
				OGRFeature of(o.m_pLayer->GetLayerDefn());
				of.SetField(odfn.GetNameRef(), a.first.c_str());
				of.SetGeometry(a.second.pGeom);
				o.m_pLayer->CreateFeature(&of);
			}
		}
		for (auto &a : stcds2)
		{
			if (a.second.z < 0) continue;
			ofilews.m_pLayer->SetSpatialFilter(nullptr);
			ofilews.m_pLayer->SetSpatialFilter(a.second.pGeom);
			int length = ofilews.m_pLayer->GetFeatureCount();
			OGRFeature *poFea;
			for (size_t i = 0; i < length; i++, OGRFeature::DestroyFeature(poFea))
			{
				poFea = ofilews.m_pLayer->GetNextFeature();
				if (poFea == nullptr)
				{
					continue;//按照getnextfeature取得的图元为和过滤器相交的. 但是总数不等于getfeaturecount
				}
				string cd = poFea->GetFieldAsString("WSCD");
				if (wscds2.find(cd) == wscds2.end())
				{
					continue;
				}
				if (poFea->GetGeometryRef()->Within(a.second.pGeom) || poFea->GetGeometryRef()->Contains(a.second.pGeom))
				{
					(*pweight)[cd].stcdweight[a.first] = wscds2.at(cd).area;
				}
				else if (poFea->GetGeometryRef()->Intersects(a.second.pGeom))
				{
					OGRGeometry *pGeom = poFea->GetGeometryRef()->Intersection(a.second.pGeom);

					if (pGeom->getGeometryType() == wkbMultiPolygon)
					{
						(*pweight)[cd].stcdweight[a.first] = ((OGRMultiPolygon *)pGeom)->get_Area();
					}
					else if (pGeom->getGeometryType() == wkbPolygon)
					{
						(*pweight)[cd].stcdweight[a.first] = ((OGRPolygon *)pGeom)->get_Area();
					}
					delete pGeom;
				}
			}
		}
		for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
		{
			geos::geom::Point *pcoord = (geos::geom::Point *)polys->getGeometryN(i)->getUserData();
			geomFact.destroyGeometry(pcoord);
		}
	}
	for (auto &a : *pweight)
	{
		double dsum = 0;
		for (auto &aa : a.second.stcdweight)
		{
			dsum += aa.second*stcds2.at(aa.first).z;
		}
		a.second.p = dsum / a.second.area;
	}

}

map<string, st2ws_ws>& calc(const string &a, int &useablest, double &d)
{
	filename = a;
	useablest = 0;
	d = 0;
	stringstream ss;
	for (auto &a : stcds2)
	{
		if (a.second.z < 0)
		{
			continue;
		}
		useablest++;
		ss << a.first;
		d = a.second.z;
	}
	cout << "useable st =" << useablest << endl;
	if (useablest == 1)
	{
		return vweight[""];
	}
	else if (useablest == 0)
	{
		return vweight[""];
	}
	else {
		string str = ss.str();
		int nfiind = (vweight.end() != vweight.find(str));
		ThiessenCalcZ2(&vweight[str], nfiind);
		return vweight[str];
	}
}

void wsfile(const string &a)
{
	_wsfile = a;
}

st2ws_ws::st2ws_ws() :area(0), p(0)
{
}

st2ws_st::st2ws_st() : x(0), y(0), z(0), pGeom(nullptr)
{
}
map<string, st2ws_ws>& rand_run(const string &r,int &u,double &d)
{
	//for (;;)
	{
		int range_max = 100, range_min = -100, nStsize = stcds2.size();
		srand((unsigned)time(NULL));
		begin_st_rain();
		for (auto &a:stcds2)
		{
			int drun = rand() % 100 - 50;
			a.second.z = drun;
		}
		return calc(r, u, d);
	}
}