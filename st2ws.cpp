// Staion2PlainRain.cpp : 定义 DLL 应用程序的导出函数。

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
#include <vector>
#include <map>
#include "../MathAlgorithmLib/MathAlgorithmLib.h"
#include "st2ws.h"
#include "json/json.h"
#include <iostream>
#include <iomanip>
#include "OGRFile.h"
using namespace Json;
#ifdef _DEBUG
#pragma comment(lib,"../debug/MathAlgorithmLib.lib")
#else
#pragma comment(lib,"../release/MathAlgorithmLib.lib")

#endif
#ifdef GDALDEBUG
#pragma comment(lib,"../debug/geos_i_361_msvc1800.lib")
#else 
#pragma comment(lib,"../release/geos_i_361_msvc1800.lib")
#endif

void ThiessenCalcZ2(int sitecnt, double *sitex, double *sitey, double *sitez, int cellcnt, geos::geom::Geometry **cells, double*cellarea,double *cellz,vector<map<int,double>> *pweight,int nfind)
{
	vector<int> sitePoly2(sitecnt);
	double *p1, *p2, *p3;
	int idxrun;
	if (!nfind)
	{
		geos::geom::Geometry **geomRun;
		geos::geom::CoordinateArraySequence seq;
		geos::geom::Envelope env;
		geos::triangulate::VoronoiDiagramBuilder builder;
		const geos::geom::GeometryFactory& geomFact(*geos::geom::GeometryFactory::getDefaultInstance());
		std::auto_ptr<geos::geom::GeometryCollection> polys;
		//sitePoly2 : index - site index ;val - cell index
		for (p3 = sitez, p1 = sitex, p2 = sitey, idxrun = 0; idxrun < sitecnt;
		p1++, idxrun++, p2++, p3++)
		{
			geos::geom::Coordinate coordrun(*p1, *p2);
			seq.add(coordrun);
			env.expandToInclude(coordrun);
		}
		for (geomRun = cells, idxrun = 0; idxrun < cellcnt; idxrun++, geomRun++)
		{
			env.expandToInclude((*geomRun)->getEnvelopeInternal());
		}
		builder.setSites(seq);
		builder.setTolerance(0);
		builder.setClipEnvelope(&env);
		polys = builder.getDiagram(geomFact);
		for (p3 = sitez, p1 = sitex, p2 = sitey, idxrun = 0; idxrun < sitecnt;
		p1++, idxrun++, p2++, p3++)
		{
			for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
			{
				const geos::geom::Point *pcoord = (const geos::geom::Point *)polys->getGeometryN(i)->getUserData();
				if (*p1 == pcoord->getX() && *p2 == pcoord->getY())
				{
					sitePoly2[idxrun] = i;
					break;
				}
			}
		}
		//流域下标-测站下标-权重
		pweight->resize(cellcnt);
		for (p2 = cellarea, p1 = cellz, geomRun = cells, idxrun = 0; idxrun < cellcnt; idxrun++, geomRun++, p1++, p2++)
		{
			*p1 = 0;
			int calc = 0;
			//两种情况 流域被一个多边形包含 或 流域与多个多边形相交
			for (std::vector<int>::iterator i = sitePoly2.begin(); i != sitePoly2.end(); i++)
			{
				const geos::geom::Polygon *a = dynamic_cast<const geos::geom::Polygon*>(polys->getGeometryN(*i));
				if (a->equals(*geomRun) || a->contains(*geomRun))
				{
					calc = 1;
					(*pweight)[idxrun][i - sitePoly2.begin()] = *p2;
				}
			}
			if (calc) continue;
			for (std::vector< int>::iterator i = sitePoly2.begin(); i != sitePoly2.end(); i++)
			{
				const geos::geom::Polygon *a = dynamic_cast<const geos::geom::Polygon*>(polys->getGeometryN(*i));
				if ((*geomRun)->intersects(a))
				{
					geos::geom::Geometry *pGeom = (*geomRun)->intersection(a);
					if (!pGeom || (geos::geom::GeometryTypeId::GEOS_POLYGON != pGeom->getGeometryTypeId() && geos::geom::GeometryTypeId::GEOS_MULTIPOLYGON != pGeom->getGeometryTypeId()))
					{
						//ofsLog << endl << "cell.cellBoundary->Intersection(&a)==0";
						throw 1;
					}
					(*pweight)[idxrun][i - sitePoly2.begin()] = pGeom->getArea();
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
	for (p2=cellz,p1=cellarea,idxrun = 0; idxrun < cellcnt;idxrun++,p1++,p2++)
	{
		double dsum = 0;
		for (auto &a:pweight->at(idxrun))
		{
			dsum += sitez[a.first] * a.second;
		}
		*p2 = dsum / *p1;
	}
}

double st2ws::great_circle_dist(double fi1/*y*/, double lam1, double fi2, double lam2)
{
	static const double _2pai = 3.1415926 / 180;
	fi1 *= _2pai;
	fi2 *= _2pai;
	lam1 *= _2pai; lam2 *= _2pai;
	double deltafi = fabs(fi1 - fi2);
	double deltalam = fabs(lam1 - lam2);
	double tao = atan2(sqrt(pow(cos(fi2)*sin(deltalam), 2) + pow(cos(fi1)*sin(fi2) - sin(fi1)*cos(fi2)*cos(deltalam), 2)) , (sin(fi1)*sin(fi2) + cos(fi1)*cos(fi2)*cos(deltalam)));
	return fabs(6378.137*tao);
}
void st2ws::set_useable_station(int stindex)
{
	styrun[useablest] = sty[stindex];
	stdrprun[useablest] = stdrp[stindex];
	stxrun[useablest] = stx[stindex];
	vstcdrun[useablest] = vstcd[stindex];
	stidxrun[useablest] = stindex;
	useablest++;
}

void st2ws::fstationinterp(station_interp _stationInterp)
{
	if (useablest == 1)
	{
		for (double &d : wsdrp) d = *stdrprun.begin();
	}
	else if (useablest == 0)
	{
		for (double &d : wsdrp) d = 0;
	}
	else {

		switch (_stationInterp)
		{
		case st2ws::KRIGING:
			KrigingCalcZ(wsx.size(), wsx.data(), wsy.data(), wsdrp.data(), useablest, stxrun.data(), styrun.data(), stdrprun.data());
			break;
		case st2ws::THIESSEN:	
		{
			stringstream ss;
			for (int i = 0; i < useablest;i++)
			{
				ss << stidxrun[i];
			} 
			string str = ss.str();
			int nfiind = (vweight.end() != vweight.find(str));
			ThiessenCalcZ2(useablest, stxrun.data(), styrun.data(), stdrprun.data(), wsx.size(), vwsgeom.data(), wsarea.data(),wsdrp.data(),&vweight[str],nfiind);
			break;
		}
		case st2ws::INVERSE:
			InverseDistCalcZ(wsx.size(), wsx.data(), wsy.data(), wsdrp.data(), useablest, stxrun.data(), styrun.data(), stdrprun.data());
			break;
		case st2ws::LINER:
			LinearityCalcZ(wsx.size(), wsx.data(), wsy.data(), wsdrp.data(), useablest, stxrun.data(), styrun.data(), stdrprun.data());
			break;
		default:
			break;
		}
	}
}

void st2ws::add_st(double x, double y, const string &cd)
{
	stx.push_back(x);
	sty.push_back(y);
	vstcd.push_back(cd);
}
void st2ws::add_ws(double x, double y,const string &cd)
{
	wsx.push_back(x);
	wsy.push_back(y);
	vwscd.push_back(cd);
}
void st2ws::add_end()
{
	stdrp.resize(stx.size());
	wsdrp.resize(wsx.size());
	stxrun.resize(stx.size());
	styrun.resize(stx.size());
	stdrprun.resize(stx.size());
	vstcdrun.resize(stx.size());
	vwsgeom.resize(wsx.size());
	stidxrun.resize(stx.size());
	wsarea.resize(wsx.size());
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
void st2ws::ws_shp(const string & spath)
{
	OGRFile ofilews(spath, OGRFile::in);
	ofilews.m_pLayer->ResetReading();
	//这个边框应包括site cell 
	for (OGRFeature *pFeature; pFeature = ofilews.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		string cd = pFeature->GetFieldAsString("WSCD");
		OGRGeometry *pGeom = pFeature->GetGeometryRef();
		SelfIntersectGeom(pGeom);
		int idx = find(vwscd.begin(), vwscd.end(), cd) - vwscd.begin();
		geos::geom::Geometry* pp=(geos::geom::Geometry*)pGeom->exportToGEOS(ofilews.geosctx());
		vwsgeom[idx] = pp;
		wsarea[idx] = pp->getArea();
	}
}
void st2ws::begin_st_rain()
{
	useablest= setrainindex = 0;
}
void st2ws::st_rain(double d)
{
	stdrp[setrainindex++] = d;
}
void st2ws::calc(useable_station u, station_interp s)
{
	fuseablestation(_usable_station=u);
	fstationinterp(_station_interp=s);
}
vector<double>& st2ws::ws_rain()
{
	return wsdrp;
}
void st2ws::fuseablestation(useable_station _useableStation)
{
	vector<int> nopst;
	size_t i;
	for (i=0;i<stdrp.size();i++)
	{
		if (stdrp[i] < 0)
		{
			if (_useableStation == DYNAMIC) continue;
			else if (_useableStation == STATIC0)  stdrp[i] = 0;
			else if (_useableStation == STATICINVERSE)
			{
				nopst.push_back(i);
				continue;
			}
		}
		set_useable_station(i);
	}
	if (_useableStation == STATICINVERSE)
	{
		int nfind = 0;
		vector<double> stdrpadd;
		vector<int> stadd;
		for (int nindx:nopst)
		{
			vector<double> stxtmp, stytmp, stdrptmp;
			for (i = 0; i < useablest; i++)
			{
				double dist = great_circle_dist(sty[nindx], stx[nindx], styrun[i], stxrun[i]);
				if (dist > MINDIST) continue;
				stxtmp.push_back(stxrun[i]);
				stytmp.push_back(styrun[i]);
				stdrptmp.push_back(stdrprun[i]);
			}
			double d=-1;
			InverseDistCalcZ(1, &stx[nindx], &sty[nindx], &d, stxtmp.size(), stxtmp.data(), stytmp.data(), stdrptmp.data());
			if (d < 0) d = 0;
			if (!stxtmp.empty()) nfind++;
			stdrpadd.push_back(d);
			stadd.push_back(nindx);
		}
		cout << nopst.size() << '\t' << nfind;
		for (i = 0; i < stadd.size(); i++)
		{
			stdrp[stadd[i]] = stdrpadd[i];
			set_useable_station(stadd[i]);
		}
	}
}

string retttt;
char * out_put_voronoi(char *spath, char * wsshp, char * stshp, char ** stcds, int stcnt, char ** wscds, int wscnt)
{
//	retttt.clear();
//	OGRFile ostfile(stshp);
//	map<string, MySites> sites;
//	int cnt(0);
//	for (int index = 0; index < stcnt; index++) { sites[stcds[index]]; }
//	for (OGRFeature *pFeature; pFeature = ostfile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
//	{
//		// 		string strsql = string("STCD = '") + stcds[index] + "'";
//		// 		ostfile.m_pLayer->SetAttributeFilter(NULL);
//		string cd = pFeature->GetFieldAsString("STCD");
//		if (sites.find(cd) == sites.end()) continue;
//
//		OGRPoint *popt = (OGRPoint *)pFeature->GetGeometryRef();
//		MySites &rsecond = sites[cd];
//		rsecond.sCD = cd;
//		rsecond.oPt.setX(popt->getX());
//		rsecond.oPt.setY(popt->getY());
//		cnt++;
//	}
//	if (cnt < 2)
//	{
//		return (char *)retttt.c_str();
//	}
//	string rspath(spath);
//	if (rspath.size() > 0 && *rspath.rbegin() != '\\' && *rspath.rbegin() != '/')
//	{
//		rspath += '\\';
//	}
//	vector<string> vv;
//	for (int i = 0; i < wscnt; i++)
//	{
//		vv.push_back(wscds[i]);
//	}
//	ThiessenInterp thiInterp;
//	thiInterp.Pts(wsshp, &vv);
//	map<string, MyCells> *pCells = thiInterp.Weight(sites, 1);
//	//output Voronoi.shp
//	OGRFile oVoronoi(rspath + "Voronoi.shp", OGRFile::out, "", "ESRI Shapefile", wkbPolygon, ostfile.m_pLayer->GetSpatialRef());
//	OGRFieldDefn oVoronoiFld("STCD", OFTString);
//	oVoronoi.addfld(oVoronoiFld);
//	for (auto &a : sites)
//	{
//		OGRFeature ofadd(oVoronoi.m_pLayer->GetLayerDefn());
//		ofadd.SetGeometry(&a.second.voroPoly);
//		ofadd.SetField("STCD", a.second.sCD.c_str());
//		oVoronoi.m_pLayer->CreateFeature(&ofadd);
//	}
//	//output VoronoiPt.shp
//	OGRFile oVoronoiPt(rspath + "VoronoiPt.shp", OGRFile::out, "", "ESRI Shapefile", wkbPoint, ostfile.m_pLayer->GetSpatialRef());
//	OGRFieldDefn oVorptSt("STCD", OFTString);
//	OGRFieldDefn oVorptWs("WSCD", OFTString);
//	OGRFieldDefn oVorptVal("VALUE", OFTReal);
//	oVoronoiPt.addfld(oVorptSt);
//	oVoronoiPt.addfld(oVorptWs);
//	oVoronoiPt.addfld(oVorptVal);
//	for (auto &a : *pCells)
//	{
//		for (auto &b : a.second.siteArea)
//		{
//			OGRFeature ofadd(oVoronoiPt.m_pLayer->GetLayerDefn());
//			ofadd.SetGeometry(&a.second.siteCenter.at(b.first));
//			ofadd.SetField("STCD", b.first.c_str());
//			ofadd.SetField("WSCD", a.first.c_str());
//			ofadd.SetField("VALUE", b.second / a.second.dF);
//			oVoronoiPt.m_pLayer->CreateFeature(&ofadd);
//		}
//	}
//
//	Json::Value val;
//	for (int index = 0; index < wscnt; index++)
//	{
//		if (pCells->find(wscds[index]) != pCells->end())
//		{
//			for (auto &a : (*pCells)[wscds[index]].siteArea)
//			{
//				Json::Value valrun;
//				valrun["WSCD"] = wscds[index];
//				string aaa = a.first;
//				valrun["STCD"] = aaa;
//				stringstream ssrun;
//				ssrun << fixed << setprecision(8) << a.second / (*pCells)[wscds[index]].dF;//jsoncpp:同样是数字 有的带引号有的不带?
//				valrun["WEIGHT"] = ssrun.str();
//				val.append(valrun);
//			}
//		}
//	}
//	Json::FastWriter wrt;
//	retttt = wrt.write(val);
	char *xml = (char *)retttt.c_str();
	return xml;
}