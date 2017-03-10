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
#include <sstream>
#include <iostream>
#include <iomanip>
#include "OGRFile.h"
#include "st2ws.h"

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
		stringstream ss;
		for (int i = 0; i < useablest; i++)
		{
			ss << stidxrun[i];
		}
		string str = ss.str();
		int nfiind = (vweight.end() != vweight.find(str));
		ThiessenCalcZ2(useablest, stxrun.data(), styrun.data(), stdrprun.data(), wsx.size(), vwsgeom.data(), wsarea.data(), wsdrp.data(), &vweight[str], nfiind);
	}
}

void st2ws::add_st(double x, double y, const string &cd)
{
	stx.push_back(x);
	sty.push_back(y);
	vstcd.push_back(cd);
}
void st2ws::add_ws(geos::geom::Geometry*p,const string &cd)
{
	vwsgeom.push_back(p);
	vwscd.push_back(cd);
	wsarea.push_back(p->getArea());
}
void st2ws::add_end()
{
	stdrp.resize(stx.size());
	wsdrp.resize(wsx.size());
	stxrun.resize(stx.size());
	styrun.resize(stx.size());
	stdrprun.resize(stx.size());
	vstcdrun.resize(stx.size());
	stidxrun.resize(stx.size());
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
			//InverseDistCalcZ(1, &stx[nindx], &sty[nindx], &d, stxtmp.size(), stxtmp.data(), stytmp.data(), stdrptmp.data());
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
