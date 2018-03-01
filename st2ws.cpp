
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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "st2ws.h"
#include "OGRFile.h"

#ifdef GDALDEBUG
#pragma comment(lib,"../debug/geos_i_361_msvc1800.lib")
#else 
#pragma comment(lib,"../release/geos_i_361_msvc1800.lib")
#endif
void SelfIntersectGeom(OGRGeometry *pGeom);
double GreatCircleDist(double fi1/*y*/, double lam1, double fi2, double lam2)
{
	static const double _2pai = 3.1415926 / 180;
	fi1 *= _2pai;
	fi2 *= _2pai;
	lam1 *= _2pai; lam2 *= _2pai;
	double deltafi = fabs(fi1 - fi2);
	double deltalam = fabs(lam1 - lam2);
	double tao = atan2(sqrt(pow(cos(fi2)*sin(deltalam), 2) + pow(cos(fi1)*sin(fi2) - sin(fi1)*cos(fi2)*cos(deltalam), 2)), (sin(fi1)*sin(fi2) + cos(fi1)*cos(fi2)*cos(deltalam)));
	return fabs(6378.137*tao);
}
void InversDistWeight(vector<double> &dist, double unknowX, double unknowY, const vector<double> &knowX, const vector<double> &knowY, double powPara = 1)
{
	vector<double> dist2;
	dist.resize(knowX.size());
	bool b(powPara > 0 && powPara != 1);
	vector<double>::iterator itDist = dist.begin();
	vector<double>::const_iterator itX = knowX.cbegin(), itY = knowY.cbegin();
	for (; itX != knowX.cend(); itDist++, itX++, itY++)
	{
		*itDist = sqrt(pow(*itX - unknowX, 2) + pow(*itY - unknowY, 2));
		if (b) *itDist = pow(*itDist, powPara);
		dist2.push_back(GreatCircleDist(*itX, *itY, unknowX, unknowY));
	}
}
void InversDistWeightZ(double &unknowZ, const vector<double> &dist, const vector<double> &knowZ)
{
	unknowZ = 0;
	double	accumDist = 0;
	vector<double>::const_iterator itDist = dist.cbegin(), itZ = knowZ.cbegin();
	for (; itDist != dist.cend(); itDist++, itZ++)
	{
		unknowZ += *itZ**itDist;
		accumDist += *itDist;
	}
	unknowZ /= accumDist;
}
void NeareastDistWeightZ(double &unknowZ, const vector<double> &dist, const vector<double> &knowZ,double near)
{
	double	nearOne = dist[0];
	int nearInd=0;
	for (int i = 0; i < dist.size(); i++)
	{
		if (dist[i] < nearOne)
		{
			nearOne = dist[i];
			nearInd = i;
		}
	}
	unknowZ = knowZ[nearInd];
	if (nearOne > near) unknowZ = 0;
}
//vector<string> vstcdrun;
//vector<int> stidxrun;

//int useablest, setrainindex;
//useable_station _usable_station;
void ThiessenCalcZ3(double *sitez, int cellcnt, double*cellarea, double *cellz, vector<map<int, double>> &weight	/*流域下标-测站下标-权重*/)
{
	double *p1, *p2;
	int idxrun;
	for (p2 = cellz, p1 = cellarea, idxrun = 0; idxrun < cellcnt; idxrun++, p1++, p2++)
	{
		double dsum = 0;
		for (auto &a : weight.at(idxrun))
		{
			dsum += sitez[a.first] * a.second;
		}
		*p2 = dsum / *p1;
	}
}
//todo 可以把voronoi多边形 和谁相交 相交的面积 保存
void ThiessenCalcWeight(int sitecnt, double *sitex, double *sitey, int cellcnt, vector<string> &cellcd2, double*cellarea,
	const string &ofilename, const string &tms, vector<map<int, double>> &weight	/*流域下标-测站下标-权重*/)
{
	vector<int> sitePoly2(sitecnt);	//sitePoly2 : index - site index ;val - cell index
	vector<OGRGeometry*> ogeos;
	double *p1, *p2, *p3;
	int idxrun;
	OGRFile ofile(ofilename);
	geos::geom::CoordinateArraySequence seq;
	geos::geom::Envelope env;
	geos::triangulate::VoronoiDiagramBuilder builder;
	const geos::geom::GeometryFactory& geomFact(*geos::geom::GeometryFactory::getDefaultInstance());
	std::auto_ptr<geos::geom::GeometryCollection> polys;
	for (p1 = sitex, p2 = sitey, idxrun = 0; idxrun < sitecnt; p1++, idxrun++, p2++)
	{
		geos::geom::Coordinate coordrun(*p1, *p2);
		seq.add(coordrun);
		env.expandToInclude(coordrun);
	}
	OGREnvelope oenv;
	ofile.m_pLayer->GetExtent(&oenv);
	env.expandToInclude(oenv.MinX, oenv.MinY);
	env.expandToInclude(oenv.MaxX, oenv.MinY);
	env.expandToInclude(oenv.MaxX, oenv.MaxY);
	env.expandToInclude(oenv.MinX, oenv.MaxY);

	builder.setSites(seq);
	builder.setTolerance(0);
	builder.setClipEnvelope(&env);
	polys = builder.getDiagram(geomFact);

	//OGRFile ovrfile("./voronoi.shp", OGRFile::out, "", "ESRI Shapefile", wkbPolygon);
	//OGRFile ostfile(tms + "voronoist.shp", OGRFile::out, "", "ESRI Shapefile", wkbPoint);
	//OGRFieldDefn odfn("STCD", OFTString);
	//ovrfile.m_pLayer->CreateField(&odfn);

	for (p1 = sitex, p2 = sitey, idxrun = 0; idxrun < sitecnt; p1++, idxrun++, p2++)
	{
		for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
		{
			const geos::geom::Point *pcoord = (const geos::geom::Point *)(polys->getGeometryN(i)->getUserData());
			if (*p1 == pcoord->getX() && *p2 == pcoord->getY())
			{
				sitePoly2[idxrun] = i;
				//OGRFeature ofadd1(ostfile.m_pLayer->GetLayerDefn());
				//ofadd.SetField(odfn.GetNameRef(),sit)
				//OGRPoint opt(pcoord->getX(), pcoord->getY());
				//ofadd1.SetGeometry(&opt);
				//ostfile.m_pLayer->CreateFeature(&ofadd1);
				break;
			}
		}
	}
	for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
	{
		auto tg = OGRGeometryFactory::createFromGEOS(ofile.geosctx(), (GEOSGeom)polys->getGeometryN(i));
		//OGRFeature ofadd(ovrfile.m_pLayer->GetLayerDefn());
		//ofadd.SetField(odfn.GetNameRef(),sit)
		//ofadd.SetGeometry(tg);
		//ovrfile.m_pLayer->CreateFeature(&ofadd);
		SelfIntersectGeom(tg);
		ogeos.push_back(tg);
	}
	//ovrfile.close();
	//ostfile.close();
	ofile.close();
	weight.resize(cellcnt);
	for (auto &a : weight)
	{
		for (int n = 0; n < sitecnt; n++)
		{
			a[n] = 0;
		}
	}
#pragma omp parallel for
	for (idxrun = 0; idxrun<sitecnt; idxrun++)
	{
		OGRGeometry *stGeom = ogeos[sitePoly2[idxrun]];
		OGRFile ofilea(ofilename);
		ofilea.m_pLayer->SetSpatialFilter(nullptr);
		ofilea.m_pLayer->SetSpatialFilter(stGeom);
		auto length = ofilea.m_pLayer->GetFeatureCount();
		OGRFeature *poFea;
		for (size_t i = 0; i < length; i++, OGRFeature::DestroyFeature(poFea))
		{
			poFea = ofilea.m_pLayer->GetNextFeature();
			if (poFea == nullptr)
			{
				continue;//按照getnextfeature取得的图元为和过滤器相交的. 但是总数不等于getfeaturecount
			}
			string cd = poFea->GetFieldAsString("WSCD");
			vector<string>::iterator iter;
			if ((iter = find(cellcd2.begin(), cellcd2.end(), cd)) == cellcd2.end())
			{
				continue;
			}
			int cellidx = iter - cellcd2.begin();
			if (poFea->GetGeometryRef()->Within(stGeom) )
			{
				weight[cellidx][idxrun] = *(cellarea + cellidx);
			}
			else if (poFea->GetGeometryRef()->Contains(stGeom))
			{
				weight[cellidx][idxrun] = ((OGRPolygon *)stGeom)->get_Area();
			}
			else if (poFea->GetGeometryRef()->Intersects(stGeom))
			{
				OGRGeometry *pGeom = poFea->GetGeometryRef()->Intersection(stGeom);

				if (pGeom->getGeometryType() == wkbMultiPolygon)
				{
					weight[cellidx][idxrun] = ((OGRMultiPolygon *)pGeom)->get_Area();
				}
				else if (pGeom->getGeometryType() == wkbPolygon)
				{
					weight[cellidx][idxrun] = ((OGRPolygon *)pGeom)->get_Area();
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
	for (int i = 0; i < cellcnt; i++)
	{
		double sum(0);
		for (auto &a : weight[i])
		{
			sum += a.second;
		}
		if (sum > cellarea[i])
		{
			cout << (sum - cellarea[i])/ cellarea[i]<< endl;
		}
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

void st2ws::fstationinterp(station_interp _stationInterp,const string &tms)
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
				ss << stidxrun[i]<<"a";
			}
			string str = ss.str();
			auto iter=vweight.find(str);
			if (iter == vweight.end())
			{
				ThiessenCalcWeight(useablest, stxrun.data(), styrun.data(), wsx.size(), vwscd, wsarea.data(), wsfile, tms, vweight[str]);
				iter = vweight.find(str);
			}
			ThiessenCalcZ3(stdrprun.data(), wsx.size(),  wsarea.data(), wsdrp.data(), iter->second);
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
void st2ws::add_ws(double x, double y, double area, const string & cd)
{
	add_ws(x, y, cd);
	wsarea.push_back(area);
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
	int needFind = 1;
	while (needFind) {
		needFind = 0;
		for (int p0 = 0, p1 = 1, p2 = 2, p3 = 3; p3 < vpts.size() && 0==needFind; ++p0, ++p1, ++p2, ++p3)
		{
			OGRLineString ls0, ls1;
			ls0.addPoint(vpts.at(p0).x, vpts.at(p0).y);
			ls0.addPoint(vpts.at(p1).x, vpts.at(p1).y);
			ls1.addPoint(vpts.at(p2).x, vpts.at(p2).y);
			ls1.addPoint(vpts.at(p3).x, vpts.at(p3).y);
			OGRGeometry *pInterGeom = ls0.Intersection(&ls1);
			OGRwkbGeometryType gt;
			OGRPoint *pInterPt = (OGRPoint *)pInterGeom;
			if (pInterGeom && pInterPt
				&& (gt = pInterGeom->getGeometryType()) == wkbPoint)
			{
				if (vpts.at(p1).x == pInterPt->getX() && vpts.at(p1).y == pInterPt->getY())
				{
					vpts.erase(vpts.begin() + p2);
				}
				else if (vpts.at(p2).x == pInterPt->getX() && vpts.at(p2).y == pInterPt->getY())
				{
					vpts.erase(vpts.begin() + p1);
				}
				else if (vpts.at(p0).x == pInterPt->getX() && vpts.at(p0).y == pInterPt->getY())
				{
					vpts.erase(vpts.begin() + p3);
				}
				else if (vpts.at(p3).x == pInterPt->getX() && vpts.at(p3).y == pInterPt->getY())
				{
					vpts.erase(vpts.begin() + p0);
				}
				else {
					OGRRawPoint pt1 = vpts.at(p1), pinter(pInterPt->getX(), pInterPt->getY());
					vpts.at(p1) = vpts.at(p2);
					vpts.at(p2) = pt1;
					vpts.insert(vpts.begin() + p1, pinter);
					vpts.insert(vpts.begin() + p3, pinter);
				}
				needFind = 1;
			}
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
		pPoly->getExteriorRing()->closeRings();
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
			pPoly->getExteriorRing()->closeRings();
		}
	}
}


typedef multimap<long, vector<OGRRawPoint> >::iterator mpIter;
const size_t SECONDPOS = 1;

bool SortPoItVe(vector<mpIter> & vePoIts)
{
	if (vePoIts.size() < 2)
	{//只有一个元素，正常逻辑不需要做顺序修改, 后面也用不到
		return false;
	}

	vector<mpIter> vecTemp(vePoIts);

	while (vecTemp.size() != 1)
	{
		vector<mpIter>::iterator iterA = vecTemp.begin();
		vector<mpIter>::iterator iterB = vecTemp.begin() + 1;
		OGRRawPoint ptA = (*iterA)->second.at(0);
		OGRRawPoint ptB = (*iterB)->second.at(0);
		if (ptA.y > ptB.y)//y越往上越大，保留
		{
			vecTemp.erase(iterB);
		}
		else if (ptA.y < ptB.y)
		{
			vecTemp.erase(iterA);
		}
		else//ptA.y == ptB.y
		{
			if (ptA.x < ptB.x)//越小越靠左，保留
			{
				vecTemp.erase(iterB);
			}
			else
			{
				vecTemp.erase(iterA);
			}
		}
	}

	vePoIts.erase(find(vePoIts.begin(), vePoIts.end(), vecTemp.at(0)));//当然能找到了
	vePoIts.insert(vePoIts.begin(), vecTemp.at(0));

	return true;
}

int UnionWsPts(const vector<vector<OGRRawPoint> > &pos, vector<vector<OGRRawPoint> > &poUnionend)
{
	multimap<long, vector<OGRRawPoint> > allPos;
	static const int iTempKey = 1;// 为了复用代码临时添加，以后修改去除 [6/26/2013 xiarl]
	for (vector<vector<OGRRawPoint> >::size_type index = 0; index < pos.size(); ++index)
	{
		allPos.insert(make_pair(iTempKey, pos.at(index)));
	}

	mpIter iterFind = allPos.begin();
	while (iterFind != allPos.end())
	{
		pair<mpIter, mpIter> findRange;
		findRange = allPos.equal_range(iterFind->first);
		iterFind = findRange.second;

		mpIter iterRun = findRange.first;
		vector<mpIter> sameIter;
		for (; iterRun != findRange.second; ++iterRun)
		{
			sameIter.push_back(iterRun);
		}

		//调序，最左上角多边形至顶,保证其为第一个元素，位置0
		SortPoItVe(sameIter);

		iterRun = sameIter.at(0);
		int indexComp = SECONDPOS;
		while (sameIter.size() > 1)
		{
			mpIter iterComp = sameIter.at(indexComp);
			vector<OGRRawPoint>::iterator iterpt;
			vector<OGRRawPoint>::iterator iterPtFind;
			for (iterpt = iterRun->second.begin();
			iterpt != iterRun->second.end(); ++iterpt)
			{
				iterPtFind = find_if(iterComp->second.begin(), iterComp->second.end(), [&](OGRRawPoint a)->bool {return a.x == iterpt->x &&a.y == iterpt->y; });
				if (iterPtFind != iterComp->second.end())//查找到相同点
				{
					break;
				}
			}
			if (iterpt != iterRun->second.end())//有找到相同点，删除iterComp
			{
				//确定会是闭合的，会有重合点？
				vector<OGRRawPoint> vecPtsTemp;
				vecPtsTemp.insert(vecPtsTemp.end(), iterPtFind, iterComp->second.end() - 1);
				if (iterPtFind != iterComp->second.begin())
				{
					vecPtsTemp.insert(vecPtsTemp.end(), iterComp->second.begin(), iterPtFind);//最后一个闭合重复点
				}

				iterRun->second.insert(iterpt, vecPtsTemp.begin(), vecPtsTemp.end());//最后一个闭合重复点				 

				sameIter.erase(find(sameIter.begin(), sameIter.end(), iterComp));
			}
			else//没找到相同点，临时跳过此位置
			{
				++indexComp;
			}

			if (indexComp == sameIter.size())//一遍过后，重置第二个位置，若sameIter 中的还不止一个，则在下个 while循环中再从第二位置来次
			{// 只剩最后一个的时候 也会进此 if，不过，没啥关系 [3/2/2013 xiarl]
				indexComp = SECONDPOS;
				// 				if (2 == sameIter.size())
				// 				{
				// 					return ERROR_DATA;
				// 				}
			}
		}

		poUnionend.push_back(sameIter.at(0)->second);//合并完成，压入
	}

	return 0;
}
int FindPolyToUni(vector<vector<OGRRawPoint>> &toFind)
{
	for (auto iter = toFind.begin(); iter != toFind.end(); iter++)
	{
		OGRPolygon opoly;
		OGRLinearRing oring;
		oring.setPoints(iter->size(), iter->data());
		oring.closeRings();
		opoly.addRing(&oring);
		for (auto iter2 = iter + 1; iter2 != toFind.end(); iter2++)
		{
			OGRPolygon opoly2;
			OGRLinearRing oring2;
			oring2.setPoints(iter2->size(), iter2->data());
			oring2.closeRings();
			opoly2.addRing(&oring2);
			if (opoly2.Intersects(&opoly))
			{
				vector<vector<OGRRawPoint>> multiPolyToUni, multiPolyUnioned;
				multiPolyToUni.push_back(*iter);
				multiPolyToUni.push_back(*iter2);
				UnionWsPts(multiPolyToUni, multiPolyUnioned);

				toFind.erase(iter2);
				toFind.erase(iter);
				toFind.push_back(multiPolyUnioned[0]);
				return 0;
			}
		}
	}
	return 1;
}
int ReadMultiPolygonPts(const string &strWsShp)
{
	OGRFile ofile(strWsShp,OGRFile::app);
	OGRLayer* poLayer = ofile.m_pLayer;
	OGRFeature * poFeature;
	
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* pGeometry = poFeature->GetGeometryRef();
		if (pGeometry == NULL)
		{
			return -1;
		}

		OGRwkbGeometryType geoType = pGeometry->getGeometryType();
		if (geoType == wkbMultiPolygon)
		{
			string strWscd = poFeature->GetFieldAsString("WSCD");
			OGRMultiPolygon* poMultiPolygon = (OGRMultiPolygon*)poFeature->GetGeometryRef();
			if (!poMultiPolygon) {
				OGRFeature::DestroyFeature(poFeature);
				continue;
			}
			poMultiPolygon->closeRings();

			vector<int> seperatePoly;
			int cnt(poMultiPolygon->getNumGeometries());


			vector<vector<OGRRawPoint>> unionedPolys;
			for (int i = 0; i < poMultiPolygon->getNumGeometries(); i++)
			{	
				vector< vector<OGRRawPoint> > multiPolygonPts;
				OGRPolygon *poPolygon = (OGRPolygon *)poMultiPolygon->getGeometryRef(i);
				int nInCount = poPolygon->getNumInteriorRings();
				for (int i = 0; i < nInCount; ++i)
				{

					vector<OGRRawPoint> RingPtsT;
					OGRLinearRing* poInRing = poPolygon->getInteriorRing(i);

					RingPtsT.resize(poInRing->getNumPoints());
					poInRing->getPoints((OGRRawPoint*)&*RingPtsT.begin());//循环获得第i个内环所有点到vecintemp[2.22.liang]

					multiPolygonPts.push_back(RingPtsT);
				}
				vector<OGRRawPoint> ExRingPtsT;
				ExRingPtsT.resize(poPolygon->getExteriorRing()->getNumPoints());
				poPolygon->getExteriorRing()->getPoints((OGRRawPoint*)&*ExRingPtsT.begin());

				multiPolygonPts.push_back(ExRingPtsT);
				vector< vector<OGRRawPoint> > uPts;
				UnionWsPts(multiPolygonPts, uPts);
				unionedPolys.push_back(uPts[0]);
			}
			for (;FindPolyToUni(unionedPolys)==0;)
			{
				;
			}
			if (unionedPolys.size()==1)
			{
				OGRPolygon		oPolygon;
				OGRLinearRing	poLinearRing;
				//流域边界点集写入shp--流域空间属性
				poLinearRing.setPoints(unionedPolys.begin()->size(), unionedPolys.begin()->data());

				poLinearRing.closeRings();
				oPolygon.addRing(&poLinearRing);
				poFeature->SetGeometry(&oPolygon);
				poLayer->SetFeature(poFeature);
			}
			else {
				OGRMultiPolygon		oPolygon;
				for (auto &a : unionedPolys)
				{
					OGRPolygon oPoly;
					OGRLinearRing	poLinearRing;
					poLinearRing.setPoints(a.size(), a.data());

					poLinearRing.closeRings();
					oPoly.addRing(&poLinearRing);
					oPolygon.addGeometry(&oPoly);
				}
				poFeature->SetGeometry(&oPolygon);
				poLayer->SetFeature(poFeature);
			}
		}
		else
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}

		OGRFeature::DestroyFeature(poFeature);
	}
	return 0;
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
		if(idx==vwscd.size()) continue;
		vwsgeom[idx] = pp;
		wsarea[idx] = pp->getArea();
	}
}
void st2ws::ws_shp_calc_area(const string & spath)
{
	OGRFile ofilews(spath, OGRFile::in);
	ofilews.m_pLayer->ResetReading();
	//这个边框应包括site cell 
	for (OGRFeature *pFeature; pFeature = ofilews.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		string cd = pFeature->GetFieldAsString("WSCD");
		OGRGeometry *pGeom = pFeature->GetGeometryRef();
		int idx = find(vwscd.begin(), vwscd.end(), cd) - vwscd.begin();
		if (idx == vwscd.size()) continue;
		auto atp = pGeom->getGeometryType();
		if (atp == wkbPolygon || atp == wkbPolygon25D)
		{
			wsarea[idx] = ((OGRPolygon*)pGeom)->get_Area();
		}
		else if (atp == wkbMultiPolygon || atp == wkbMultiPolygon25D)
		{
			wsarea[idx] = ((OGRMultiPolygon*)pGeom)->get_Area();
		}
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

vector<string>& st2ws::wscds()
{
	return vwscd;
}

void st2ws::calc(useable_station u, station_interp s, const string &tmsuffix)
{
	fuseablestation(_usable_station=u);
	fstationinterp(_station_interp=s,tmsuffix);
}
vector<double>& st2ws::ws_rain()
{
	return wsdrp;
}

void st2ws::fuseablestation(useable_station _useableStation)
{
	vector<int> nopst;
	size_t i;
	if (_useableStation == STATICNEAR)
	{
		vector<double> stxtmp, stytmp, stdrptmp;
		for (i = 0; i < stdrp.size(); i++)
		{
			if (stdrp[i] < 0) continue;
			stxtmp.push_back(stx[i]);
			stytmp.push_back(sty[i]);
			stdrptmp.push_back(stdrp[i]);
		}
		for (i = 0; i < stdrp.size(); i++)
		{
			if (stdrp[i] < 0)
			{
				if (!stxtmp.empty()) {
					vector<double> dist;
					InversDistWeight(dist, stx.at(i), sty.at(i), stxtmp, stytmp);
					NeareastDistWeightZ(stdrp.at(i), dist, stdrptmp,nearVal);
				}
				if (stdrp.at(i) < 0) stdrp.at(i) = 0;
			}
		}
	}
	for (i=0;i<stdrp.size();i++)
	{
		if (stdrp[i] < 0)
		{
			if (_useableStation == DYNAMIC) continue;
			else if (_useableStation == STATIC0)  stdrp[i] = 0;
		}
		set_useable_station(i);
	}
}


void st2ws_wsshp(char *s)
{
	string sMark(s);
	sMark += ".st2ws.mark";
	ifstream ifs(sMark);
	if (ifs) return;
	ReadMultiPolygonPts(s);

	OGRFile ofilews(s, OGRFile::app);
	for (OGRFeature *pFeature; pFeature = ofilews.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
	{
		OGRGeometry *pg = pFeature->GetGeometryRef();
		SelfIntersectGeom(pg);
		pFeature->SetGeometry(pg);
		ofilews.m_pLayer->SetFeature(pFeature);
	}

	string sql = string("CREATE SPATIAL INDEX ON ") + ofilews.m_pLayer->GetName();
	ofilews.m_poDS->ExecuteSQL(sql.c_str(), nullptr, nullptr);
	ofstream ofs(sMark);
	if (!ofs) throw errno;
}
