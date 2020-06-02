#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/geom/Point.h>
#include <geos/precision/GeometryPrecisionReducer.h>
#include <geos/geom/PrecisionModel.h>
#include<geos/geom/LineString.h>
#include<geos_c.h>
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

#ifdef _DEBUG
#pragma comment(lib,"gdald.lib")
#pragma comment(lib,"libgeosd.lib")
#pragma comment(lib,"geos_cd.lib")
#else 
#pragma comment(lib,"gdal.lib")
#pragma comment(lib,"libgeos.lib")
#pragma comment(lib,"geos_c.lib")
#endif
#pragma comment(lib,"Ws2_32.lib")
#pragma comment(lib,"legacy_stdio_definitions.lib")

/*
https://stackoverflow.com/questions/31473553/is-there-a-way-to-convert-a-self-intersecting-polygon-to-a-multipolygon-in-jts
*/
void OpenRings(vector<OGRRawPoint> &vpts)
{
	double xdiff=vpts.begin()->x-vpts.rbegin()->x,ydiff=vpts.begin()->y-vpts.rbegin()->y;
	if(xdiff==0&&ydiff==0)
	{
		vpts.erase(vpts.begin());
	}
}
int SelfIntersectPt(vector<OGRRawPoint>& vpts)
{
	int needFind = 1, mod = 0;
	while (needFind) {
		needFind = 0;
		for (int p0 = 0, p1 = 1, p2 = 2, p3 = 3; p3 < vpts.size() && 0 == needFind; ++p0, ++p1, ++p2, ++p3)
		{
			OGRLineString ls0, ls1;
			ls0.addPoint(vpts.at(p0).x, vpts.at(p0).y);
			ls0.addPoint(vpts.at(p1).x, vpts.at(p1).y);
			ls1.addPoint(vpts.at(p2).x, vpts.at(p2).y);
			ls1.addPoint(vpts.at(p3).x, vpts.at(p3).y);
			OGRGeometry* pInterGeom = ls0.Intersection(&ls1);
			OGRwkbGeometryType gt;
			OGRPoint* pInterPt = (OGRPoint*)pInterGeom;
			if (pInterGeom && pInterPt
				&& (gt = pInterGeom->getGeometryType()) == wkbPoint)
			{
				mod = 1;
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
			delete pInterGeom;
		}
	}
	return mod;
}
int SelfIntersectGeom(OGRGeometry* pGeom)
{
	int mod;
	if (pGeom->getGeometryType() == wkbPolygon)
	{
		OGRPolygon* pPoly = (OGRPolygon*)pGeom;
		vector<OGRRawPoint> vpts(pPoly->getExteriorRing()->getNumPoints());
		pPoly->getExteriorRing()->getPoints(vpts.data());
		OpenRings(vpts);
		mod = SelfIntersectPt(vpts);
		pPoly->getExteriorRing()->setPoints(vpts.size(), vpts.data());
		pPoly->getExteriorRing()->closeRings();
	}
	else if (pGeom->getGeometryType() == wkbMultiPolygon)
	{
		OGRMultiPolygon* pMp = (OGRMultiPolygon*)pGeom;
		for (int i = 0; i < pMp->getNumGeometries(); ++i)
		{
			OGRPolygon* pPoly = (OGRPolygon*)pMp->getGeometryRef(i);
			vector<OGRRawPoint> vpts(pPoly->getExteriorRing()->getNumPoints());
			pPoly->getExteriorRing()->getPoints(vpts.data());
			OpenRings(vpts);
			mod = SelfIntersectPt(vpts);
			pPoly->getExteriorRing()->setPoints(vpts.size(), vpts.data());
			pPoly->getExteriorRing()->closeRings();
		}
	}
	return mod;
}

OGRGeometry* ReducePrecision(const OGRGeometry *oPoly,double dPrecision)
{
	GEOSGeom hGeosGeom = NULL;
    OGRGeometry *poOGRProduct = NULL;

    GEOSContextHandle_t hGEOSCtxt = OGRGeometry::createGEOSContext();
    hGeosGeom = oPoly->exportToGEOS(hGEOSCtxt);
    if( hGeosGeom != NULL )
    {
		typedef geos::geom::GeometryFactory GeometryFactory;
		geos::geom::PrecisionModel pm_fixed_(dPrecision);	
		auto factory_fixed_(GeometryFactory::create(&pm_fixed_, 0));
		geos::precision::GeometryPrecisionReducer reducer_(*factory_fixed_);
		std::unique_ptr<geos::geom::Geometry> hGeosProduct= reducer_.reduce(*(geos::geom::Geometry*)hGeosGeom);

        GEOSGeom_destroy_r( hGEOSCtxt, hGeosGeom );
		
        if( hGeosProduct.get())
        {
            poOGRProduct = OGRGeometryFactory::createFromGEOS(hGEOSCtxt, (GEOSGeom)(hGeosProduct.get()));
            if( poOGRProduct != NULL && oPoly->getSpatialReference() != NULL )
                poOGRProduct->assignSpatialReference(oPoly->getSpatialReference());
            //poOGRProduct = OGRGeometry::OGRGeometryRebuildCurves(this, NULL, poOGRProduct);
            GEOSGeom_destroy_r( hGEOSCtxt, (GEOSGeom)hGeosProduct.get() );
			hGeosProduct.release();
        }
    }
    OGRGeometry::freeGEOSContext(hGEOSCtxt);
    return poOGRProduct;
}
OGRGeometry* GeomIntersection(OGRGeometry *lPoly,OGRGeometry *rPoly,double dPrecision)
{
	OGRGeometry* pInterGeom =  lPoly->Intersection(rPoly);
	if(!pInterGeom)
	{
		//cout<<"intersection failed.try reduce precision."<<endl;
		OGRGeometry *l =ReducePrecision(lPoly,dPrecision),*r = ReducePrecision(rPoly,dPrecision);
		pInterGeom=l->Intersection(r);
		delete l;
		delete r;
		if(!pInterGeom)
		{
			//cout<<"intersection with reduced precision failed.try buffer zero."<<endl;
			OGRGeometry *p1= lPoly->Buffer(0),*p2=rPoly->Buffer(0);//缓冲
			pInterGeom=p1->Intersection(p2);
			delete p1;
			delete p2;
		}
	}
	if (!pInterGeom)
	{
		cout << "geominter" << endl;
	}
	return pInterGeom;
}
int main(int argc,char*argv[])
{
	const string& stFile(argv[1]);
	const string& ofilename(argv[2]);
	vector<map<int, double>> weight;
	vector<string> vstcdrun, cellcd2;
	vector<double> cellarea;
	vector<int> sitePoly2;	//sitePoly2 : index - site index ;val - cell index
	vector<OGRGeometry*> ogeos;
	double* p1, * p2;
	int idxrun, sitecnt, cellcnt;

	GDALReg r;
	{
		OGRFile ofile(ofilename, OGRFile::app);
		cellcnt = ofile.m_pLayer->GetFeatureCount();
		for (OGRFeature* pFeature; pFeature = ofile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
		{
			OGRGeometry* pg = pFeature->GetGeometryRef();
			if (SelfIntersectGeom(pg)) {
				pFeature->SetGeometry(pg);
				ofile.m_pLayer->SetFeature(pFeature);
			}
			cellcd2.push_back(pFeature->GetFieldAsString("WSCD"));
			OGRGeometry* pGeom = pFeature->GetGeometryRef();
			auto atp = pGeom->getGeometryType();
			if (atp == wkbPolygon || atp == wkbPolygon25D)
			{
				cellarea.push_back(((OGRPolygon*)pGeom)->get_Area());
			}
			else if (atp == wkbMultiPolygon || atp == wkbMultiPolygon25D)
			{
				cellarea.push_back(((OGRMultiPolygon*)pGeom)->get_Area());
			}
		}
		string sql = string("CREATE SPATIAL INDEX ON ") + ofile.m_pLayer->GetName();
		ofile.m_poDS->ExecuteSQL(sql.c_str(), nullptr, nullptr);
	}
	OGRFile ofile(ofilename);
	geos::geom::CoordinateArraySequence seq;
	geos::geom::Envelope env;
	geos::triangulate::VoronoiDiagramBuilder builder;
	const geos::geom::GeometryFactory& geomFact(*geos::geom::GeometryFactory::getDefaultInstance());
	std::unique_ptr<geos::geom::GeometryCollection> polys;	
	vector<geos::geom::Coordinate> cords;
	{
		OGRFile owsfile(stFile);
		for (OGRFeature* pFeature; pFeature = owsfile.m_pLayer->GetNextFeature(); OGRFeature::DestroyFeature(pFeature))
		{
			OGRPoint* oPt = (OGRPoint*)pFeature->GetGeometryRef();
			geos::geom::Coordinate coordrun(oPt->getX(), oPt->getY());
			if (find_if(cords.begin(), cords.end(), [&](geos::geom::Coordinate a)->bool {
				return a.equals2D(coordrun);
				}) != cords.end())
			{
				continue;
			}
			string cd = pFeature->GetFieldAsString("STCD");
			vstcdrun.push_back(cd);
			seq.add(coordrun);
			env.expandToInclude(coordrun);
			cords.push_back(coordrun);
		}
		sitePoly2.resize(vstcdrun.size());
		sitecnt = vstcdrun.size();
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
	{
		for (idxrun = 0;idxrun!=cords.size(); idxrun++)
		{
			auto oPt = cords.at(idxrun);
			for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
			{
				const geos::geom::Point *pcoord = (const geos::geom::Point *)(polys->getGeometryN(i)->getUserData());
				if (oPt.x == pcoord->getX() && oPt.y == pcoord->getY())
				{
					sitePoly2[idxrun] = i;
					break;
				}
			}
		}
		sitePoly2.resize(vstcdrun.size());
	}
	for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
	{
		auto tg = OGRGeometryFactory::createFromGEOS(ofile.geosctx(), (GEOSGeom)polys->getGeometryN(i));
		SelfIntersectGeom(tg);
		ogeos.push_back(tg);
	}
	vector<map<int, OGRRawPoint>> weight2;
	weight2.resize(ofile.m_pLayer->GetFeatureCount());
	weight.resize(ofile.m_pLayer->GetFeatureCount());
	CPLSetConfigOption("CPL_LOG", "l.txt");
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
				continue;//
			}
			string cd = poFea->GetFieldAsString("WSCD");
			vector<string>::iterator iter;
			if ((iter = find(cellcd2.begin(), cellcd2.end(), cd)) == cellcd2.end())
			{
				continue;
			}
			int cellidx = iter - cellcd2.begin();
			OGRPoint oCenterPt;
			poFea->GetGeometryRef()->Centroid(&oCenterPt);
			if (poFea->GetGeometryRef()->Within(stGeom) )
			{
				#pragma omp critical
				{
					weight[cellidx][idxrun] = cellarea[ cellidx];
					weight2[cellidx][idxrun].x = oCenterPt.getX();
					weight2[cellidx][idxrun].y = oCenterPt.getY();
				}
			}
			else if (poFea->GetGeometryRef()->Contains(stGeom))
			{
				#pragma omp critical
				{
					weight[cellidx][idxrun] = ((OGRPolygon *)stGeom)->get_Area();
					weight2[cellidx][idxrun].x = oCenterPt.getX();
					weight2[cellidx][idxrun].y = oCenterPt.getY();
				}
			}
			else if (poFea->GetGeometryRef()->Intersects(stGeom))
			{
				OGRGeometry *pGeom =GeomIntersection(poFea->GetGeometryRef(),stGeom,1e5);
				if (pGeom->getGeometryType() == wkbMultiPolygon)
				{
					#pragma omp critical
					{
						weight[cellidx][idxrun] = ((OGRMultiPolygon *)pGeom)->get_Area();
						OGRPoint opt;
						((OGRMultiPolygon*)pGeom)->Centroid(&opt);
						weight2[cellidx][idxrun].x = opt.getX();
						weight2[cellidx][idxrun].y = opt.getY();
					}
				}
				else if (pGeom->getGeometryType() == wkbPolygon)
				{
					#pragma omp critical
					{
						weight[cellidx][idxrun] = ((OGRPolygon *)pGeom)->get_Area();
						OGRPoint opt;
						((OGRPolygon*)pGeom)->Centroid(&opt);
						weight2[cellidx][idxrun].x = opt.getX();
						weight2[cellidx][idxrun].y = opt.getY();
					}
				}
				else if (pGeom->getGeometryType() == wkbGeometryCollection)
				{
					#pragma omp critical
					{
						weight[cellidx][idxrun] = ((OGRGeometryCollection *)pGeom)->get_Area();
						OGRPoint opt;
						((OGRGeometryCollection*)pGeom)->Centroid(&opt);
						weight2[cellidx][idxrun].x = opt.getX();
						weight2[cellidx][idxrun].y = opt.getY();
					}
				}
				else
				{
#pragma omp critical
					{
						cout << "unexpected type:"<< pGeom->getGeometryType() << endl;
						throw 1;
					}
				}
				delete pGeom;
			}
		}
		ofilea.m_pLayer->SetSpatialFilter(nullptr);
	}
	CPLSetConfigOption("CPL_LOG", nullptr);
	for (int i = 0; i < cellcnt; i++)
	{
		double sum(0);
		for (auto &a : weight[i])
		{
			sum += a.second;
		}
		double dtest(fabs(sum - cellarea[i])/ cellarea[i]);
		if(dtest>0.05) std::cout<<"??"<<i<<std::endl;
	}

	OGRFile ovrfile("Voronoi.shp", OGRFile::out, "", "ESRI Shapefile", wkbPolygon);
	OGRFile ostfile("VoronoiPt.shp", OGRFile::out, "", "ESRI Shapefile", wkbPoint);
	OGRFieldDefn ostcd("STCD", OFTString), ostval("VALUE", OFTReal),owscd("WSCD",OFTString);
	ostfile.m_pLayer->CreateField(&ostcd);
	ostfile.m_pLayer->CreateField(&ostval);
	ostfile.m_pLayer->CreateField(&owscd);
	ovrfile.m_pLayer->CreateField(&ostcd);
	
	for (idxrun = 0; idxrun < sitecnt; idxrun++)
	{
		OGRGeometry* stGeom = ogeos[sitePoly2[idxrun]];
		OGRFeature oFeature(ovrfile.m_pLayer->GetLayerDefn());
		oFeature.SetField(ostcd.GetNameRef(), vstcdrun.at(idxrun).c_str());
		oFeature.SetGeometry(stGeom);
		ovrfile.m_pLayer->CreateFeature(&oFeature);

		int cellidx = 0;
		for (auto& aCellWeight : weight)
		{
			auto iterFind = aCellWeight.find(idxrun);
			if (iterFind != aCellWeight.end())
			{
				OGRFeature oFeaturePt(ostfile.m_pLayer->GetLayerDefn());
				oFeaturePt.SetField(ostcd.GetNameRef(), vstcdrun.at(idxrun).c_str());
				oFeaturePt.SetField(owscd.GetNameRef(), cellcd2.at(cellidx).c_str());
				oFeaturePt.SetField(ostval.GetNameRef(), iterFind->second/cellarea[cellidx]);
				
				OGRPoint opt(weight2[cellidx].at(idxrun).x, weight2[cellidx].at(idxrun).y);
				oFeaturePt.SetGeometry(&opt);
				ostfile.m_pLayer->CreateFeature(&oFeaturePt);
			}
			cellidx++;
		}
	}
	ovrfile.close();
	ostfile.close();
	for(auto &a:ogeos)
	{
		OGRGeometryFactory::destroyGeometry(a);
	}
	ogeos.clear();
	for (std::size_t i = 0; i < polys->getNumGeometries(); ++i)
	{
		geos::geom::Point *pcoord = (geos::geom::Point *)polys->getGeometryN(i)->getUserData();
		geomFact.destroyGeometry(pcoord);
	}
}