
//#include "Files.h"

#include "OGRFile.h"
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;

#pragma comment(lib,"..\\lib\\gdal_i.lib")
int GDALReg::m_nCnt=0;
//char* pszOldEnc, *pszOldUTF8;
GDALReg::GDALReg()
{	
	m_nCnt++;
	if(m_nCnt>1) return;
	GDALAllRegister();
	OGRRegisterAll();
	// backup old value
	//const char* pszOldValTmp = CPLGetConfigOption("SHAPE_ENCODING", NULL);
	//pszOldEnc = pszOldValTmp ? CPLStrdup(pszOldValTmp) : NULL;
	//const char* pszOldValTmp2 = CPLGetConfigOption("GDAL_FILENAME_IS_UTF8", NULL);
	//pszOldUTF8 = pszOldValTmp2 ? CPLStrdup(pszOldValTmp2) : NULL;

	
 	CPLSetConfigOption("SHAPE_ENCODING", "");
 	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
// 	CPLSetConfigOption("OGR_XLSX_HEADERS","FORCE");
// override with new value
//	CPLSetConfigOption(pszKey, pszNewVal);
	// do something useful
	// restore old value

}

GDALReg::~GDALReg()
{
	m_nCnt--;
	if (m_nCnt==0)
	{
		//CPLSetConfigOption("SHAPE_ENCODING", pszOldEnc);
		//CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", pszOldUTF8);
		//CPLFree(pszOldUTF8);
		//CPLFree(pszOldEnc);
		OGRCleanupAll();
	}
}

OGRFile::OGRFile( string strName, openmode mode, string strLayerName,string drvname,OGRwkbGeometryType geotype, OGRSpatialReference *psrs )
{
	init(strName,mode,strLayerName,drvname,geotype,psrs);
}

OGRFile::operator bool()
{
	return m_poDS!=nullptr && m_pLayer!=nullptr;
}

OGRFile::~OGRFile()
{
	close();
}

void OGRFile::close()
{
	if(m_poDS)
	{
		GDALClose(m_poDS);

	}	
	m_poDS=nullptr;
	m_pLayer=nullptr;
}

void OGRFile::open(openmode op)
{
	close();
	m_op=op;

	if((m_poDS=(GDALDataset *)GDALOpenEx(m_name.c_str(),(op==OGRFile::app?GDAL_OF_UPDATE:GDAL_OF_READONLY)||GDAL_OF_VECTOR,nullptr,nullptr,nullptr)) && (m_pDrv=m_poDS->GetDriver())) m_pLayer=m_poDS->GetLayer(0);
}

void OGRFile::erase()
{
	close();
	if(m_pDrv) m_pDrv->QuietDelete(m_name.c_str());
}

int OGRFile::equals(const OGRPoint &lhs,const OGRPoint &rhs)//geos没有这种选项可以设?
{
	return (fabs(lhs.getX()-rhs.getX())<1e-3 && fabs(lhs.getY()-rhs.getY())<1e-3);
}

// int OGRFile::copy( string strname )
// {
// 	open();
// 	if (*this && m_pDrv)
// 	{
// 		m_pDrv->DeleteDataSource(strname.c_str());
// 		OGRDataSource *pnewds= m_pDrv->CopyDataSource(m_poDS,strname.c_str());
// 		OGRDataSource::DestroyDataSource(pnewds);
// 		OGRFile ofnew(strname,OGRFile::in);
// 		if(ofnew) return 0;
// 	}
// 	return 1;
// }

void OGRFile::init( string strName, openmode mode/*=in*/, string strLayerName/*=""*/,string drvname/*=""*/,OGRwkbGeometryType geotype/*=wkbUnknown*/, OGRSpatialReference *psrs )
{
	m_poDS=nullptr;
	m_pLayer=nullptr;
	m_pDrv=nullptr;
	m_name=strName;
	m_drvname=drvname;
	m_layername=strLayerName;
	m_op=mode;
	m_gt=geotype;
	m_psrs = psrs;
	if(strName.empty()) return;
	if (mode==out)
	{
		create(m_psrs);
	} 
	else
	{
		open(mode);
	}		
}

void OGRFile::create(OGRSpatialReference *p)
{
	close();
	m_pDrv = GetGDALDriverManager()->GetDriverByName(m_drvname.c_str());
	if (m_pDrv != nullptr)
	{
		m_pDrv->QuietDelete(m_name.c_str());
		m_poDS = m_pDrv->Create(m_name.c_str(), 0, 0, 0, GDT_Unknown, NULL);
		if(m_poDS)	m_pLayer = m_poDS->CreateLayer("out", p, m_gt, NULL);
	}
}
bool OGRFile::addfld(OGRFieldDefn &ofld)
{
	bool bret=false;
	if (*this)
	{
		if(-1==this->m_pLayer->GetLayerDefn()->GetFieldIndex(ofld.GetNameRef())) 
			bret=(OGRERR_NONE == m_pLayer->CreateField(&ofld));
		else bret=true;
	}	
	return bret;
}
//等价python的positionalongline  至少有两个点 dlen比vlspt的长度短
int OGRFile::insertatlen( vector<OGRPoint> &vlspt,double dlen,OGRPoint *pInsertPt )
{
	OGRLineString ols;
	ols.addPoint(&vlspt.at(0));
	ols.addPoint(&vlspt.at(1));

	int nindex;
	double dlenrun;	
	OGRPoint opt;
	if (ols.get_Length()>dlen)
	{
		ols.Value(dlen,&opt);		
		nindex=1;
	}
	else
	{
		for (nindex=1;(dlenrun=ols.get_Length())<=dlen;) 
		{
			/*if(++nindex>=vlspt.size()) break;*/
			ols.addPoint(&vlspt.at(++nindex));
		}
		OGRLineString olstmp;
		olstmp.addPoint(&vlspt.at(nindex-1));
		olstmp.addPoint(&vlspt.at(nindex));
		double ddist=dlen-(dlenrun-olstmp.get_Length());
		olstmp.Value(ddist,&opt);
	}
	vlspt.insert(vlspt.begin()+nindex,opt);	
	if(pInsertPt!=nullptr) *pInsertPt=opt;
	return nindex;
}
void OGRFile::makelonger2( OGRLineString &olslonger,const OGRLineString &ols,double len/*=1500*/ )
{
	OGRPoint o1,o2;
	ols.StartPoint(&o1);
	ols.EndPoint(&o2);

	double a1=o1.getX(),b1=o1.getY(),a2=o2.getX(),b2=o2.getY(),k=(b1-b2)/(a1-a2),x1,x2,y1,y2;
	a1=(a1+a2)/2,b1=(b1+b2)/2;
	if (fabs(a1-a2)<1e-6)
	{
		x1=a1,x2=a1,y1=b1-len,y2=b1+len;
	}
	else
	{
		x1=a1+len/sqrt(1+k*k),y1=k*(x1-a1)+b1,x2=a1-len/sqrt(1+k*k),y2=k*(x2-a1)+b1;
	}

	olslonger.empty();
	olslonger.addPoint(x1,y1);
	olslonger.addPoint(x2,y2);
}
void OGRFile::makelonger( OGRLineString &olslonger,const OGRLineString &ols,double len/*=1500*/ )
{
	len=ols.get_Length();
	OGRPoint o1,o2;
	ols.StartPoint(&o1);
	ols.EndPoint(&o2);

	double a1=o1.getX(),b1=o1.getY(),a2=o2.getX(),b2=o2.getY(),k=(b1-b2)/(a1-a2),x1,x2,y1,y2;
	a1=(a1+a2)/2,b1=(b1+b2)/2;
	if (fabs(a1-a2)<1e-6)
	{
		x1=a1,x2=a1,y1=b1-len,y2=b1+len;
	}
	else
	{
		x1=a1+len/sqrt(1+k*k),y1=k*(x1-a1)+b1,x2=a1-len/sqrt(1+k*k),y2=k*(x2-a1)+b1;
	}

	olslonger.empty();
	olslonger.addPoint(x1,y1);
	olslonger.addPoint(x2,y2);
}

void OGRFile::line2pt( vector<OGRPoint> &vpt,const OGRLineString &ols )
{
	vpt.resize(ols.getNumPoints());
	for(int i=0;i!=ols.getNumPoints();i++)
	{
		ols.getPoint(i,&vpt[i]);
	}
}

void OGRFile::pt2line( OGRLineString &ols,const vector<OGRPoint> &vpt )
{
	ols.empty();
	for_each(vpt.begin(),vpt.end(),[&](OGRPoint pt)->void {ols.addPoint(&pt);});
}
//ols-圈 olsr-两个点的长线 pPt-参数所指交点 交点在圈上存在或者需要插入？
int OGRFile::putinterptonline( OGRLineString &ols,const OGRLineString &olsr,const OGRGeometry *pGeomInter,int nGeomIndex,OGRPoint *pInPt )
{
	OGRPoint *pPt=nullptr;
	OGRwkbGeometryType gt=pGeomInter->getGeometryType();
	if (gt==wkbPoint || gt==wkbPoint25D)
	{
		pPt=(OGRPoint *)pGeomInter;
	}
	else if (gt==wkbMultiPoint || gt==wkbMultiPoint25D)
	{
		pPt=(OGRPoint *)(((OGRMultiPoint *)pGeomInter)->getGeometryRef(nGeomIndex));
	}
	if(pPt==nullptr) return -1;

	vector<OGRPoint> vpttmp;
	line2pt(vpttmp,ols);
	for (int nIndes=0;nIndes!=vpttmp.size()-1;nIndes++)
	{
		OGRLineString olstmp;//圈上两点
		olstmp.addPoint(&vpttmp.at(nIndes));
		olstmp.addPoint(&vpttmp.at(nIndes+1));
		if (fabs(pPt->getX()-vpttmp.at(nIndes).getX())<1e-3 && fabs(pPt->getY()-vpttmp.at(nIndes).getY())<1e-3)
		{
			if(pInPt) *pInPt=vpttmp.at(nIndes);
			return nIndes;
		}
		else if (fabs(pPt->getX()-vpttmp.at(nIndes+1).getX())<1e-3 && fabs(pPt->getY()-vpttmp.at(nIndes+1).getY())<1e-3)
		{
			if(pInPt) *pInPt=vpttmp.at(nIndes);
			return nIndes+1;
		}
		else if (olsr.Intersects(&olstmp) )
		{
			OGRPoint *ptmp=(OGRPoint *)olsr.Intersection(&olstmp);
			if (fabs(pPt->getX()-ptmp->getX())<1e-3 && fabs(pPt->getY()-ptmp->getY())<1e-3)
			{
				vpttmp.insert(vpttmp.begin()+nIndes+1,*pPt);
				pt2line(ols,vpttmp);
				if(pInPt) *pInPt=*pPt;
				return nIndes+1;
			}
		}
	}
	return -1;
}

int OGRFile::spatialfilter_plpt(vector<int> &vfids, const OGRPoint &opt)
{
	if (!*this)
	{
		return 1;
	}
	m_pLayer->SetSpatialFilter(nullptr);
	m_pLayer->SetSpatialFilter((OGRPoint *)&opt);
	
	long *pints=nullptr;/*m_pLayer->GetFeatureCountSF(0);*/
	vector<int> vfidsret;
	while (pints && *pints!=-1)
	{
		vfidsret.push_back(*pints);
		pints++;
	}
	for(vector<int>::iterator iter=vfidsret.begin();iter!=vfidsret.end();iter++)
	{
		OGRFeature *pfea=m_pLayer->GetFeature(*iter);
		OGRLineString *pls=((OGRPolygon *)pfea->GetGeometryRef())->getExteriorRing();
		vector<OGRRawPoint> opttmp(pls->getNumPoints());
		pls->getPoints(&*opttmp.begin());
		if (find_if(opttmp.begin(),opttmp.end(),[&](OGRRawPoint orpt)->bool {return fabs(orpt.x-opt.getX())<1e-3 && fabs(orpt.y-opt.getY())<1e-3;})!=opttmp.end())
		{
			vfids.push_back(*iter);
		}
		OGRFeature::DestroyFeature(pfea);
	}
	return vfids.size();
}

int OGRFile::putinterptonline( vector<OGRPoint> &nIndex,OGRLineString &ols,const OGRLineString &olsr )
{
	OGRMultiPoint *pmp=(OGRMultiPoint *)ols.Intersection(&olsr);
	if(pmp==nullptr) return -1;

	int icnt=pmp->getNumGeometries();
	if (icnt < 0) return -1;
	//取出所有交点
	//vector<OGRPoint> vpttmp,vInterPts;
	for (int i=0;i<icnt;i++)
	{
		OGRPoint opt;
		int ntmp=putinterptonline(ols,olsr,pmp,i,&opt);
		if(ntmp==-1) return -1;
		nIndex.push_back(opt);
	}

	return nIndex.size();
}

int OGRFile::insertatrate( vector<OGRPoint> &vlspt,double drate,OGRPoint *pInsertPt/*=nullptr*/ )
{
	OGRLineString ols;
	pt2line(ols,vlspt);
	return insertatlen(vlspt,drate*ols.get_Length(),pInsertPt);
}

int OGRFile::insertatrate( vector<OGRPoint> &vlspt,map<double,map<int,OGRPoint>> &drate,bool bcalc )
{
	OGRLineString ols;
	pt2line(ols,vlspt);
	const double dlen=ols.get_Length();

	vector<OGRPoint>::iterator nPtIndex=vlspt.begin();
	double dLenRun=0;
	for (map<double,map<int,OGRPoint>>::iterator itRate=drate.begin();itRate!=drate.end();itRate++)
	{
		double dThisLen=itRate->first*dlen;
		for (;dLenRun<dThisLen;nPtIndex++)
		{
			dLenRun+=nPtIndex->Distance(&*(nPtIndex+1));
		}
		OGRLineString olstmp;
		OGRPoint opttmp, &optl=*(nPtIndex-1),&optr=*nPtIndex;

		olstmp.addPoint(&optl);
		olstmp.addPoint(&optr);
		olstmp.Value(olstmp.get_Length()-(dLenRun-dThisLen),&opttmp);
		if (bcalc)
		{
			double D=optl.Distance(&optr),d=optl.Distance(&opttmp),h1=optl.getZ(),h2=optr.getZ();
			opttmp.setZ(d/D*(h2-h1)+h1);
		}
		nPtIndex=vlspt.insert(nPtIndex,opttmp);
		dLenRun=dThisLen;
	}
	return 0;
}

int OGRFile::insertatrate(vector<OGRPoint> &vptinsert, const vector<OGRPoint> &vlsptr, map<double, map<int, OGRPoint>> &drate)
{
	OGRLineString ols;
	pt2line(ols,vlsptr);
	const double dlen=ols.get_Length();

	vector<OGRPoint> vlspt(vlsptr.begin(), vlsptr.end());
	vector<OGRPoint>::iterator nPtIndex=vlspt.begin();
	vptinsert.clear();
	vptinsert.push_back(*vlsptr.cbegin());
	double dLenRun=0;
	for (map<double,map<int,OGRPoint>>::iterator itRate=drate.begin();itRate!=drate.end();itRate++)
	{
		double dThisLen=itRate->first*dlen;
		for (;dLenRun<dThisLen;nPtIndex++)
		{
			dLenRun+=nPtIndex->Distance(&*(nPtIndex+1));
		}
		OGRLineString olstmp;
		OGRPoint opttmp, &optl=*(nPtIndex-1),&optr=*nPtIndex;

		olstmp.addPoint(&optl);
		olstmp.addPoint(&optr);
		olstmp.Value(olstmp.get_Length()-(dLenRun-dThisLen),&opttmp);

		double D=optl.Distance(&optr),d=optl.Distance(&opttmp),h1=optl.getZ(),h2=optr.getZ();
		opttmp.setZ(d/D*(h2-h1)+h1);

		nPtIndex=vlspt.insert(nPtIndex,opttmp);
		vptinsert.push_back(opttmp);
		dLenRun=dThisLen;
	}
	vptinsert.push_back(*vlsptr.crbegin());
	return 0;
}

int OGRFile::indexofptonline( const OGRPoint &opt,const OGRLineString &ols )
{
	for (int i=0;i<ols.getNumPoints();i++)
	{
		OGRPoint opttmp;
		ols.getPoint(i,&opttmp);
		if (fabs(opttmp.getX()-opt.getX())<1e-3 && fabs(opttmp.getY()-opt.getY())<1e-3)
		{
			return i;
		}
	}
	return -1;
}

void OGRFile::clonepoint(OGRPoint &odst, const OGRPoint &osrc)
{
	odst.setX(osrc.getX());
	odst.setY(osrc.getY());
	odst.setZ(osrc.getZ());
	odst.assignSpatialReference(osrc.getSpatialReference());
}

