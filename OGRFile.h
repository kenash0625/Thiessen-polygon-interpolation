// #include "..\\include\\gdal_priv.h"
// #include "..\\include\\gdal_alg.h"

#include "..\\include\\ogrsf_frmts.h"
//#include "..\\include\\gdal.h"
#include <set>
#include <string>
#include <map>
#include <vector>
using namespace std;
class OGRFile;

class GDALReg 
{
	friend class OGRFile;
	static int m_nCnt;
	GDALReg();
	~GDALReg();
};
//表示一个shp文件 excel的一个sheet
class OGRFile:private GDALReg
{
	OGRFile(const OGRFile &rhs);
	OGRFile operator=(const OGRFile &rhs);
public:
	enum openmode{in=0,app,out};
	openmode m_op;
	string m_name,m_drvname,m_layername;
	OGRwkbGeometryType m_gt;
	GDALDriver *m_pDrv;
	GDALDataset *m_poDS;
	OGRSpatialReference *m_psrs;
	OGRLayer *m_pLayer;
public:	
	OGRFile(string strName="", openmode mode=in, string strLayerName="",string drvname="",OGRwkbGeometryType geotype=wkbUnknown,OGRSpatialReference *psrs=nullptr);
	void init(string strName="", openmode mode=in, string strLayerName="",string drvname="",OGRwkbGeometryType geotype=wkbUnknown,OGRSpatialReference *psrs=nullptr);
	~OGRFile();
	void close();
	void open(openmode op=in);
	void erase();
	void create(OGRSpatialReference *p=nullptr);
	int copy(string strname);
	operator bool();

	//vector 改变大小时会调用元素的析构函数 有poSRS时会delete两次poSRS-OGRFeature::DsestroyFeature

	static int equals(const OGRPoint &lhs,const OGRPoint &rhs);
	static void makelonger(OGRLineString &olslonger,const OGRLineString &ols,double len=1500);
	
	static int insertatlen(vector<OGRPoint> &vlspt,double dlen,OGRPoint *pInsertPt=nullptr);
	
	static int putinterptonline( OGRLineString &ols,const OGRLineString &olsr,const OGRGeometry *pGeomInter,int nGeomIndex=-1,OGRPoint *pInPt=nullptr );
	static int putinterptonline( vector<OGRPoint> &nIndex,OGRLineString &ols,const OGRLineString &olsr );
	static int indexofptonline(const OGRPoint &opt,const OGRLineString &ols);
	//给两个点的线段插入一些点
	static int insertatrate(vector<OGRPoint> &vlspt,map<double,map<int,OGRPoint>> &drate,bool bcalc=false);
	static int insertatrate(vector<OGRPoint> &vptinsert,const vector<OGRPoint> &vlspt,map<double,map<int,OGRPoint>> &drate);
	static int insertatrate(vector<OGRPoint> &vlspt,double drate,OGRPoint *pInsertPt=nullptr);

	static void line2pt(vector<OGRPoint> &vpt,const OGRLineString &ols);
	static void pt2line(OGRLineString &ols,const vector<OGRPoint> &vpt);
	bool addfld(OGRFieldDefn &ofld);

	int spatialfilter_plpt(vector<int> &vfids, const OGRPoint &opt);
	static void makelonger2(OGRLineString &olslonger, const OGRLineString &ols, double len=1500);
	//ogrpoint合成的复制构造函数 不会增加投影引用计数 析构时会减少并delete多次
	//spref 打开文件时为1
	static void clonepoint(OGRPoint &odst, const OGRPoint &osrc);
};
struct MyPoint{
	double x;
	double y;
	bool operator<(const MyPoint &rhs) const
	{
		if(x < rhs.x) return true;
		else if (x == rhs.x) return y < rhs.y;
		else return false;
	}
	bool operator==(const MyPoint &rhs) const
	{
		return (fabs(x-rhs.x)<0.001 && fabs(y-rhs.y)<0.001);
	}
};
