// #include "..\\include\\gdal_priv.h"
// #include "..\\include\\gdal_alg.h"

#include "ogrsf_frmts.h"
//#include "..\\include\\gdal.h"
#include <string>
using namespace std;
class OGRFile;

class GDALReg 
{
	//friend class OGRFile;
public:
	//static int m_nCnt;
	static GEOSContextHandle_t m_geo;

	GDALReg();
	~GDALReg();
};
//表示一个shp文件 excel的一个sheet
class OGRFile//:private GDALReg
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
	void create(OGRSpatialReference *p=nullptr);
	operator bool();
	GEOSContextHandle_t geosctx();
};
