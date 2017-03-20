
//#include "Files.h"

#include "OGRFile.h"
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;
#ifdef GDALDEBUG
#pragma comment(lib,"../debug/gdal_i.lib")
#else 
#pragma comment(lib,"../release/gdal_i.lib")
#endif
int GDALReg::m_nCnt=0;

GEOSContextHandle_t GDALReg::m_geo=nullptr;

//char* pszOldEnc, *pszOldUTF8;
GDALReg::GDALReg()
{	
	m_nCnt++;
	if(m_nCnt>1) return;
	m_geo=OGRGeometry::createGEOSContext();
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
		OGRGeometry::freeGEOSContext(m_geo);
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

GEOSContextHandle_t OGRFile::geosctx()
{
	return m_geo;
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
