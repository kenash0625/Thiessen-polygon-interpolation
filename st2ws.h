#include <vector>
#include <string>
#include <map>
#ifdef ST2WSEXPORT 
#define ST2WS __declspec(dllexport)
#else 
#define ST2WS __declspec(dllimport)
#endif
using namespace std;
namespace geos
{
	namespace geom
	{
		class Geometry;
	}
}

struct ST2WS st2ws
{	
	enum useable_station { DYNAMIC, STATIC0, STATICNEAR, };
	enum station_interp { KRIGING, THIESSEN, INVERSE, LINER,THIESSEN2 };

	void add_st(double x,double y,const string & cd);
	void add_ws(double x,double y,const string & cd);
	void add_ws(double x, double y, double area, const string & cd);
	void add_end();
	void ws_shp(const string &spath);
	void ws_shp_calc_area(const string &spath);

	void begin_st_rain();
	void st_rain(double);
	vector<string> &wscds();
	void calc(useable_station u,station_interp s,const string &tmsuffix="");
	vector<double> &ws_rain();
	void add_st_and_rain(map<string,vector<double>> &strain,int rainindex);
//private:
	void fuseablestation(useable_station );
	void fstationinterp(station_interp,const string &);
	double great_circle_dist(double fi1/*y*/, double lam1, double fi2, double lam2);
	void set_useable_station(int stindex);
	vector<double> wsx, wsy, wsdrp, stx, sty, stdrp,wsarea;
	vector<string> vwscd, vstcd;
	vector<double> stxrun, styrun, stdrprun;
	vector<string> vstcdrun;
	vector<int> stidxrun;
	vector<geos::geom::Geometry*> vwsgeom;

	map<string,vector<map<int, double>>> vweight;
	int useablest,setrainindex;
	useable_station _usable_station;
	station_interp _station_interp;
	string wsfile;
	double nearVal;
	
};

extern "C" {
//一个shp 处理一次  去掉自交多边形 和 合并能合并的 多多边形 建 空间索引
ST2WS void st2ws_wsshp(char *s);
}