#include <vector>
#include <string>
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
	enum useable_station { DYNAMIC, STATIC0, STATICINVERSE, };
	enum station_interp { KRIGING, THIESSEN, INVERSE, LINER, };
	//添加所有流域测站 调用一次 
	void add_st(double x,double y,const string & cd);
	void add_ws(double x,double y,const string & cd);
	void add_end();
	void ws_shp(const string &spath);
	//循环设置所有测站雨量
	void begin_st_rain();
	void st_rain(double);
	void calc(useable_station u,station_interp s);
	vector<double> &ws_rain();
//private:
	void fuseablestation(useable_station );
	void fstationinterp(station_interp);
	double great_circle_dist(double fi1/*y*/, double lam1, double fi2, double lam2);
	void set_useable_station(int stindex);
	vector<double> wsx, wsy, wsdrp, stx, sty, stdrp,wsarea;
	vector<string> vwscd, vstcd;
	vector<double> stxrun, styrun, stdrprun;
	vector<string> vstcdrun;
	vector<int> stidxrun;
	vector<geos::geom::Geometry*> vwsgeom;
	//使用这些站点时 各个流域的weight
	map<string,vector<map<int, double>>> vweight;
	int useablest,setrainindex;
	useable_station _usable_station;
	station_interp _station_interp;
	
};
ST2WS char *  out_put_voronoi(char *spath, char * wsshp, char * stshp, char ** stcds, int stcnt, char ** wscds, int wscnt);
