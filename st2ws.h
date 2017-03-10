#include <vector>
#include <string>
#include <map>
using namespace std;
namespace geos
{
	namespace geom
	{
		class Geometry;
	}
}

struct  st2ws
{	
	enum useable_station { DYNAMIC, STATIC0, STATICINVERSE, };
	enum station_interp { KRIGING, THIESSEN, INVERSE, LINER, };
	//添加所有流域测站 调用一次 
	void add_st(double x,double y,const string & cd);
	void add_ws(geos::geom::Geometry*,const string & cd);
	void add_end();
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
