#include <vector>
#include <string>
#include <map>
#include <list>
#include <mutex>
#include <condition_variable>
using namespace std;
namespace geos
{
	namespace geom
	{
		class Geometry;
	}
}
struct thi_weight_area
{
	int useablest;
	vector<int> stidx;
	vector<int> voropoly;
	vector<OGRRawPoint> voropts;
	vector<map<int, double>> weights;
};
struct  st2ws
{	
	enum useable_station { DYNAMIC  };
	enum station_interp {  THIESSEN  };
	//添加所有流域测站 调用一次 
	void add_st(double x,double y,const string & cd);
	void add_ws(geos::geom::Geometry*,const string & cd);
	void add_end();
	//循环设置所有测站雨量
	void begin_st_rain();
	void st_rain(double);
	void calc(useable_station u,station_interp s);
	vector<double> &ws_rain();
	
	vector<int> stidxrun;
	int useablest;
	thi_weight_area *weightrun;
	vector<geos::geom::Geometry*> vwsgeom;
	void rand_run();
private:
	void fuseablestation(useable_station );
	void fstationinterp(station_interp);
	void thiessen_calcz(int sitecnt, double *sitex, double *sitey, double *sitez, int cellcnt, geos::geom::Geometry **cells, double*cellarea, double *cellz, thi_weight_area *pweight);
	void set_useable_station(int stindex);
	thi_weight_area* find_weight_area();
public:
	vector<double>  wsdrp, stx, sty, stdrp,wsarea;
	vector<string> vwscd, vstcd;
	vector<double> stxrun, styrun, stdrprun;
	vector<string> vstcdrun;


	//使用这些站点时 各个流域的weight
	list<thi_weight_area> tho;
	int setrainindex;
	useable_station _usable_station;
	station_interp _station_interp;
	
};
