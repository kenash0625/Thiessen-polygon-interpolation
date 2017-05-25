#include <vector>
#include <string>
#include <map>

using namespace std;
class OGRGeometry;


struct st2ws_ws {
	st2ws_ws();
	map<string, double> stcdweight;
	double p, area;
}; 
struct st2ws_st {
	st2ws_st();
	double x, y, z;
	OGRGeometry *pGeom;
};
//添加所有流域测站 调用一次 
//要查坐标
void add_st(string &stcd, double x, double y);
//放到vector
void add_ws(string &wscds);
void set_ws_shp(const string &spath);;
//循环设置所有测站雨量
void begin_st_rain();
void st_rain(const string &, double);
map<string, st2ws_ws> & calc(const string &a, int &useblest, double &);
map<string, st2ws_ws> & rand_run(const string &r,int &u,double &d);
void wsfile(const string &);
void SelfIntersectGeom(const string &fn);
map<string, st2ws_st> &get_stcds2();

