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
