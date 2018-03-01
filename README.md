# fantastic-guacamole
In order to achieve accurate estimation of the spatial distribution of rainfall, it is
necessary to use interpolation methods,this is the area weighted Thiessen polygon method.
input:station rain 
output watershed rain

using  gdal2.1.0 geos3.6.1

Edited geos-3.6.1\src\triangulate\quadedge\QuadEdgeSubdivision.cpp line 588 ,added the following :
	
		cellPoly->setUserData(reinterpret_cast<void*>(geomFact.createPoint(c)));
		return cellPoly;

Edited geos-3.6.1\src\triangulate\VoronoiDiagramBuilder.cpp line 136 ,moved one line:

		else if(clipEnv.intersects(g->getEnvelopeInternal()))
		{
			result.reset( clipPoly->intersection(g) );
		}
		
		if(result.get() && !result->isEmpty() )
		{result->setUserData(((Geometry*)g)->getUserData()); // moved
			clipped->push_back(result.release());
		}

How to use:
1.st2ws_wsshp    remove self-intersect polygons and multipolygons
2.
//add watershed. use one string to represent one watershed  
st2ws::add_ws
//add station.use one string to represent one station. do not add station with same x,y
st2ws::add_st 
//allocate memort
st2ws::add_end
st2ws::wsfile = wsShp
//calculate watershed area
st2ws::ws_shp_calc_area(wsShp);
//set station rain
st2ws::begin_st_rain
st2ws::st_rain
//calculate watershe rain
st2ws::calc
//get watershed rain
st2ws.ws_rain


