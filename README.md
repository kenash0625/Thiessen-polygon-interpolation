# fantastic-guacamole
In order to achieve accurate estimation of the spatial distribution of rainfall, it is
necessary to use interpolation methods,this is the area weighted Thiessen polygon method.
input:point layer and polygon layer 
output:thiessen polygon and weight

dependencies  gdal2.3.2 geos3.6.3

Edited geos-3.6.3\src\triangulate\quadedge\QuadEdgeSubdivision.cpp line 588 ,added the following :
	
		cellPoly->setUserData(reinterpret_cast<void*>(geomFact.createPoint(c)));
		return cellPoly;

Edited geos-3.6.3\src\triangulate\VoronoiDiagramBuilder.cpp line 136 ,moved one line:

		else if(clipEnv.intersects(g->getEnvelopeInternal()))
		{
			result.reset( clipPoly->intersection(g) );
		}
		
		if(result.get() && !result->isEmpty() )
		{result->setUserData(((Geometry*)g)->getUserData()); // moved
			clipped->push_back(result.release());
		}

Usage:(edit paths)
D:\MyFirstProject\fantastic-guacamole.git\trunk\Debug\ConsoleApplication2.exe D:\MyFirstProject\fantastic-guacamole.git\trunk\Debug\test\jczd.shp D:\MyFirstProject\fantastic-guacamole.git\trunk\Debug\test\wata.shp

