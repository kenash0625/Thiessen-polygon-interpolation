# fantastic-guacamole
给出站点坐标 和 流域边界
计算站点对应的泰森多边形 和 与每个流域的相交的多边形面积所占流域面积的权重
界面展示

In order to achieve accurate estimation of the spatial distribution of rainfall, it is
necessary to use interpolation methods,this is the area weighted Thiessen polygon method.

using wxWidgets gdal2.1.0 geos3.6.1

Edited geos-3.6.1\src\triangulate\quadedge\QuadEdgeSubdivision.cpp line 588 ,added the following :
	geos::geom::Point *pc = geomFact.createPoint(c);
	cellPoly->setUserData(reinterpret_cast<void*>(pc));
	return cellPoly;