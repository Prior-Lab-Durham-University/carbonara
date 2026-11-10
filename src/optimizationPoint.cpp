#include "optimizationPoint.h"

optimizationPoint::optimizationPoint(std::vector<point> &pointList,std::vector<std::vector<point> > &kderivList,std::vector<std::vector<point> > &tderivList,std::vector<std::vector<point> > &LderivList){
points = pointList;
kderivs = kderivList;
tderivs = tderivList;
Lderivs = LderivList;
}

std::vector<point> optimizationPoint::getPoints(){
return points;
}

std::vector<std::vector<point> > optimizationPoint::getKderivs(){
return kderivs;
}

std::vector<std::vector<point> > optimizationPoint::getTderivs(){
return tderivs;
}

std::vector<std::vector<point> > optimizationPoint::getLderivs(){
return Lderivs;
}
