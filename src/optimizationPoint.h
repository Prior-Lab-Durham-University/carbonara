#ifndef OPTPOINT_H
#define OPTPOINT_H

#include "point.h"
#include <vector>

class optimizationPoint{
public:
optimizationPoint(std::vector<point> &pointList,std::vector<std::vector<point> > &kderivList,std::vector<std::vector<point> > &tderivList,std::vector<std::vector<point> > &LderivList);
std::vector<point> getPoints();
std::vector<std::vector<point> > getKderivs();
std::vector<std::vector<point> > getTderivs();
std::vector<std::vector<point> > getLderivs();
private:
std::vector<point> points;
std::vector<std::vector<point> > kderivs;
std::vector<std::vector<point> > tderivs;
std::vector<std::vector<point> > Lderivs;
};

#endif
