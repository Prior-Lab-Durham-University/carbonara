#ifndef BINARY_FIND_H
#define BINARY_FIND_H

#include <algorithm>
#include "point.h"

class binaryFind{
	public:
	binaryFind();
	void checkInSec(double lowerSec, double upperSec,double& zval);
	bool checkInSecInitial(double lowerSec, double upperSec,double& zv);
	std::pair<int,int> getContainingPair(std::vector<point>& pointList,int& lowerIndex,int& upperIndex,double& zv);
	private:
	bool isIn;	
	int lowerIndex,upperIndex;
	std::pair<int,int> boundingIndicies;
	double zval;
};

#endif
