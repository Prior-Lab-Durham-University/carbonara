#ifndef MUTUAL_WIND2_H
#define MUTUAL_WIND2_H

#include <cmath>
#include <complex>
#include "point.h"
#include "binaryFind.h"
#include <iostream>
#include <fstream>


class windingStorage{
public:
  windingStorage(std::vector<point> &sublist1In,std::vector<point> &sublist2In,bool &getWindIn,int &sigma1In,int &sigma2IN);
  windingStorage(){};
  std::vector<point> getSublist1();
  std::vector<point> getSublist2();
  bool getIsWind();
  int getSigma1();
  int getSigma2();
private:  
  std::vector<point> sublist1;
  std::vector<point> sublist2;
  bool getWind;
  int sigma1,sigma2;
};

class mutualWind2{
	public:
	mutualWind2();
        point interp(point& lower,point& upper,double zval);
        point getInterpolatedTop(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList);
        std::pair<double,double> getExtremaZ(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList);
        int quadrant(point &pointdif);
        point getTurningAngleQuad(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList);
        double angle(double xdif,double ydif);
        void branchTracker(int& currBranch,int& prevBranch,int& windingSum,double &rotDirec);
        bool checkAlligned(std::vector<point> &sublist1,std::vector<point> &sublist2);
        std::pair<int,point> getNextIndex(std::vector<point> &sublist1,std::vector<point> &sublist2,std::pair<int,int> &indexPair,int &prevQuad);
        void AppendExtremumJoin(std::vector<point> &sublist1,std::vector<point> &sublist2,point &newS1minExtremum,point &newS2maxExtremum,point &newJointExtremum);
        void AppendExtremumJoinClosed(std::vector<point> &sublist1,std::vector<point> &sublist2,point &newS1minExtremum,point &newS2maxExtremum,point &newJointExtremum);
	void checkAppendExtrema(std::vector<point> &sublist,point &minExtremum,point &maxExtremum);
	std::vector<std::vector<point> > interpolateForWinding(std::vector<point> &curveList,std::vector<point> &tanlist,std::vector<std::pair<int,int> > turningList);
	std::vector<windingStorage>  getMutualSections(std::vector<std::vector<point> > &superSublist1,std::vector<std::vector<point> > &superSublist2,int &firstIndex);
        int getIntegerWindingOfPairJoined(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp);
        int getIntegerWindingOfPairJoinedClosedEnd(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp);
        int getIntegerWindingOfPair(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp);
	 int getIntegerWindingOfPairLink(std::vector<point> &curveList1,std::vector<point> &curveList2,int &s1minInp,int &s1maxInp,int &s2minInp,int &s2maxInp);
	 int getIntegerWindingOfPairLinkTesting(windingStorage &w);
        int getIntegerWindingsStar(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList);
  void extendAboveSpan(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2);
  void extendBelowSpan(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2);
	void extendBothPosNeg(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2);
	void extendBothNegPos(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2);
  std::vector<int>   getIntegerWindingsStarGen(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList);
  std::vector<int> getIntegerWindingsClosed(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList);
  std::vector<int>  getIntegerWindingsStarGenStandard(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList,std::vector<double> &startAngles,std::vector<double> &endAngles);
  int getWindingInteger(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<int> &turnList1,std::vector<int> &turnList2);
};
#endif 
