#ifndef LOCAL_WRITHE_H
#define LOCAL_WRITHE_H

#include "point.h"
#include "getSections.h"
#include "mutualWind2.h"
#include <complex>

class localWrithe{
public:
  localWrithe();
  void localDensity(point& p1,point& p2,point& p3,point& p4);
  void rotateListToMaximum(std::vector<point> &list);
  void derivList(std::vector<point>& list);
  point interpolateTang(int lowIndex,double t);
  void reconstructCurve(std::vector<point>& list);
  point interp(point& lower,point& upper,double zval);
  void clearVecs();
  std::vector<point> derivListReturn(std::vector<point>& list);
  void localDensity(point& tangent,point& tderiv);
  void addInterpolatedTan(point& tan1,double& den1, point tan2,double& den2,double& integral);
  void addInterpolatedTanLocal(point& tan1,double& den1, point tan2,double& den2,double& integral);
  void FillDensityList();
  void SimpsonsRule(int startval,int endval);
  double getWrithePl(std::vector<point>& list,std::vector<int> turnlist);
  void getListTurnPts(std::vector<point> list);
  std::vector<int> getListTurnPtsReturn(std::vector<point> list);
  double angle(double xdif,double ydif);
  void checkTurn(point& n1,int n1pos,point& n2,int n2pos,std::vector<int> &turnlist,bool &wasPrev);
  void checkTurnWinding(point& nminus,point& nplus,point &ncurr,int &n,std::vector<int> &turnlist);
  std::vector<std::pair<int,int> > makeTurnForNonLocal();
  std::vector<std::pair<int,int> > makeTurnForNonLocal(std::vector<int> &turnListIn);
  double getWrithe(std::vector<point>& list);
  double getWritheLocal(std::vector<point>& list);
  std::vector<double> getWritheGen(std::vector<point>& list);
  std::vector<double> getWritheGenSmooth(std::vector<point>& list);
  std::vector<double> getWritheGenClosed(std::vector<point>& list);
  std::vector<double> getWritheGenStandard(std::vector<point>& list);
  std::vector<double> getWritheGenStandardSmooth(std::vector<point>& list);
  double getEndAngles(std::vector<point> &curvelist1,std::vector<point> &curvelist2,std::vector<int> & turnlist1,std::vector<int> & turnlist2);
  double getWinding(std::vector<point>& list1,std::vector<point>& list2);
  double getWindingSmooth(std::vector<point>& list1,std::vector<point>& list2);
  point dbold(std::vector<point>& pointList,int size,int i,int j);
  point dboldLk(std::vector<point>& pointList1,std::vector<point>& pointList2,int size1,int size2,int i,int j);
  point te(std::vector<point>& pointList,int size,int i);
  double mu(std::vector<point>& pointList,int size,int i,int k,int m, int j);
  double wij(std::vector<point>& pointList,int size,int i,int j);
  double lkij(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i,int j);
  double acosC(double temp);
  double DI(std::vector<point>& pointList);
  double DIAbs(std::vector<point>& pointList);
  //std::vector<double>  DI(std::vector<point>& pointList);
  double DIClosed(std::vector<point>& pointList);
  double muLink(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i, int j);
  double DIClosedLk(std::vector<point>& pointList,std::vector<point>& pointList2);
  std::vector<std::pair<std::pair<int,int>,double> > DIDownSample(std::vector<std::vector<point> >& pointListIn);
  std::vector<std::pair<std::pair<int,int>,double> > DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn);
  void  DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
  void  DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
private:
  double den,lambda;
  double writhepl;
  std::vector<point> dlist;
  std::vector<point> tanlist;
  std::vector<point> tderiv;
  std::vector<double> densityList;
  std::vector<int> turnlist;
  std::vector<int> turnPts;
  std::vector<double> arcLengthLst;
};

#endif
