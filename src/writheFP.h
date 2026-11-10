#ifndef WRITHEFP_H
#define WRITHEFP_H

#include "point.h"
#include "ktlMoleculeRandom.h"
#include "skmt.h"

class writheFP{
public:
  writheFP();
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
  double DIDownSampleSingle(std::vector<std::vector<point> >& pointListIn);
  double DIDownSampleAbsSingle(std::vector<std::vector<point> >& pointListIn);
  double DIDownSampleAbsSingleSKMT(ktlMolecule&molin);
  void  DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
  void  DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename);
  std::pair<double,double> sumWr(int& iv,int& jv,std::vector<std::vector<double> >& wijVals);
  std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > > getWritheFingerprints(std::vector<point>& pointList);
  double compare(std::vector<point> & sec1,std::vector<point> & sec2,int& maxSz);
  bool checkOverlap(std::pair<std::pair<int,int>,std::pair<int,int> >& p1,std::pair<std::pair<int,int>,std::pair<int,int> >& p2);
  bool checkOverlap(std::pair<int,int>& p1,std::pair<int,int>& p2);
  std::pair<bool,double>  checkOverlapLargest(std::pair<int,int>& p1,std::pair<int,int>& p2);
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > compareFingerPrints(std::vector<std::vector<point> >& fp1In,int& len1,std::vector<std::vector<point> >& fp2In,int& len2,double& cutOff);
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > finalSort(std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > & prSet);
  std::pair<int,std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > findBiggestOverlap(std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > >& separatedCombo,std::vector<int>& perm);
  std::vector<std::vector<int> > getAllSubsets(std::vector<int>& set);
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
