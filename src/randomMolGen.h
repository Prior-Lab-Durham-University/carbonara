#ifndef RANDOM_MOL
#define RANDOM_MOL


#include "point.h"
#include <iostream>
#include <random>
#include "polyHelix.h"
#include <algorithm>
#include <iterator> 
#include <chrono>

class randomMol{
 public:
  randomMol(){};
  randomMol(double &rminIn,double &rmaxIn,double &lminIn);
  void setParams(double &rminIn,double &rmaxIn,double &lminIn);
  std::pair<double,double> kapTau(point &pt1,point &pt2,point &pt3,point &pt4);
  point binorm(point &T1,point &T2);
  point getNextPoint(point &pt1,point &pt2,point &pt3,double &r,double &kap,double &tau);
  point getNextPointTest(point &pt1,point &pt2,point &pt3,double &r,double &kap,double &tau);
  double interpVel2D(double &x0,double &y0,double &x,double &y,std::vector<std::vector<double> > &vx,double &d1,double &d2,int &nx,int &ny);
  double interpVel1D(double &x0,double &p,std::vector<double> &vx,double &d1,int &nx);
  double bisection(double (randomMol::*func)(double),double &val,int n,double rangeLower,double rangeUpper,double &prevVal);
  double bisectionDoubleArg(double (randomMol::*func)(double,double),double &val,int n,double rangeLower,double rangeUpper,double secondArg,double &prevVal);
  double cumDistMixed(double x);
  double cumDistNegBeta(double x);
  double cumDistPosBeta(double x);
  double cumDistAlpha(double x);
  double cumDistMixedStrand(double x);
  double cumDistNegBetaStrand(double x);
  double cumDistPosBetaStrand(double x);
  double cumDistAlphaStrand(double x);
  double cumDistLinkerToHelix(double x);
  double cumDistHelixToLinker(double x);
  double cumDistStrandToHelix(double x);
  double cumDistHelixToStrand(double x);
  double cumDistStrandToLinker(double x);
  double cumDistLinkerToStrand(double x);
  double cumXGYMixed(double Y,double kapValIn);
  double cumXGYNegBeta(double Y,double kapValIn);
  double cumXGYPosBeta(double Y,double kapValIn);
  double cumXGYAlpha(double Y,double kapValIn);
  double cumXGYMixedStrand(double Y,double kapValIn);
  double cumXGYNegBetaStrand(double Y,double kapValIn);
  double cumXGYPosBetaStrand(double Y,double kapValIn);
  double cumXGYAlphaStrand(double Y,double kapValIn);
  double cumXGYLinkerToHelix(double Y,double kapValIn);
  double cumXGYHelixToLinker(double Y,double kapValIn);
  double cumXGYStrandToHelix(double Y,double kapValIn);
  double cumXGYHelixToStrand(double Y,double kapValIn);
  double cumXGYStrandToLinker(double Y,double kapValIn);
  double cumXGYLinkerToStrand(double Y,double kapValIn);
  std::pair<double,double> getKapTauProb(int &type,std::default_random_engine &generator); 
  std::pair<double,double> getProbFromKapTau(int &type,std::default_random_engine &generator,std::pair<double,double> &kapTauPair,double probChange);
  std::vector<point> getTanForHelix(double &k,double &t,double &l,double &xv,double &yv,double &zv,point &p);
  point parallelTransport(point &tan1,point &tan2,point &norm1);
  std::vector<point> makeRandomSectionWithDist(point &stPt,point &stTan,point &stNorm,int &npts,std::default_random_engine &generator,int &loopOrStrand,bool &suceeded);
  std::vector<point> alterSectionWithDist(std::vector<point> currPts,int &strandOrLinker,double radOfSearch,std::default_random_engine &generator,bool &suceeded);
  std::vector<point>  blendLoopToHelix(point &ptloop1,point &ptloop2,point &ptloop3,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,bool &suceeded,int &index);
  std::vector<point> blendLoopToHelixAlter(point &ptloop1,point &ptloop2,point &ptloop3,point &ptloop4,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,double radOfSearch,bool &suceeded,int &index);
  std::vector<point>  blendLoopToHelixStrand(point &ptloop1,point &ptloop2,point &ptloop3,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,bool &suceeded);
  std::vector<point>  blendHelixToLoop(point &pthel1,point &pthel2,point &pthel3,int &nl,std::default_random_engine &generator,int loopOrStrand,bool &suceeded);
  std::vector<point> blendLoopWithStrand(point &ptloop1,point &ptloop2,point &ptloop3,int &nHel,std::default_random_engine &generator,int strandOrLoop,bool &suceeded);
  std::vector<point> getFrameHelix(std::vector<point> &coords);
  std::vector<point> blendHelixToHelix(std::vector<point> &prevSec,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,bool &suceeded,int &index);
  std::vector<point> blendHelixToHelixAlter(point &ptloop1,point &ptloop2,point &ptloop3,point &ptloop4,double k,double t,double al,double cl,int &nHel,std::default_random_engine &generator,double radOfSearch,bool &suceeded,int &index);
  std::vector<std::vector<point> > makeRandomMolecule(std::vector<std::pair<std::string,int> > &molDat,point &sp,bool &suceeded);
  void tanToSec(point &sp,std::vector<point> &prevSec);
  double getMaxDistChange(std::vector<point> &oldPts,std::vector<point> &newPoints);
  double reshapeMol(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,point &sp,bool &suceeded);
  double reshapeMolHelixSet(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,int noHelices,point &sp,bool &suceeded);
  double reshapeMolLoopThenHelixSet(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,int noHelices,point &sp,bool &suceeded);
  double reshapeMolSmallVariation(std::vector<std::pair<std::string,int> > &molDat,std::vector<std::vector<point> > &molPos,int &index,point &sp,double variationSize,bool &suceeded);
 private:
  std::vector<std::vector<double> > mixedYfromXInterpolants;
  std::vector<std::vector<double> > posBetaYfromXInterpolants;
  std::vector<std::vector<double> > negBetaYfromXInterpolants;
  std::vector<std::vector<double> > alphaYfromXInterpolants;
  std::vector<std::vector<double> > helixToLinkerYfromXInterpolants;
  std::vector<std::vector<double> > linkerToHelixYfromXInterpolants;
  std::vector<std::vector<double> > strandToLinkerYfromXInterpolants;
  std::vector<std::vector<double> > linkerToStrandYfromXInterpolants;
  std::vector<double> mixedXCumInterpolants;
  std::vector<double> posBetaXCumInterpolants;
  std::vector<double> negBetaXCumInterpolants;
  std::vector<double> alphaXCumInterpolants;
  std::vector<double> helixToLinkerXCumInterpolants;
  std::vector<double> linkerToHelixXCumInterpolants;
  std::vector<double> strandToLinkerXCumInterpolants;
  std::vector<double> linkerToStrandXCumInterpolants;
  std::vector<std::vector<double> > mixedStrandYfromXInterpolants;
  std::vector<std::vector<double> > posBetaStrandYfromXInterpolants;
  std::vector<std::vector<double> > negBetaStrandYfromXInterpolants;
  std::vector<std::vector<double> > alphaStrandYfromXInterpolants;
  std::vector<std::vector<double> > helixToStrandYfromXInterpolants;
  std::vector<std::vector<double> > strandToHelixYfromXInterpolants;
  std::vector<double> mixedStrandXCumInterpolants;
  std::vector<double> posBetaStrandXCumInterpolants;
  std::vector<double> negBetaStrandXCumInterpolants;
  std::vector<double> alphaStrandXCumInterpolants;
  std::vector<double> helixToStrandXCumInterpolants;
  std::vector<double> strandToHelixXCumInterpolants;
  std::uniform_real_distribution<double> distribution;
  //std::vector<std::vector<point> > frameSet;
  std::vector<point> dummyFrame;
  double kapL;double kapU;double tauL;double tauU;
  int nkap,ntau,nGridSize;
  double dkap,dtau,kapVal,pAlpha,pposBeta,pnegBeta,pMixed,pAlphaStrand,pposBetaStrand,pnegBetaStrand,pMixedStrand,lmin,rmin,rmax;
  double alphaKap,alphaTau,alphaAl,alphaCl,betaKap,betaTau,betaAl,betaCl;
};

#endif
