#ifndef POLY_HELIX_H
#define POLY_HELIX_H

#include <cmath>
#include "optimizationPoint.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

class polyHelix{
public : 
  polyHelix(std::vector<double> &kL,std::vector<double> thL,std::vector<double> &lL);
  polyHelix(std::vector<double> &kL,std::vector<double> thL,std::vector<double> &lL,std::vector<point> &d1vec,std::vector<point> &nvec,std::vector<point> &bvec,std::vector<point> &rvec,std::vector<double> &sdvec,std::vector<int> &nptsVec);
  polyHelix(){};
  void updateFrenetMatParameters(double &kappa,double &tau, double &arclength);
  void updateFrenetMatParametersKderiv(double &k,double &t, double &s);
  void updateFrenetMatParametersTderiv(double &k,double &t, double &s);
  void updateFrenetMatParametersLderiv(double &k,double &t, double &s,double &L);
  void updateYvec(double &k,double &t, double &s);
  void updateYvecGeneral(double &k,double &t,point &d,point &n,point &b,point &r,double &s);
  void updateYvecKderiv(double &k,double &t, double &s,int &derivNo);
  void updateYvecTderiv(double &k,double &t, double &s,int &derivNo);
  void updateYvecLderiv(double &k,double &t, double &s,double &L,int &derivNo);
  void updateYvecWithKderiv(double &k,double &t, double &s,int &derivNo);
  void updateYvecWithTderiv(double &k,double &t, double &s,int &derivNo);
  void updateYvecWithLderiv(double &k,double &t, double &s,int &derivNo);
  void updateYvecWithLderiv(double &k,double &t, double &s,double &L,int &derivNo);
  void storeConfig();
  point getTangent();
  point getNormal();
  point getBinormal();
  std::vector<point> getFrame();
  point getCoord();
  void chartOneHelix(int &index,int &npts,const char* filename);
  void chartHelix(int &npts,const char* filename);
  std::vector<point> returnCurve(int npts);
  double getTotalLength(int &upperIndex);
  std::vector<point> discretisePolyHelix(int npts);
  std::vector<point> discretisePolyHelixList(std::vector<double>& parameterList);
  std::vector<point> discretisePolyHelixListGen();
  optimizationPoint discretisePolyHelixAndDerivatives(int npts);
  void writeCurrentHelix(int &open,const char* filename);
  void writeCurrentHelixDeriv(const char* filename);
private:
  double a,asq,acub,aquin,ahex,c1,c2,c3,a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double d1,n1,b1,d2,n2,b2,d3,n3,b3,r1,r2,r3,d1c,n1c,b1c,d2c,n2c,b2c,d3c,n3c,b3c,r1c,r2c,r3c;
  double a11kd,a12kd,a13kd,a21kd,a22kd,a23kd,a31kd,a32kd,a33kd,c1kd,c2kd,c3kd;
  double a11td,a12td,a13td,a21td,a22td,a23td,a31td,a32td,a33td,c1td,c2td,c3td;
  double a11Ld,a12Ld,a13Ld,a21Ld,a22Ld,a23Ld,a31Ld,a32Ld,a33Ld,c1Ld,c2Ld,c3Ld;
  std::vector<double> d1kd,n1kd,b1kd,d2kd,n2kd,b2kd,d3kd,n3kd,b3kd,r1kd,r2kd,r3kd,d1ckd,n1ckd,b1ckd,d2ckd,n2ckd,b2ckd,d3ckd,n3ckd,b3ckd,r1ckd,r2ckd,r3ckd;
  std::vector<double> d1td,n1td,b1td,d2td,n2td,b2td,d3td,n3td,b3td,r1td,r2td,r3td,d1ctd,n1ctd,b1ctd,d2ctd,n2ctd,b2ctd,d3ctd,n3ctd,b3ctd,r1ctd,r2ctd,r3ctd;
  std::vector<double> d1Ld,n1Ld,b1Ld,d2Ld,n2Ld,b2Ld,d3Ld,n3Ld,b3Ld,r1Ld,r2Ld,r3Ld,d1cLd,n1cLd,b1cLd,d2cLd,n2cLd,b2cLd,d3cLd,n3cLd,b3cLd,r1cLd,r2cLd,r3cLd;
  double stepSize;
  int molSize;
  std::vector<double> thetalist;
  std::vector<double> philist;
  std::vector<double> klist;
  std::vector<double> tlist;
  std::vector<double> lList;
  std::vector<point> ptlist;
  std::vector<point> tanlist;
  std::vector<point> normlist;
  std::vector<point> binormlist;
  std::vector<point> fullPtList;
  std::vector<point> d1list;
  std::vector<point> nlist;
  std::vector<point> blist;
  std::vector<point> rlist;
  std::vector<double> sdlist;
  std::vector<int> nptslist;
  std::vector<std::vector<point> > fullKderivList;
  std::vector<std::vector<point> > fullTderivList;
  std::vector<std::vector<point> > fullLderivList;
};

#endif
