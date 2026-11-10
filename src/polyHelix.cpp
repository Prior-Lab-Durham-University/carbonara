#include "polyHelix.h"

polyHelix::polyHelix(std::vector<double> &kL,std::vector<double> thL,std::vector<double> &lL){
  klist = kL;
  int paramLength = klist.size();
  tlist = thL;
  int paramLength2 = tlist.size();
  lList = lL;
  int paramLength3 = lList.size();
  try{
    if(paramLength != paramLength2 ||paramLength2 !=  paramLength3){
      throw "curvature torsion and length lists are not the same size, re-edit the files, because of this ";
    }
  }catch(const char* msg){
    std::cerr <<msg<<std::endl;
    std::terminate();
  }
  d1c = 1; 
  d2c = 0;
  d3c = 0;
  n1c = 0;
  n2c = 1;
  n3c = 0;
  b1c = 0;
  b2c = 0;
  b3c = 1;
  r1c = 0;
  r2c = 0;
  r3c = 0;
  // now initilise the derivative y vetcors
  d1kd.resize(paramLength);d1kd.resize(paramLength);n1kd.resize(paramLength);b1kd.resize(paramLength);d2kd.resize(paramLength);n2kd.resize(paramLength);b2kd.resize(paramLength);d3kd.resize(paramLength);n3kd.resize(paramLength);b3kd.resize(paramLength);r1kd.resize(paramLength);r2kd.resize(paramLength);r3kd.resize(paramLength);d1ckd.resize(paramLength);n1ckd.resize(paramLength);b1ckd.resize(paramLength);d2ckd.resize(paramLength);n2ckd.resize(paramLength);b2ckd.resize(paramLength);d3ckd.resize(paramLength);n3ckd.resize(paramLength);b3ckd.resize(paramLength);r1ckd.resize(paramLength);r2ckd.resize(paramLength);r3ckd.resize(paramLength);
d1td.resize(paramLength);d1td.resize(paramLength);n1td.resize(paramLength);b1td.resize(paramLength);d2td.resize(paramLength);n2td.resize(paramLength);b2td.resize(paramLength);d3td.resize(paramLength);n3td.resize(paramLength);b3td.resize(paramLength);r1td.resize(paramLength);r2td.resize(paramLength);r3td.resize(paramLength);d1ctd.resize(paramLength);n1ctd.resize(paramLength);b1ctd.resize(paramLength);d2ctd.resize(paramLength);n2ctd.resize(paramLength);b2ctd.resize(paramLength);d3ctd.resize(paramLength);n3ctd.resize(paramLength);b3ctd.resize(paramLength);r1ctd.resize(paramLength);r2ctd.resize(paramLength);r3ctd.resize(paramLength);
d1Ld.resize(paramLength);d1Ld.resize(paramLength);n1Ld.resize(paramLength);b1Ld.resize(paramLength);d2Ld.resize(paramLength);n2Ld.resize(paramLength);b2Ld.resize(paramLength);d3Ld.resize(paramLength);n3Ld.resize(paramLength);b3Ld.resize(paramLength);r1Ld.resize(paramLength);r2Ld.resize(paramLength);r3Ld.resize(paramLength);d1cLd.resize(paramLength);n1cLd.resize(paramLength);b1cLd.resize(paramLength);d2cLd.resize(paramLength);n2cLd.resize(paramLength);b2cLd.resize(paramLength);d3cLd.resize(paramLength);n3cLd.resize(paramLength);b3cLd.resize(paramLength);r1cLd.resize(paramLength);r2cLd.resize(paramLength);r3cLd.resize(paramLength);
  molSize =klist.size();
  for(int l=0;l<molSize;l++){
    d1ckd[l] = 1;d1ctd[l] = 1;d1cLd[l] = 1;  
    d2ckd[l] = 0;d2ctd[l] = 0;d2cLd[l] = 0; 
    d3ckd[l] = 0;d3ctd[l] = 0;d2cLd[l] = 0;
    n1ckd[l] = 0;n1ctd[l] = 0;n1cLd[l] = 0;
    n2ckd[l] = 1;n2ctd[l] = 1;n2cLd[l] = 1;
    n3ckd[l] = 0;n3ctd[l] = 0;n3cLd[l] = 0;
    b1ckd[l] = 0;b1ctd[l] = 0;b1cLd[l] = 0;
    b2ckd[l] = 0;b2ctd[l] = 0;b2cLd[l] = 0;
    b3ckd[l] = 1;b3ctd[l] = 1;b3cLd[l] = 1;
    r1ckd[l] = 0;r1ctd[l] = 0;r1cLd[l] = 0;
    r2ckd[l] = 0;r2ctd[l] = 0;r2cLd[l] = 0;
    r3ckd[l] = 0;r3ctd[l] = 0;r3cLd[l] = 0;
  }
  // finally we resize the derivative vetcors 
  fullKderivList.resize(molSize);
  fullTderivList.resize(molSize);
  fullLderivList.resize(molSize);
}

/*Second constructor for the case in which tangent's do not match*/


polyHelix::polyHelix(std::vector<double> &kL,std::vector<double> thL,std::vector<double> &lL,std::vector<point> &d1vec,std::vector<point> &nvec,std::vector<point> &bvec,std::vector<point> &rvec,std::vector<double> &sdvec,std::vector<int> &nptsVec){
  klist = kL;
  int paramLength = klist.size();
  tlist = thL;
  int paramLength2 = tlist.size();
  lList = lL;
  int paramLength3 = lList.size();
  d1list = d1vec;
  nlist = nvec;
  blist = bvec;
  rlist = rvec;
  sdlist = sdvec;
  nptslist = nptsVec;
  try{
    if(paramLength != paramLength2 ||paramLength2 !=  paramLength3){
      throw "curvature torsion and length lists are not the same size, re-edit the files, because of this ";
    }
  }catch(const char* msg){
    std::cerr <<msg<<std::endl;
    std::terminate();
  }
}


void polyHelix::updateFrenetMatParameters(double &k,double &t, double &s){
  asq = k*k + t*t;
  a = std::sqrt(asq);
  acub =  asq*a;
  double sinas=std::sin(a*s);
  double cosas=std::cos(a*s);
  a11 = (t*t + k*k*cosas)/asq;
  a12 = k*sinas/a;
  a13 = k*t*(1-cosas)/asq;
  a21 = -a12;
  a22 = cosas;
  a23 = t*sinas/a;
  a31 =  a13;
  a32 = -a23;
  a33 = (k*k + t*t*cosas )/asq;
  c1 = (a*s*t*t + k*k*sinas)/acub;
  c2 = k*(1-cosas)/asq;
  c3 = k*t*(a*s -sinas)/acub;
}

void polyHelix::updateFrenetMatParametersKderiv(double &k,double &t, double &s){
  asq = k*k + t*t;
  a = std::sqrt(asq);
  acub =  asq*a;
  aquin =  asq*asq;
  ahex =  asq*acub;
  double sinas=std::sin(a*s);
  double cosas=std::cos(a*s);
  a11kd = (2*k*t*t*(-1.0 +cosas) - s*a*k*k*k*sinas)/aquin;
  a12kd = (s*a*k*k*cosas+t*t*sinas)/acub;
  a13kd = ((k-t)*t*(k+t)*(cosas-1.0) +s*a*k*k*t*sinas)/aquin;
  a21kd = -a12kd;
  a22kd = -s*k*sinas/a;
  a23kd = k*t*(a*s*cosas-sinas)/acub;
  a31kd =  a13kd;
  a32kd = -a23kd;
  a33kd =  -k*t*t*(-2.0 +2.0*cosas +s*a*sinas)/aquin;
  c1kd = k*(-2.0*s*a*t*t+s*a*k*k*cosas -(k*k-2.0*t*t)*sinas)/ahex;
  c2kd = ((k-t)*(k+t)*(-1.0 +cosas)+s*a*k*k*sinas)/aquin;
  c3kd = -t*(s*a*(k-t)*(k+t)+s*a*k*k*cosas +(-2.0*k*k + t*t)*sinas)/ahex;
}

void polyHelix::updateFrenetMatParametersTderiv(double &k,double &t, double &s){
  asq = k*k + t*t;
  a = std::sqrt(asq);
  acub =  asq*a;
  aquin =  asq*asq;
  ahex =  asq*acub;
  double sinas=std::sin(a*s);
  double cosas=std::cos(a*s);
  a11td = -k*k*t*(-2.0 +2.0*cosas +s*a*sinas)/aquin;
  a12td = k*t*(a*s*cosas-sinas)/acub;
  a13td = (-(k-t)*k*(k+t)*(cosas-1.0) +s*a*k*t*t*sinas)/aquin;
  a21td = -a12td;
  a22td = -s*t*sinas/a;
  a23td = (t*t*a*s*cosas+k*k*sinas)/acub;
  a31td =  a13td;
  a32td = -a23td;
  a33td =  (2*k*k*t*(-1.0 + cosas) -s*a*t*t*t*sinas)/aquin;
  c1td = k*k*t*(s*a*(2.0 + cosas) -3.0*sinas)/ahex;
  c2td = k*t*(-2.0 +2.0*cosas+a*s*sinas)/aquin;
  c3td = -k*(s*a*(-k*k+t*t)+s*a*t*t*cosas +(k*k-2.0*t*t)*sinas)/ahex;
}

void polyHelix::updateFrenetMatParametersLderiv(double &k,double &t, double &s,double &L){
  asq = k*k + t*t;
  a = std::sqrt(asq);
  double sinas=std::sin(a*s);
  double cosas=std::cos(a*s);
  a11Ld = -s*k*k*sinas/(L*a);
  a12Ld = s*k*cosas/L;
  a13Ld = s*k*t*sinas/(a*L);
  a21Ld = -a12Ld;
  a22Ld = -s*a*sinas/L;
  a23Ld = s*t*cosas/L;
  a31Ld =  a13td;
  a32Ld = -a23td;
  a33Ld =  -s*t*t*sinas/(a*L);
  c1Ld = s*(t*t + k*k*cosas)/(asq*L);
  c2Ld = s*k*sinas/(a*L);
  c3Ld = -s*k*t*(-1 + cosas)/(asq*L);
}



void polyHelix::updateYvec(double &k,double &t, double &s){
  updateFrenetMatParameters(k,t,s);
  d1 = a13*b1c + a11*d1c + a12*n1c;
  n1 = a23*b1c + a21*d1c + a22*n1c;
  b1 = a33*b1c + a31*d1c + a32*n1c;
  d2 = a13*b2c + a11*d2c + a12*n2c;
  n2 = a23*b2c + a21*d2c + a22*n2c;
  b2 = a33*b2c + a31*d2c + a32*n2c;
  d3 = a13*b3c + a11*d3c + a12*n3c;
  n3 = a23*b3c + a21*d3c + a22*n3c;
  b3 = a33*b3c + a31*d3c + a32*n3c;
  r1 = r1c + b1c*c3 + d1c*c1 + c2*n1c;
  r2 = r2c + b2c*c3 + d2c*c1 + c2*n2c;
  r3 = r3c + b3c*c3 + d3c*c1 + c2*n3c;
}

void polyHelix::updateYvecGeneral(double &k,double &t,point &d,point &n,point &b,point &r,double &s){
  updateFrenetMatParameters(k,t,s);
  d1 = a13*b.getX() + a11*d.getX() + a12*n.getX();
  n1 = a23*b.getX() + a21*d.getX() + a22*n.getX();
  b1 = a33*b.getX() + a31*d.getX() + a32*n.getX();
  d2 = a13*b.getY() + a11*d.getY() + a12*n.getY();
  n2 = a23*b.getY() + a21*d.getY() + a22*n.getY();
  b2 = a33*b.getY() + a31*d.getY() + a32*n.getY();
  d3 = a13*b.getZ() + a11*d.getZ() + a12*n.getZ();
  n3 = a23*b.getZ() + a21*d.getZ() + a22*n.getZ();
  b3 = a33*b.getZ() + a31*d.getZ() + a32*n.getZ();
  r1 = r.getX() + b.getX()*c3 + d.getX()*c1 + c2*n.getX();
  r2 = r.getY() + b.getY()*c3 + d.getY()*c1 + c2*n.getY();
  r3 = r.getZ() + b.getZ()*c3 + d.getZ()*c1 + c2*n.getZ();
}

point polyHelix::getTangent(){
  point p(d1,d2,d3);
  return p;
};

point polyHelix::getNormal(){
  point p(n1,n2,n3);
  return p;
};

point polyHelix::getBinormal(){
  point p(b1,b2,b3);
  return p;
};

std::vector<point> polyHelix::getFrame(){
  point t(d1,d2,d3);point n(n1,n2,n3);point b(b1,b2,b3);
  std::vector<point> fr;fr.push_back(t);fr.push_back(n);fr.push_back(b);
  return fr;
}

point polyHelix::getCoord(){
  point p(r1,r2,r3);
  return p;
};

void polyHelix::updateYvecKderiv(double &k,double &t, double &s,int &derivNo){
  updateFrenetMatParametersKderiv(k,t,s);
  d1kd[derivNo] = a13kd*b1ckd[derivNo] + a11kd*d1ckd[derivNo] + a12kd*n1ckd[derivNo];
  n1kd[derivNo] = a23kd*b1ckd[derivNo] + a21kd*d1ckd[derivNo] + a22kd*n1ckd[derivNo];
  b1kd[derivNo] = a33kd*b1ckd[derivNo] + a31kd*d1ckd[derivNo] + a32kd*n1ckd[derivNo];
  d2kd[derivNo] = a13kd*b2ckd[derivNo] + a11kd*d2ckd[derivNo] + a12kd*n2ckd[derivNo];
  n2kd[derivNo] = a23kd*b2ckd[derivNo] + a21kd*d2ckd[derivNo] + a22kd*n2ckd[derivNo];
  b2kd[derivNo] = a33kd*b2ckd[derivNo] + a31kd*d2ckd[derivNo] + a32kd*n2ckd[derivNo];
  d3kd[derivNo] = a13kd*b3ckd[derivNo] + a11kd*d3ckd[derivNo] + a12kd*n3ckd[derivNo];
  n3kd[derivNo] = a23kd*b3ckd[derivNo] + a21kd*d3ckd[derivNo] + a22kd*n3ckd[derivNo];
  b3kd[derivNo] = a33kd*b3ckd[derivNo] + a31kd*d3ckd[derivNo] + a32kd*n3ckd[derivNo];
  //r1kd[derivNo] = r1ckd[derivNo] + b1ckd[derivNo]*c3kd + d1ckd[derivNo]*c1kd + c2kd*n1ckd[derivNo];
  //r2kd[derivNo] = r2ckd[derivNo] + b2ckd[derivNo]*c3kd + d2ckd[derivNo]*c1kd + c2kd*n2ckd[derivNo];
  //r3kd[derivNo] = r3ckd[derivNo] + b3ckd[derivNo]*c3kd + d3ckd[derivNo]*c1kd + c2kd*n3ckd[derivNo];
  r1kd[derivNo] = b1ckd[derivNo]*c3kd + d1ckd[derivNo]*c1kd + c2kd*n1ckd[derivNo];
  r2kd[derivNo] = b2ckd[derivNo]*c3kd + d2ckd[derivNo]*c1kd + c2kd*n2ckd[derivNo];
  r3kd[derivNo] = b3ckd[derivNo]*c3kd + d3ckd[derivNo]*c1kd + c2kd*n3ckd[derivNo];
}

void polyHelix::updateYvecTderiv(double &k,double &t, double &s,int &derivNo){
  updateFrenetMatParametersTderiv(k,t,s);
  d1td[derivNo] = a13td*b1ctd[derivNo] + a11td*d1ctd[derivNo] + a12td*n1ctd[derivNo];
  n1td[derivNo] = a23td*b1ctd[derivNo] + a21td*d1ctd[derivNo] + a22td*n1ctd[derivNo];
  b1td[derivNo] = a33td*b1ctd[derivNo] + a31td*d1ctd[derivNo] + a32td*n1ctd[derivNo];
  d2td[derivNo] = a13td*b2ctd[derivNo] + a11td*d2ctd[derivNo] + a12td*n2ctd[derivNo];
  n2td[derivNo] = a23td*b2ctd[derivNo] + a21td*d2ctd[derivNo] + a22td*n2ctd[derivNo];
  b2td[derivNo] = a33td*b2ctd[derivNo] + a31td*d2ctd[derivNo] + a32td*n2ctd[derivNo];
  d3td[derivNo] = a13td*b3ctd[derivNo] + a11td*d3ctd[derivNo] + a12td*n3ctd[derivNo];
  n3td[derivNo] = a23td*b3ctd[derivNo] + a21td*d3ctd[derivNo] + a22td*n3ctd[derivNo];
  b3td[derivNo] = a33td*b3ctd[derivNo] + a31td*d3ctd[derivNo] + a32td*n3ctd[derivNo];
  // r1td[derivNo] = r1ctd[derivNo] + b1ctd[derivNo]*c3td + d1ctd[derivNo]*c1td + c2td*n1ctd[derivNo];
  // r2td[derivNo] = r2ctd[derivNo] + b2ctd[derivNo]*c3td + d2ctd[derivNo]*c1td + c2td*n2ctd[derivNo];
  // r3td[derivNo] = r3ctd[derivNo] + b3ctd[derivNo]*c3td + d3ctd[derivNo]*c1td + c2td*n3ctd[derivNo];
  r1td[derivNo] =  b1ctd[derivNo]*c3td + d1ctd[derivNo]*c1td + c2td*n1ctd[derivNo];
  r2td[derivNo] =  b2ctd[derivNo]*c3td + d2ctd[derivNo]*c1td + c2td*n2ctd[derivNo];
  r3td[derivNo] =  b3ctd[derivNo]*c3td + d3ctd[derivNo]*c1td + c2td*n3ctd[derivNo];
}

void polyHelix::updateYvecLderiv(double &k,double &t, double &s,double &L,int &derivNo){
  updateFrenetMatParametersLderiv(k,t,s,L);
  d1Ld[derivNo] = a13Ld*b1cLd[derivNo] + a11Ld*d1cLd[derivNo] + a12Ld*n1cLd[derivNo];
  n1Ld[derivNo] = a23Ld*b1cLd[derivNo] + a21Ld*d1cLd[derivNo] + a22Ld*n1cLd[derivNo];
  b1Ld[derivNo] = a33Ld*b1cLd[derivNo] + a31Ld*d1cLd[derivNo] + a32Ld*n1cLd[derivNo];
  d2Ld[derivNo] = a13Ld*b2cLd[derivNo] + a11Ld*d2cLd[derivNo] + a12Ld*n2cLd[derivNo];
  n2Ld[derivNo] = a23Ld*b2cLd[derivNo] + a21Ld*d2cLd[derivNo] + a22Ld*n2cLd[derivNo];
  b2Ld[derivNo] = a33Ld*b2cLd[derivNo] + a31Ld*d2cLd[derivNo] + a32Ld*n2cLd[derivNo];
  d3Ld[derivNo] = a13Ld*b3cLd[derivNo] + a11Ld*d3cLd[derivNo] + a12Ld*n3cLd[derivNo];
  n3Ld[derivNo] = a23Ld*b3cLd[derivNo] + a21Ld*d3cLd[derivNo] + a22Ld*n3cLd[derivNo];
  b3Ld[derivNo] = a33Ld*b3cLd[derivNo] + a31Ld*d3cLd[derivNo] + a32Ld*n3cLd[derivNo];
  //r1Ld[derivNo] = r1cLd[derivNo] + b1cLd[derivNo]*c3Ld + d1cLd[derivNo]*c1Ld + c2Ld*n1cLd[derivNo];
  //r2Ld[derivNo] = r2cLd[derivNo] + b2cLd[derivNo]*c3Ld + d2cLd[derivNo]*c1Ld + c2Ld*n2cLd[derivNo];
  //r3Ld[derivNo] = r3cLd[derivNo] + b3cLd[derivNo]*c3Ld + d3cLd[derivNo]*c1Ld + c2Ld*n3cLd[derivNo];
  r1Ld[derivNo] = b1cLd[derivNo]*c3Ld + d1cLd[derivNo]*c1Ld + c2Ld*n1cLd[derivNo];
  r2Ld[derivNo] = b2cLd[derivNo]*c3Ld + d2cLd[derivNo]*c1Ld + c2Ld*n2cLd[derivNo];
  r3Ld[derivNo] = b3cLd[derivNo]*c3Ld + d3cLd[derivNo]*c1Ld + c2Ld*n3cLd[derivNo];
}

void polyHelix::updateYvecWithKderiv(double &k,double &t, double &s,int &derivNo){
  updateFrenetMatParameters(k,t,s);
  d1kd[derivNo] = a13*b1ckd[derivNo] + a11*d1ckd[derivNo] + a12*n1ckd[derivNo];
  n1kd[derivNo] = a23*b1ckd[derivNo] + a21*d1ckd[derivNo] + a22*n1ckd[derivNo];
  b1kd[derivNo] = a33*b1ckd[derivNo] + a31*d1ckd[derivNo] + a32*n1ckd[derivNo];
  d2kd[derivNo] = a13*b2ckd[derivNo] + a11*d2ckd[derivNo] + a12*n2ckd[derivNo];
  n2kd[derivNo] = a23*b2ckd[derivNo] + a21*d2ckd[derivNo] + a22*n2ckd[derivNo];
  b2kd[derivNo] = a33*b2ckd[derivNo] + a31*d2ckd[derivNo] + a32*n2ckd[derivNo];
  d3kd[derivNo] = a13*b3ckd[derivNo] + a11*d3ckd[derivNo] + a12*n3ckd[derivNo];
  n3kd[derivNo] = a23*b3ckd[derivNo] + a21*d3ckd[derivNo] + a22*n3ckd[derivNo];
  b3kd[derivNo] = a33*b3ckd[derivNo] + a31*d3ckd[derivNo] + a32*n3ckd[derivNo];
  r1kd[derivNo] = r1ckd[derivNo] + b1ckd[derivNo]*c3 + d1ckd[derivNo]*c1 + c2*n1ckd[derivNo];
  r2kd[derivNo] = r2ckd[derivNo] + b2ckd[derivNo]*c3 + d2ckd[derivNo]*c1 + c2*n2ckd[derivNo];
  r3kd[derivNo] = r3ckd[derivNo] + b3ckd[derivNo]*c3 + d3ckd[derivNo]*c1 + c2*n3ckd[derivNo];
}

void polyHelix::updateYvecWithTderiv(double &k,double &t, double &s,int &derivNo){
  updateFrenetMatParameters(k,t,s);
  d1td[derivNo] = a13*b1ctd[derivNo] + a11*d1ctd[derivNo] + a12*n1ctd[derivNo];
  n1td[derivNo] = a23*b1ctd[derivNo] + a21*d1ctd[derivNo] + a22*n1ctd[derivNo];
  b1td[derivNo] = a33*b1ctd[derivNo] + a31*d1ctd[derivNo] + a32*n1ctd[derivNo];
  d2td[derivNo] = a13*b2ctd[derivNo] + a11*d2ctd[derivNo] + a12*n2ctd[derivNo];
  n2td[derivNo] = a23*b2ctd[derivNo] + a21*d2ctd[derivNo] + a22*n2ctd[derivNo];
  b2td[derivNo] = a33*b2ctd[derivNo] + a31*d2ctd[derivNo] + a32*n2ctd[derivNo];
  d3td[derivNo] = a13*b3ctd[derivNo] + a11*d3ctd[derivNo] + a12*n3ctd[derivNo];
  n3td[derivNo] = a23*b3ctd[derivNo] + a21*d3ctd[derivNo] + a22*n3ctd[derivNo];
  b3td[derivNo] = a33*b3ctd[derivNo] + a31*d3ctd[derivNo] + a32*n3ctd[derivNo];
  r1td[derivNo] = r1ctd[derivNo] + b1ctd[derivNo]*c3 + d1ctd[derivNo]*c1 + c2*n1ctd[derivNo];
  r2td[derivNo] = r2ctd[derivNo] + b2ctd[derivNo]*c3 + d2ctd[derivNo]*c1 + c2*n2ctd[derivNo];
  r3td[derivNo] = r3ctd[derivNo] + b3ctd[derivNo]*c3 + d3ctd[derivNo]*c1 + c2*n3ctd[derivNo];
}

void polyHelix::updateYvecWithLderiv(double &k,double &t, double &s,double &L,int &derivNo){
  updateFrenetMatParameters(k,t,s);
  d1Ld[derivNo] = a13*b1cLd[derivNo] + a11*d1cLd[derivNo] + a12*n1cLd[derivNo];
  n1Ld[derivNo] = a23*b1cLd[derivNo] + a21*d1cLd[derivNo] + a22*n1cLd[derivNo];
  b1Ld[derivNo] = a33*b1cLd[derivNo] + a31*d1cLd[derivNo] + a32*n1cLd[derivNo];
  d2Ld[derivNo] = a13*b2cLd[derivNo] + a11*d2cLd[derivNo] + a12*n2cLd[derivNo];
  n2Ld[derivNo] = a23*b2cLd[derivNo] + a21*d2cLd[derivNo] + a22*n2cLd[derivNo];
  b2Ld[derivNo] = a33*b2cLd[derivNo] + a31*d2cLd[derivNo] + a32*n2cLd[derivNo];
  d3Ld[derivNo] = a13*b3cLd[derivNo] + a11*d3cLd[derivNo] + a12*n3cLd[derivNo];
  n3Ld[derivNo] = a23*b3cLd[derivNo] + a21*d3cLd[derivNo] + a22*n3cLd[derivNo];
  b3Ld[derivNo] = a33*b3cLd[derivNo] + a31*d3cLd[derivNo] + a32*n3cLd[derivNo];
  r1Ld[derivNo] = r1cLd[derivNo] + b1cLd[derivNo]*c3 + d1cLd[derivNo]*c1 + c2*n1cLd[derivNo];
  r2Ld[derivNo] = r2cLd[derivNo] + b2cLd[derivNo]*c3 + d2cLd[derivNo]*c1 + c2*n2cLd[derivNo];
  r3Ld[derivNo] = r3cLd[derivNo] + b3cLd[derivNo]*c3 + d3cLd[derivNo]*c1 + c2*n3cLd[derivNo];
}


void polyHelix::storeConfig(){
  point tangent(d1,d2,d3);
  point normal(n1,n2,n3);
  point binormal(b1,b2,b3);
  point p(r1,r2,r3);
  tanlist.push_back(tangent);
  normlist.push_back(normal);
  binormlist.push_back(binormal);
  ptlist.push_back(p);
}

void polyHelix::chartOneHelix(int &index,int &npts,const char* filename){
  double stepsize = lList[index]/npts;
  double arclength;
  int j=0;
  for(int i =0;i<=npts; i++){
    arclength  = stepsize*i;
    updateYvec(klist[index],tlist[index],arclength);
    storeConfig();
    } 
  writeCurrentHelix(index,filename);
  d1c = d1; 
  d2c = d2;
  d3c = d3;
  n1c = n1;
  n2c = n2;
  n3c = n3;
  b1c = b1;
  b2c = b2;
  b3c = b3;
  r1c = r1;
  r2c = r2;
  r3c = r3;
}

void polyHelix::chartHelix(int &npts,const char* filename){
  for(int j = 0;j<molSize;j++){
    chartOneHelix(j,npts,filename);
  } 
}

std::vector<point> polyHelix::returnCurve(int npts){
  return ptlist; 
}

double polyHelix::getTotalLength(int &upperIndex){
  double tl=0;
  for(int i=0;i<upperIndex;i++){
    tl = tl + lList[i];
  }
  return tl;
}

std::vector<point>  polyHelix::discretisePolyHelix(int npts){
  double totalLength = getTotalLength(molSize);
  double stepsize = totalLength/(npts-1);
  double currLength=0;
  double prevLength=0; 
  double lengthdif=0;
  int index = 0;
  bool isDone=false;
  for(int i=0;i<npts;i++){
    if(currLength <= lList[index]+prevLength){
      lengthdif = currLength-prevLength;
      updateYvec(klist[index],tlist[index],lengthdif);
      point p(r1,r2,r3);
      fullPtList.push_back(p);
      currLength = currLength + stepsize;
    }
    else{
      updateYvec(klist[index],tlist[index],lList[index]);
      index ++;
      if(index<molSize){
      prevLength = getTotalLength(index);
      d1c = d1; 
      d2c = d2;
      d3c = d3;
      n1c = n1;
      n2c = n2;
      n3c = n3;
      b1c = b1;
      b2c = b2;
      b3c = b3;
      r1c = r1;
      r2c = r2;
      r3c = r3;
	lengthdif = currLength-prevLength;
	updateYvec(klist[index],tlist[index],lengthdif);
	point p(r1,r2,r3);
	fullPtList.push_back(p);
	currLength = currLength + stepsize;}{
      }	
     }
  }  
  return fullPtList;
}

/*optimizationPoint polyHelix::discretisePolyHelixAndDerivatives(int npts){
  double totalLength = getTotalLength(molSize);
  double stepsize = totalLength/(npts-1);
  double currLength=0;
  double prevLength=0; 
  double lengthdif=0;
  int index = 0;
  for(int i=0;i<npts;i++){
    if(currLength <= lList[index]+prevLength){
      lengthdif = currLength-prevLength;
      // update the curve itself
      updateYvec(klist[index],tlist[index],lengthdif);
      point p(r1,r2,r3);
      fullPtList.push_back(p);
      // For sections less than <index> the derivatives are zero
      for(int l=0;l<index;l++){
      updateYvecWithKderiv(klist[index],tlist[index],lengthdif,l);
      updateYvecWithTderiv(klist[index],tlist[index],lengthdif,l);
      updateYvecWithLderiv(klist[index],tlist[index],lengthdif,lList[index],l);
      point pkderiv(r1kd[l],r2kd[l],r3kd[l]);
      fullKderivList[l].push_back(pkderiv);
      point ptderiv(r1td[l],r2td[l],r3td[l]);
      fullTderivList[l].push_back(ptderiv);
      point pLderiv(r1Ld[l],r2Ld[l],r3Ld[l]);
      fullLderivList[l].push_back(pLderiv);
      }
      updateYvecKderiv(klist[index],tlist[index],lengthdif,index);
      updateYvecTderiv(klist[index],tlist[index],lengthdif,index);
      updateYvecLderiv(klist[index],tlist[index],lengthdif,lList[index],index);
      point pkderiv(r1kd[index],r2kd[index],r3kd[index]);
      fullKderivList[index].push_back(pkderiv);
      point ptderiv(r1td[index],r2td[index],r3td[index]);
      fullTderivList[index].push_back(ptderiv);
      point pLderiv(r1Ld[index],r2Ld[index],r3Ld[index]);
      fullLderivList[index].push_back(pLderiv);
      // for sections after <index> we also update but not with the dervative matrix, just the usual frenet matrix
      if(index <= molSize-2){
      for(int k=index+1;k<molSize;k++){
      point pkderiv(0,0,0);
      fullKderivList[k].push_back(pkderiv);
      point ptderiv(0,0,0);
      fullTderivList[k].push_back(ptderiv);
      point pLderiv(0,0,0);
      fullLderivList[k].push_back(pLderiv);
      }
      }
      currLength = currLength + stepsize;
    }
    else{
      updateYvec(klist[index],tlist[index],lList[index]);
      d1c = d1; 
      d2c = d2;
      d3c = d3;
      n1c = n1;
      n2c = n2;
      n3c = n3;
      b1c = b1;
      b2c = b2;
      b3c = b3;
      r1c = r1;
      r2c = r2;
      r3c = r3;
      for(int l=0;l<index;l++){
      updateYvecWithKderiv(klist[index],tlist[index],lList[index],l);
      d1ckd[l] = d1kd[l];
      d2ckd[l] = d2kd[l];
      d3ckd[l] = d3kd[l]; 
      n1ckd[l] = n1kd[l];
      n2ckd[l] = n2kd[l];
      n3ckd[l] = n3kd[l];
      b1ckd[l] = b1kd[l];
      b2ckd[l] = b2kd[l];
      b3ckd[l] = b3kd[l];
      r1ckd[l] = r1kd[l];
      r2ckd[l] = r2kd[l];
      r3ckd[l] = r3kd[l];
      updateYvecWithTderiv(klist[index],tlist[index],lList[index],l);
      d1ctd[l] = d1td[l];
      d2ctd[l] = d2td[l];
      d3ctd[l] = d3td[l]; 
      n1ctd[l] = n1td[l];
      n2ctd[l] = n2td[l];
      n3ctd[l] = n3td[l];
      b1ctd[l] = b1td[l];
      b2ctd[l] = b2td[l];
      b3ctd[l] = b3td[l];
      r1ctd[l] = r1td[l];
      r2ctd[l] = r2td[l];
      r3ctd[l] = r3td[l];      
      updateYvecWithLderiv(klist[index],tlist[index],lList[index],lList[index],l);
      d1cLd[l] = d1Ld[l];
      d2cLd[l] = d2Ld[l];
      d3cLd[l] = d3Ld[l]; 
      n1cLd[l] = n1Ld[l];
      n2cLd[l] = n2Ld[l];
      n3cLd[l] = n3Ld[l];
      b1cLd[l] = b1Ld[l];
      b2cLd[l] = b2Ld[l];
      b3cLd[l] = b3Ld[l];
      r1cLd[l] = r1Ld[l];
      r2cLd[l] = r2Ld[l];
      r3cLd[l] = r3Ld[l];  
      }
      updateYvecKderiv(klist[index],tlist[index],lList[index],index);
      d1ckd[index] = d1kd[index];
      d2ckd[index] = d2kd[index];
      d3ckd[index] = d3kd[index]; 
      n1ckd[index] = n1kd[index];
      n2ckd[index] = n2kd[index];
      n3ckd[index] = n3kd[index];
      b1ckd[index] = b1kd[index];
      b2ckd[index] = b2kd[index];
      b3ckd[index] = b3kd[index];
      r1ckd[index] = r1kd[index];
      r2ckd[index] = r2kd[index];
      r3ckd[index] = r3kd[index];
      updateYvecTderiv(klist[index],tlist[index],lList[index],index);
      d1ctd[index] = d1td[index];
      d2ctd[index] = d2td[index];
      d3ctd[index] = d3td[index]; 
      n1ctd[index] = n1td[index];
      n2ctd[index] = n2td[index];
      n3ctd[index] = n3td[index];
      b1ctd[index] = b1td[index];
      b2ctd[index] = b2td[index];
      b3ctd[index] = b3td[index];
      r1ctd[index] = r1td[index];
      r2ctd[index] = r2td[index];
      r3ctd[index] = r3td[index];
      updateYvecLderiv(klist[index],tlist[index],lList[index],lList[index],index);
      d1cLd[index] = d1Ld[index];
      d2cLd[index] = d2Ld[index];
      d3cLd[index] = d3Ld[index]; 
      n1cLd[index] = n1Ld[index];
      n2cLd[index] = n2Ld[index];
      n3cLd[index] = n3Ld[index];
      b1cLd[index] = b1Ld[index];
      b2cLd[index] = b2Ld[index];
      b3cLd[index] = b3Ld[index];
      r1cLd[index] = r1Ld[index];
      r2cLd[index] = r2Ld[index];
      r3cLd[index] = r3Ld[index];
      if(index <= molSize-2){
      for(int k=index+1;k<molSize;k++){
      updateYvecWithKderiv(klist[index],tlist[index],lList[index],k);
      d1ckd[k] = d1kd[k];
      d2ckd[k] = d2kd[k];
      d3ckd[k] = d3kd[k]; 
      n1ckd[k] = n1kd[k];
      n2ckd[k] = n2kd[k];
      n3ckd[k] = n3kd[k];
      b1ckd[k] = b1kd[k];
      b2ckd[k] = b2kd[k];
      b3ckd[k] = b3kd[k];
      r1ckd[k] = r1kd[k];
      r2ckd[k] = r2kd[k];
      r3ckd[k] = r3kd[k];
      updateYvecWithTderiv(klist[index],tlist[index],lList[index],k);
      d1ctd[k] = d1td[k];
      d2ctd[k] = d2td[k];
      d3ctd[k] = d3td[k]; 
      n1ctd[k] = n1td[k];
      n2ctd[k] = n2td[k];
      n3ctd[k] = n3td[k];
      b1ctd[k] = b1td[k];
      b2ctd[k] = b2td[k];
      b3ctd[k] = b3td[k];
      r1ctd[k] = r1td[k];
      r2ctd[k] = r2td[k];
      r3ctd[k] = r3td[k];  
      updateYvecWithLderiv(klist[index],tlist[index],lList[index],lList[index],k);
      d1cLd[k] = d1Ld[k];
      d2cLd[k] = d2Ld[k];
      d3cLd[k] = d3Ld[k]; 
      n1cLd[k] = n1Ld[k];
      n2cLd[k] = n2Ld[k];
      n3cLd[k] = n3Ld[k];
      b1cLd[k] = b1Ld[k];
      b2cLd[k] = b2Ld[k];
      b3cLd[k] = b3Ld[k];
      r1cLd[k] = r1Ld[k];
      r2cLd[k] = r2Ld[k];
      r3cLd[k] = r3Ld[k];       
      }
      }
      index ++;
      if(index<molSize){
      prevLength = getTotalLength(index);
      lengthdif = currLength-prevLength;
      updateYvec(klist[index],tlist[index],lengthdif);
      point p(r1,r2,r3);
      fullPtList.push_back(p);
      for(int l=0;l<index;l++){
      updateYvecWithKderiv(klist[index],tlist[index],lengthdif,l);
      updateYvecWithTderiv(klist[index],tlist[index],lengthdif,l);
      updateYvecWithLderiv(klist[index],tlist[index],lengthdif,lList[index],l);
      point pkderiv(r1kd[l],r2kd[l],r3kd[l]);
      fullKderivList[l].push_back(pkderiv);
      point ptderiv(r1td[l],r2td[l],r3td[l]);
      fullTderivList[l].push_back(ptderiv);
      point pLderiv(r1Ld[l],r2Ld[l],r3Ld[l]);
      fullLderivList[l].push_back(pLderiv);
      }
      updateYvecKderiv(klist[index],tlist[index],lengthdif,index);
      updateYvecTderiv(klist[index],tlist[index],lengthdif,index);
      updateYvecLderiv(klist[index],tlist[index],lengthdif,lList[index],index);
      point pkderiv(r1kd[index],r2kd[index],r3kd[index]);
      fullKderivList[index].push_back(pkderiv);
      point ptderiv(r1td[index],r2td[index],r3td[index]);
      fullTderivList[index].push_back(ptderiv);
      point pLderiv(r1Ld[index],r2Ld[index],r3Ld[index]);
      fullLderivList[index].push_back(pLderiv);
      // for sections after <index> we also update but not with the dervative matrix, just the usual frenet matrix
      if(index <= molSize-2){
      for(int k=index+1;k<molSize;k++){
      point pkderiv(0,0,0);
      fullKderivList[k].push_back(pkderiv);
      point ptderiv(0,0,0);
      fullTderivList[k].push_back(ptderiv);
      point pLderiv(0,0,0);
      fullLderivList[k].push_back(pLderiv);
      }
      }
      currLength = currLength + stepsize;
      }
    }
  }  
  return optimizationPoint(fullPtList,fullKderivList,fullTderivList,fullLderivList);
}
*/

void polyHelix::writeCurrentHelix(int &open,const char* filename){
  std::cout<<"tried to open ?"<<fullPtList.size()<<"\n";
  if(fullPtList.size()>0){
    std::ofstream outfile;
    if(open ==0){
      outfile.open(filename);
    }else{
      outfile.open(filename,std::ofstream::app);
    }
    if(outfile.fail()){
      std::cout<<"file for polyhelix data has failed to open";
    }else{
      for(int i=0;i<fullPtList.size();i++){
	//	outfile <<tanlist[i].getX()<<" "<<normlist[i].getX()<<" "<<binormlist[i].getX()<<" "<<tanlist[i].getY()<<" "<<normlist[i].getY()<<" "<<binormlist[i].getY()<<" "<<tanlist[i].getZ()<<" "<<normlist[i].getZ()<<" "<<binormlist[i].getZ()<<" "<<fullPtList[i].getX()<<" "<<fullPtList[i].getY()<<" "<<fullPtList[i].getZ()<<"\n"*/
	outfile <<fullPtList[i].getX()<<" "<<fullPtList[i].getY()<<" "<<fullPtList[i].getZ()<<"\n";
      }
      //outfile <<"\n";
      outfile.close();
      /*
      ptlist.clear();
      tanlist.clear();
      normlist.clear();
      binormlist.clear();*/
    }
    std::cout<<" polyhelix has been charted and written to \n"<<filename<<"\n";
  }else{
    std::cout<<"No polyhelix has been charted\n";
  }
}

void polyHelix::writeCurrentHelixDeriv(const char* filename){
  char dat[10] = ".dat";
    std::ofstream outfile;
    std::vector<std::ofstream*> kderivFiles;
    std::vector<std::ofstream*> tderivFiles;
    std::vector<std::ofstream*> LderivFiles;
    for(int i=0;i<molSize;i++){
      std::string istr;
      std::ostringstream convert;   // stream used for the conversion
      convert << i+1;      // insert the textual representation of 'Number' in the characters in the stream
      istr = convert.str(); 
      char kfn[50] = "kderiv";
      char tfn[50] = "tderiv";
      char Lfn[50] = "Lderiv";
      strcat(kfn,filename);
      strcat(tfn,filename);
      strcat(Lfn,filename);
      strcat(kfn,istr.c_str());
      strcat(tfn,istr.c_str());
      strcat(Lfn,istr.c_str());
      strcat(kfn,dat);
      strcat(tfn,dat);
      strcat(Lfn,dat);
      kderivFiles.push_back(new std::ofstream(kfn,std::ios::out));
      tderivFiles.push_back(new std::ofstream(tfn,std::ios::out));
      LderivFiles.push_back(new std::ofstream(Lfn,std::ios::out));
    }
    for(int j=0;j<molSize;j++){
      std::vector<point> kps = fullKderivList[j];
      std::vector<point> tps = fullTderivList[j];
      std::vector<point> Lps = fullLderivList[j];
        for(int i=0;i<fullPtList.size();i++){
	  *(kderivFiles[j]) <<kps[i].getX()<<" "<<kps[i].getY()<<" "<<kps[i].getZ()<<"\n";
	  *(tderivFiles[j]) <<tps[i].getX()<<" "<<tps[i].getY()<<" "<<tps[i].getZ()<<"\n";
	  *(LderivFiles[j]) <<Lps[i].getX()<<" "<<Lps[i].getY()<<" "<<Lps[i].getZ()<<"\n";
        }
   }
      for(int i=0;i<molSize;i++){
      kderivFiles[i]->close();
      tderivFiles[i]->close();
      LderivFiles[i]->close();
    }
}




// take a list  of points on the domain [0,1] and evaluates the curve at these points

std::vector<point>  polyHelix::discretisePolyHelixList(std::vector<double>& parameterList){
  double totalLength = getTotalLength(molSize);
  double currLength=0;
  double prevLength=0; 
  double lengthdif=0;;
  int index = 0;
  for(int i=0;i<parameterList.size();i++){
    if(currLength <= lList[index]+prevLength){
      lengthdif = currLength-prevLength;
      updateYvec(klist[index],tlist[index],lengthdif);
      point p(r1,r2,r3);
      fullPtList.push_back(p);
      currLength = parameterList[i]*totalLength;
    }
    else{
      prevLength = getTotalLength(index);
      updateYvec(klist[index],tlist[index],lList[index]);
      index ++;
      if(index<molSize){
      d1c = d1; 
      d2c = d2;
      d3c = d3;
      n1c = n1;
      n2c = n2;
      n3c = n3;
      b1c = b1;
      b2c = b2;
      b3c = b3; 
      r1c = r1;
      r2c = r2;
      r3c = r3;
      lengthdif = currLength-prevLength;
      updateYvec(klist[index],tlist[index],lengthdif);
      point p(r1,r2,r3);
      fullPtList.push_back(p);
      currLength = parameterList[i]*totalLength;	
      }
    }
  }  
  return fullPtList;
}

std::vector<point>  polyHelix::discretisePolyHelixListGen(){
  double currLength=0;
  double sval;
  double prevLength=0; 
  double lengthdif=0;;
  int index = 0;
  for(int i=0;i<klist.size();i++){
    //std::cout<<"kappa vals "<<klist[i]<<"\n";
    //std::cout<<"n mols  vals "<<nptslist[i]<<"\n";
    if(i==0){
      for(int j=0;j<nptslist[i];j++){
	sval = sdlist[0]*j;
	updateYvecGeneral(klist[0],tlist[0],d1list[0],nlist[0],blist[0],rlist[0],sval);
	point p(r1,r2,r3);
	fullPtList.push_back(p);
      }
    }
    else{
      for(int j=1;j<nptslist[i];j++){
	  sval = sdlist[i]*j;
	  updateYvecGeneral(klist[i],tlist[i],d1list[i],nlist[i],blist[i],rlist[i],sval);
	  point p(r1,r2,r3);
	  fullPtList.push_back(p);
      }
    }
  }
  return fullPtList;
}
