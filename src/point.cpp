#include "point.h"

point::point(double x,double y,double z){
      X =x;
      Y =y;
      Z =z;                                  
}

point::point(){};

point::point(std::string& triplet){
	std::stringstream ss( triplet );
    double d=0.0;
    int iv=1;
    std::vector<double> dv;
    while(ss>> d){
      dv.push_back(d);
    }
    X=dv[0];Y=dv[1];Z=dv[2];
}

double point::getX(){
      return X;       
}

double point::getY(){
      return Y;       
}

double point::getZ(){
       return Z;
}

void point::setX(double xv){
      X=xv;       
}

void point::setY(double yv){
      Y=yv;       
}

void point::setZ(double zv){
      Z=zv;
}

double point::length(){
  return std::sqrt(X*X+Y*Y +Z*Z); 
}

point point::sum(point& b){
	return point(b.getX() +X,b.getY()+Y,b.getZ()+Z);
}

point point::dif(point& b){
	return point(X-b.getX(),Y- b.getY(),Z - b.getZ());
}


point  point::cross(point& p2){
		return point(Y*p2.getZ()-Z*p2.getY(),Z*p2.getX()-X*p2.getZ(),X*p2.getY()-Y*p2.getX()); 
}
double  point::dotprod(point& p2){
		return X*p2.getX()+Y*p2.getY()+Z*p2.getZ(); 
}

double  point::scalarTriple(point&p1,point& p2,point& p3){
  point cp = p2.cross(p3);
  return p1.dotprod(cp);
}


void point::scalarMult(double a){
	X = X*a;
	Y = Y*a;
	Z = Z*a;
}

bool point::checkEqual(point& p2){
	bool isEqual = false;
	if(X == p2.getX()&& Y == p2.getY()&& Z == p2.getZ()){
		isEqual = true;
	}
	return isEqual;
}

void point::normalise(){
	norm = sqrt(X*X + Y*Y +Z*Z);
	if(norm != 0){
	X = X/norm;
	Y = Y/norm;
	Z = Z/norm;
	}
}

void point::znormalise(){
  if(std::abs(Z)>0.00000001){
	X = X/Z;
	Y = Y/Z;
	Z = 1.0;
     }
}

double point::pairDist(point &p){
  double xd = X -p.getX();
  double yd = Y -p.getY();
  double zd = Z -p.getZ();
  if(std::abs(xd) <0.00000001 && std::abs(yd) <0.00000001 && std::abs(zd) <0.00000001){
    return false;
  }else{
    return true;
  }
}

bool point::isNonzero(){
  if(std::abs(X) <0.00000001 && std::abs(Y) <0.00000001 && std::abs(Z) <0.00000001){
    return false;
  }else{
    return true;
  }
}

void point::printPoint(){
  std::cout<<X<<" "<<Y<<" "<<Z<<"\n";
}

/*point point::operator+(point &p){
  point pout;
  pout.setX(X + p.getX());
  pout.setY(Y + p.getY());
  pout.setZ(Z + p.getZ());
  return pout;
  }*/

point point::operator+(point p){
  point pout;
  pout.setX(X + p.getX());
  pout.setY(Y + p.getY());
  pout.setZ(Z + p.getZ());
  return pout;
  }

/*point point::operator-(point &p){
  point pout;
  pout.setX(X - p.getX());
  pout.setY(Y - p.getY());
  pout.setZ(Z - p.getZ());
  return pout;
  }*/

point point::operator-(point p){
  point pout;
  pout.setX(X - p.getX());
  pout.setY(Y - p.getY());
  pout.setZ(Z - p.getZ());
  return pout;
}

point point::operator*(double d){
  point pout;
  pout.setX(d*X);
  pout.setY(d*Y);
  pout.setZ(d*Z);
  return pout;
}

point point::operator/(double d){
  point pout;
  pout.setX(X/d);
  pout.setY(Y/d);
  pout.setZ(Z/d);
  return pout;
};

/*point point::operator*(double d){
  point pout;
  pout.setX(d*X);
  pout.setY(d*Y);
  pout.setZ(d*Z);
  return pout;
  }*/

double point::eDist(point &p2){
  double xd = X-p2.getX();
  double yd = Y-p2.getY();
  double zd = Z-p2.getZ();
  return std::sqrt(xd*xd + yd*yd + zd*zd);
}




