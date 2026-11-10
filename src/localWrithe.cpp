#include "localWrithe.h"

const double PI  =3.141592653589793238463;

localWrithe::localWrithe(){
	writhepl =0;
}

void localWrithe::clearVecs(){
  dlist.clear();
  tanlist.clear();
  tderiv.clear();
  densityList.clear();
  turnlist.clear();
  turnPts.clear();
  arcLengthLst.clear();
}

point localWrithe::interp(point& lower,point& upper,double zval){
  double t;
  if(std::abs(upper.getZ() -lower.getZ()) > 0.0001){
  t = (zval -lower.getZ())/(upper.getZ() -lower.getZ());
  }else{
    t = 0.5;
  }
  double xval = lower.getX() + t*(upper.getX()-lower.getX());
  double yval = lower.getY() + t*(upper.getY()-lower.getY());
  point p(xval,yval,zval);
  return p;
}

void localWrithe::rotateListToMaximum(std::vector<point> &list){
   int maxZIndex=0;
   double maxZ = -1000000.0;
   for(int i=0;i<list.size();i++){
     if(list[i].getZ()>maxZ){
       maxZIndex =i;
       maxZ = list[i].getZ();
     }
   }
   // now rotate the list so that the min point is the first in the list
   std::rotate(list.begin(),list.begin()+maxZIndex,list.end());
}

void localWrithe::derivList(std::vector<point>& list){
	int npts = list.size();
	dlist.push_back(list[1].dif(list[0]));
	arcLengthLst.push_back(dlist[0].length());
	dlist[0].normalise();
	tanlist.push_back(dlist[0]);
	for(int i=1;i<npts-1;i++){
		dlist.push_back(list[i+1].dif(list[i-1]));
		arcLengthLst.push_back(0.5*dlist[i].length());
		dlist[i].scalarMult(0.5);
		dlist[i].normalise();
		tanlist.push_back(dlist[i]);
	}
	dlist.push_back(list[npts-1].dif(list[npts-2]));
	arcLengthLst.push_back(0.5*dlist[npts-1].length());
	dlist[npts-1].normalise();	
	tanlist.push_back(dlist[npts-1]);
	// now the tangent derivative 
	tderiv.push_back(tanlist[1].dif(tanlist[0]));
	tderiv[0].scalarMult(2);
	for(int i=1;i<npts-1;i++){
		tderiv.push_back(tanlist[i+1].dif(tanlist[i-1]));
		tderiv[i].scalarMult(0.5);
	}
	tderiv.push_back(tanlist[npts-1].dif(tanlist[npts-2]));	
	tderiv[npts-1].scalarMult(2);
}

point localWrithe::interpolateTang(int lowIndex,double t){
   double tcub =t*t*t;double tsq =t*t;
   point tan= tanlist[lowIndex]*(2.0*tcub-3.0*tsq +1.0) + tderiv[lowIndex]*(tcub-2.0*tsq + t) + tanlist[lowIndex+1]*(-2.0*tcub+3.0*tsq)+ tderiv[lowIndex+1]*(tcub-tsq);
   tan.normalise();
    return tan;
}


void localWrithe::reconstructCurve(std::vector<point>& list){
 double coeff= 1/12.0;
 for(int i=0;i<(tanlist.size()-1);i++){
   point k1 = tanlist[i];
   point k2 = interpolateTang(i,0.25)*4.0;
   point k3 = interpolateTang(i,0.5)*2.0;
   point k4 = interpolateTang(i,0.75)*4.0;
   point k5 = tanlist[i+1];
   list[i+1]=list[i]+(k1+k2+k3+k4+k5)*coeff*arcLengthLst[i];
   //now get 
 }
}

std::vector<point> localWrithe::derivListReturn(std::vector<point>& list){
  std::vector<point> dlistTemp;
  std::vector<point> tanlistTemp;
  int npts = list.size();
  //std::cout<<dlist.size()<<" "<<arcLengthLst.size()<<" "<<tanlist.size()<<"\n";
  dlist.push_back(list[1].dif(list[0]));
  arcLengthLst.push_back(dlist[0].length());
  dlist[0].normalise();
  tanlist.push_back(dlist[0]);
  for(int i=1;i<npts-1;i++){
    dlist.push_back(list[i+1].dif(list[i-1]));
    arcLengthLst.push_back(0.5*dlist[i].length());
    dlist[i].scalarMult(0.5);
    dlist[i].normalise();
    tanlist.push_back(dlist[i]);
  }
  dlist.push_back(list[npts-1].dif(list[npts-2]));
  arcLengthLst.push_back(0.5*dlist[npts-1].length());
  dlist[npts-1].normalise();	
  tanlist.push_back(dlist[npts-1]);
  // now the tangent derivative 
  tderiv.push_back(tanlist[1].dif(tanlist[0]));
  tderiv[0].scalarMult(2);
  for(int i=1;i<npts-1;i++){
    tderiv.push_back(tanlist[i+1].dif(tanlist[i-1]));
    tderiv[i].scalarMult(0.5);
  }
  tderiv.push_back(tanlist[npts-1].dif(tanlist[npts-2]));	
  tderiv[npts-1].scalarMult(2);
  return tanlist;
}


/*void localWrithe::checkTurn(point& n1,int n1pos,point& n2,int n2pos){
      double prod = n1.getZ()*n2.getZ();
      if(prod<=0){
         turnlist.push_back(n1pos);
         turnlist.push_back(n2pos);
      }
      }*/


void localWrithe::checkTurn(point& n1,int n1pos,point& n2,int n2pos,std::vector<int> &turnlist,bool &wasPrev){
  double prod = n1.getZ()*n2.getZ();
  if(prod<0.0){
    turnlist.push_back(n1pos);
    turnlist.push_back(n2pos);
  }else if(prod == 0 && wasPrev == false && n1.isNonzero()==true &&  n2.isNonzero()==true ){
    turnlist.push_back(n1pos);
    turnlist.push_back(n2pos);
    wasPrev =true;
  }else{
    wasPrev = false;
  }
}

void localWrithe::checkTurnWinding(point& nminus,point& nplus,point &ncurr,int &n,std::vector<int> &turnlist){
  //std::cout<<nminus.getZ()<<" "<<ncurr.getZ()<<" "<<nplus.getZ()<<"\n";
  double prod = (nplus.getZ()-ncurr.getZ())*(ncurr.getZ()-nminus.getZ());
  //std::cout<<prod<<"\n";
  if(prod<0.0){
    //std::cout<<"here "<<n<<"\n";
    turnlist.push_back(n);
  }
}


void localWrithe::getListTurnPts(std::vector<point> list){
            turnlist.push_back(0);
	    bool wasPrev  =false;
	    int npts = list.size();
	    for(int i=0;i<list.size()-1;i++){
	      checkTurn(list[i],i,list[i+1],i+1,turnlist,wasPrev);   
            }
            turnlist.push_back(npts-1);
}



std::vector<int> localWrithe::getListTurnPtsReturn(std::vector<point> list){
  std::vector<int> turnListOut;
  turnListOut.push_back(0);
  bool wasPrev  =false;
  int npts = list.size();
  for(int i=1;i<list.size()-1;i++){
   checkTurnWinding(list[i-1],list[i+1],list[i],i,turnListOut);   
  }
  turnListOut.push_back(npts-1);
  return turnListOut;
}


double localWrithe::angle(double xdif,double ydif){
	if(xdif==0 && ydif ==0){
	return 0;	
	}
	std::complex<double> imag(xdif,ydif);
	double angv = std::arg(imag);
	return angv;
}

void localWrithe::addInterpolatedTan(point& tan1,double& den1, point tan2,double& den2,double& integral){
	double trap1,trap2,t,turnAng;
	if(tan1.getZ()!= tan2.getZ()){
		t= -tan1.getZ()/(tan2.getZ()-tan1.getZ());
		//std::cout<<"t is "<<t<<"\n";
		if(t<-1||t>1){t=0.5;}
		double xnew = tan1.getX() + t*(tan2.getX()-tan1.getX());
		double ynew = tan1.getY() + t*(tan2.getY()-tan1.getY());
		turnAng = angle(xnew,ynew);
		//std::cout<<turnAng<<"\n";
		// multiply by -1 if a maxium
		if(tan1.getZ()>0 || tan2.getZ()<0){
			turnAng = -turnAng;
		}
	}else{
		t=0.5;
		double xnew = tan1.getX() + t*(tan2.getX()-tan1.getX());
		double ynew = tan1.getY() + t*(tan2.getY()-tan1.getY());
		turnAng = angle(xnew,ynew);
		//std::cout<<"out o bounds "<<turnAng<<"\n";
		// multiply by -1 if a maxium
		if(tan1.getZ()>0 || tan2.getZ()<0){
			turnAng = -turnAng;
		}
		trap1=0;trap2=0;
	}
	if(t!=0){
		double denmid = den1 + t*(-den2 - den1);
		trap1 = ((den1 + denmid)*t)/2;
		trap2 = ((den2 - denmid)*(1-t))/2;
	}else{
		trap1 = 0;trap2=0;
	}
	integral = integral + trap1 + trap2 +2*turnAng; 
}

void localWrithe::addInterpolatedTanLocal(point& tan1,double& den1, point tan2,double& den2,double& integral){
	double trap1,trap2,t,turnAng;
	double denmid = den1 + t*(-den2 - den1);
	trap1 = ((den1 + denmid)*t)/2;
	trap2 = ((den2 - denmid)*(1-t))/2;
	integral = integral + trap1 + trap2;
}

void localWrithe::localDensity(point& tangent,point& tderiv){
	lambda = tangent.getZ();
	if(lambda>0){
		den = (tangent.getX()*tderiv.getY() - tangent.getY()*tderiv.getX())/(1+std::abs(lambda));
	}
	else{
		den = (-1.0)*(tangent.getX()*tderiv.getY() - tangent.getY()*tderiv.getX())/(1+std::abs(lambda));
	}
}

void localWrithe::FillDensityList(){
	int npts = tanlist.size();
	for(int i=0;i<npts;i++){
		localDensity(tanlist[i],tderiv[i]);
		densityList.push_back(den);
	}
}

void localWrithe::SimpsonsRule(int startval,int endval){
	if(endval-startval>7){
		//if 8 or more data points use Simpson's rule
		writhepl = writhepl+ (17.0/48.0)*densityList[startval];
		writhepl = writhepl+ (59.0/48.0)*densityList[startval+1];
		writhepl = writhepl+ (43.0/48.0)*densityList[startval+2];
		writhepl = writhepl+ (49.0/48.0)*densityList[startval+3];
		for(int i=startval+4;i<endval-3;i++){
			writhepl = writhepl+ densityList[i];
		}
		writhepl = writhepl+ (49.0/48.0)*densityList[endval-3];;
		writhepl = writhepl+ (43.0/48.0)*densityList[endval-2];
		writhepl = writhepl+ (59.0/48.0)*densityList[endval-1];
		writhepl = writhepl+ (17.0/48.0)*densityList[endval];
	}else{
		// else use Trapezium rule
		writhepl = writhepl+ 0.25*densityList[startval];
		for(int i=startval+1;i<endval;i++){
			writhepl = writhepl+ 0.5*densityList[i];
		}
		writhepl = writhepl+ 0.25*densityList[endval];
	}
}

/*void localWrithe::makeTurnForNonLocal(std::vector<point>& list){
	turnPts.push_back(0);
	std::vector<point>::iterator it;
	it = list.begin();
	for(int i = 1;i<turnlist.size()-1;i=i+2){
		double t= -tanlist[turnlist[i]].getZ()/(tanlist[turnlist[i+1]].getZ()-tanlist[turnlist[i]].getZ());
		double xnew = list[turnlist[i]].getX() + t*(list[turnlist[i+1]].getX()-list[turnlist[i]].getX());
		double ynew = list[turnlist[i]].getY() + t*(list[turnlist[i+1]].getY()-list[turnlist[i]].getY());
		double znew = list[turnlist[i]].getZ() + t*(list[turnlist[i+1]].getZ()-list[turnlist[i]].getZ());
		point newpoint(xnew,ynew,znew);
		list.insert(it+turnlist[i+1],newpoint);
		turnPts.push_back(turnlist[i+1]);
		for(int j = i+2;j< turnlist.size();j++){
			turnlist[j] = turnlist[j]+1;
		}
	}
	turnPts.push_back(turnlist.back());
}*/

std::vector<std::pair<int,int> > localWrithe::makeTurnForNonLocal(){
  std::vector<std::pair<int,int> > turnPairs;
  for(int i = 1;i<turnlist.size()-2;i=i+2){
	  std::pair<int,int> turnPair;
	  turnPair.first =turnlist[i];
	  turnPair.second = turnlist[i+1];
	  turnPairs.push_back(turnPair);
	}
  return turnPairs;
}

std::vector<std::pair<int,int> > localWrithe::makeTurnForNonLocal(std::vector<int> &turnListIn){
  std::vector<std::pair<int,int> > turnPairs;
  for(int i = 1;i<turnListIn.size()-2;i=i+2){
	  std::pair<int,int> turnPair;
	  turnPair.first =turnListIn[i];
	  turnPair.second = turnListIn[i+1];
	  turnPairs.push_back(turnPair);
	}
  return turnPairs;
}


double localWrithe::getWrithe(std::vector<point>& list){
	derivList(list);
	FillDensityList();
	getListTurnPts(tanlist);
	// first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	writhepl = (1/(2*PI))*writhepl;
	//smooth out the curve
	reconstructCurve(list);
	// now get the integer winding
	std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
	int integerWindingTotal =mw.getIntegerWindingsStar(list,tanlist,pairedTurnList);
	clearVecs();
	return writhepl + 2.0*double(integerWindingTotal);
}


double localWrithe::getWritheLocal(std::vector<point>& list){
	derivList(list);
	FillDensityList();
	getListTurnPts(tanlist);
	// first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTanLocal(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	return (1/(2*PI))*writhepl;
}


std::vector<double> localWrithe::getWritheGen(std::vector<point>& list){
	derivList(list);
	//first remove any sharp tangent changes
	double prod;
	std::vector<std::pair<int,int> > tppairs;
	bool isFirst=true;
	std::pair<int,int> pr;
	for(int i=0;i<tanlist.size()-2;i++){
	  prod = (tanlist[i+2].getZ()-tanlist[i+1].getZ())*(tanlist[i+1].getZ()-tanlist[i].getZ());
	  if(isFirst){
	     if(prod<0){
	       isFirst=false;
	       pr.first=i+1;
	     }
	  }else{
	    //here we have foudn one tp
	    if(prod>0){
	      isFirst=true;
	      pr.second=i+2;
	      tppairs.push_back(pr);
	    }
	   }
	}
	for(int id=0;id<tppairs.size();id++){
	  if(tppairs[id].second-tppairs[id].first>2){
	    point p1 = tanlist[tppairs[id].first];point p2 = tanlist[tppairs[id].second];
	    int dif = tppairs[id].second-tppairs[id].first-1;
	    double ds = 1.0/(double(dif)-1.0);
	    for(int i=0;i<dif;i++){
	      point pnew = p1 + (p2-p1)*ds*double(i);
	      tanlist[tppairs[id].first + i] = pnew;
	     }
	  };
	};
	FillDensityList();
	 getListTurnPts(tanlist);
	// first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	double actualWrithepl = writhepl*(1.0/(2.0*PI));
	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	writhepl = (1/(2*PI))*writhepl;
	// now get the integer winding
	std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
	std::vector<int> integerWinding =mw.getIntegerWindingsStarGen(list,tanlist,pairedTurnList);
	clearVecs();
	double writhepnl=0.0;
	for(int i=0;i<integerWinding.size();i++){
	  writhepnl = writhepnl + 2.0*double(integerWinding[i]); 
	}
	std::vector<double> output;
	output.push_back(actualWrithepl);
	output.push_back(writhepnl + (writhepl-actualWrithepl));
	output.push_back(writhepl+writhepnl);
	for(int i=0;i<integerWinding.size();i++){
	  output.push_back(integerWinding[i]); 
	}
	return output;
}

std::vector<double> localWrithe::getWritheGenClosed(std::vector<point>& list){
	rotateListToMaximum(list);
	derivList(list);
	//first remove any sharp tangent changes
	double prod;
	std::vector<std::pair<int,int> > tppairs;
	bool isFirst=true;
	std::pair<int,int> pr;
	for(int i=0;i<tanlist.size()-2;i++){
	  prod = (tanlist[i+2].getZ()-tanlist[i+1].getZ())*(tanlist[i+1].getZ()-tanlist[i].getZ());
	  if(isFirst){
	     if(prod<0){
	       isFirst=false;
	       pr.first=i+1;
	     }
	  }else{
	    //here we have foudn one tp
	    if(prod>0){
	      isFirst=true;
	      pr.second=i+2;
	      tppairs.push_back(pr);
	    }
	   }
	}
	for(int id=0;id<tppairs.size();id++){
	  if(tppairs[id].second-tppairs[id].first>2){
	    point p1 = tanlist[tppairs[id].first];point p2 = tanlist[tppairs[id].second];
	    int dif = tppairs[id].second-tppairs[id].first-1;
	    double ds = 1.0/(double(dif)-1.0);
	    for(int i=0;i<dif;i++){
	      point pnew = p1 + (p2-p1)*ds*double(i);
	      tanlist[tppairs[id].first + i] = pnew;
	     }
	  };
	};
	FillDensityList();
	reconstructCurve(list);
	// redo the turn pt fining
	dlist.clear();
	 tanlist.clear();
	 tderiv.clear();
	 densityList.clear();
	 turnlist.clear();
	 turnPts.clear();
	 arcLengthLst.clear();
	 derivList(list);
	 FillDensityList();
	 getListTurnPts(tanlist);
	// first add the main quaratures due to the indiviudal secttions separted by turning points 
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	double actualWrithepl = writhepl*(1.0/(2.0*PI));
	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	writhepl = (1/(2*PI))*writhepl;
	// now get the integer winding
	std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
	std::vector<int> integerWindings =mw.getIntegerWindingsClosed(list,tanlist,pairedTurnList);
	// get end anngle sign +1 or -1
	point tan1 = tanlist[0];
	point tan2 = tanlist[tanlist.size()-1];
	//std::cout<<tanlist.size()<<"\n";
	//tan1.printPoint();
	//tan2.printPoint();
	double turnAng,t;
	//std::cout<<tan1.getZ()<<" "<<tan2.getZ()<<"\n";
	t= -tan1.getZ()/(tan2.getZ()-tan1.getZ());
	//std::cout<<"t is "<<t<<"\n";
	if(t<-1||t>1){t=0.5;}
	double xnew = tan1.getX() + t*(tan2.getX()-tan1.getX());
	double ynew = tan1.getY() + t*(tan2.getY()-tan1.getY());
	turnAng = angle(-1.0*xnew,-1.0*ynew);
	//std::cout<<" turn angle "<<turnAng<<"\n";
	/*if(turnAng<-0.0001){
	  turnAng = 2.0*PI+turnAng;
	  }*/
	turnAng = -1.0*turnAng/(PI);;
	//std::cout<<"angle writhe "<<writhepl-actualWrithepl<<"\n";
	clearVecs();
	double writhepnl=0.0; 
	for(int i=0;i<integerWindings.size();i++){
	  writhepnl = writhepnl + 2.0*double(integerWindings[i]); 
	}
	// add the final turning angle
	writhepnl = writhepnl + turnAng;
	// add the final local writhe density
	std::vector<double> output;
	output.push_back(actualWrithepl);
	output.push_back(writhepnl + (writhepl-actualWrithepl));
	output.push_back(writhepl+writhepnl);
	for(int i=0;i<integerWindings.size();i++){
	  output.push_back(integerWindings[i]); 
	}
	return output;
}


std::vector<double> localWrithe::getWritheGenSmooth(std::vector<point>& list){
  derivList(list);
	FillDensityList();
	//smooth the curve
	//first remove any sharp tangent changes
	double prod;
	std::vector<std::pair<int,int> > tppairs;
	bool isFirst=true;
	std::pair<int,int> pr;
	for(int i=0;i<tanlist.size()-2;i++){
	  prod = (tanlist[i+2].getZ()-tanlist[i+1].getZ())*(tanlist[i+1].getZ()-tanlist[i].getZ());
	  if(isFirst){
	     if(prod<0){
	       isFirst=false;
	       pr.first=i+1;
	     }
	  }else{
	    //here we have foudn one tp
	    if(prod>0){
	      isFirst=true;
	      pr.second=i+2;
	      tppairs.push_back(pr);
	    }
	   }
	}
	for(int id=0;id<tppairs.size();id++){
	  if(tppairs[id].second-tppairs[id].first>2){
	    point p1 = tanlist[tppairs[id].first];point p2 = tanlist[tppairs[id].second];
	    int dif = tppairs[id].second-tppairs[id].first-1;
	    double ds = 1.0/(double(dif)-1.0);
	    for(int i=0;i<dif;i++){
	      point pnew = p1 + (p2-p1)*ds*double(i);
	      tanlist[tppairs[id].first + i] = pnew;
	     }
	  };
	};
	reconstructCurve(list);
	/*std::ofstream myfile;
	myfile.open("curveFiles/testCurveSmooth.dat");
        for(int i=0;i<list.size();i++){
	  myfile<<list[i].getX()<<" "<<list[i].getY()<<" "<<list[i].getZ()<<"\n";
	}
	myfile.close();*/
	// redo the turn pt fining
	dlist.clear();
	tanlist.clear();
	tderiv.clear();
	densityList.clear();
	turnlist.clear();
	turnPts.clear();
	arcLengthLst.clear();
	//search for any local zig-zag
	derivList(list);
	FillDensityList();
	getListTurnPts(tanlist);
	//check for points of inflection
	// first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	 double actualWrithepl = writhepl*(1.0/(2.0*PI));
	 // add any required interpolation arouns the turning points;
 	for(int j=1;j<turnlist.size()-1;j=j+2){
	   addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	 }
	writhepl = (1.0/(2.0*PI))*writhepl;
	 // now get the integer winding
	 std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
  std::vector<int> integerWinding =mw.getIntegerWindingsStarGen(list,tanlist,pairedTurnList);
	clearVecs();
	double writhepnl=0.0;
	for(int i=0;i<integerWinding.size();i++){
	  writhepnl = writhepnl + 2.0*double(integerWinding[i]); 
	}
	std::vector<double> output;
	output.push_back(actualWrithepl);
	output.push_back(writhepnl + (writhepl-actualWrithepl));
	output.push_back(writhepl+writhepnl);
	for(int i=0;i<integerWinding.size();i++){
	  output.push_back(integerWinding[i]); 
	}
	return output;
	return output;
}




std::vector<double> localWrithe::getWritheGenStandard(std::vector<point>& list){
	derivList(list);
	//first remove any sharp tangent changes
	double prod;
	std::vector<std::pair<int,int> > tppairs;
	bool isFirst=true;
	std::pair<int,int> pr;
	for(int i=0;i<tanlist.size()-2;i++){
	  prod = (tanlist[i+2].getZ()-tanlist[i+1].getZ())*(tanlist[i+1].getZ()-tanlist[i].getZ());
	  if(isFirst){
	     if(prod<0){
	       isFirst=false;
	       pr.first=i+1;
	     }
	  }else{
	    //here we have foudn one tp
	    if(prod>0){
	      isFirst=true;
	      pr.second=i+2;
	      tppairs.push_back(pr);
	    }
	   }
	}
	for(int id=0;id<tppairs.size();id++){
	  if(tppairs[id].second-tppairs[id].first>2){
	    point p1 = tanlist[tppairs[id].first];point p2 = tanlist[tppairs[id].second];
	    int dif = tppairs[id].second-tppairs[id].first-1;
	    double ds = 1.0/(double(dif)-1.0);
	    for(int i=0;i<dif;i++){
	      point pnew = p1 + (p2-p1)*ds*double(i);
	      tanlist[tppairs[id].first + i] = pnew;
	     }
	  };
	};
	FillDensityList();
	getListTurnPts(tanlist);
	// first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	double actualWrithepl = writhepl*(1.0/(2.0*PI));	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	writhepl = (1.0/(2.0*PI))*writhepl;
	// now get the integer winding
	std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
	std::vector<double> startAngs;
	std::vector<double> endAngs;
	std::vector<int> integerWindings =mw.getIntegerWindingsStarGenStandard(list,tanlist,pairedTurnList,startAngs,endAngs);
	clearVecs();
	double writhepnl=0.0;
	for(int i=0;i<integerWindings.size();i++){
	  writhepnl = writhepnl + 2.0*double(integerWindings[i]); 
	}
	double wrNoAngles = writhepl+writhepnl;
	// add the start angles
	for(int i=0;i<startAngs.size();i++){
	  writhepnl = writhepnl + (1/(PI))*startAngs[i]; 
	}
	// add the end angles
	for(int i=0;i<endAngs.size();i++){
	  writhepnl = writhepnl + (1/(PI))*endAngs[i]; 
	}
	std::vector<double> output;
	output.push_back(actualWrithepl);
	output.push_back(writhepnl + (writhepl-actualWrithepl));
	output.push_back(writhepl+writhepnl);
	for(int i=0;i<integerWindings.size();i++){
	  output.push_back(integerWindings[i]); 
	}
	return output;
}


std::vector<double> localWrithe::getWritheGenStandardSmooth(std::vector<point>& list){
	derivList(list);
	FillDensityList();
	//smooth the curve
	//first remove any sharp tangent changes
	double prod;
	std::vector<std::pair<int,int> > tppairs;
	bool isFirst=true;
	std::pair<int,int> pr;
	for(int i=0;i<tanlist.size()-2;i++){
	  prod = (tanlist[i+2].getZ()-tanlist[i+1].getZ())*(tanlist[i+1].getZ()-tanlist[i].getZ());
	  if(isFirst){
	     if(prod<0){
	       isFirst=false;
	       pr.first=i+1;
	     }
	  }else{
	    //here we have foudn one tp
	    if(prod>0){
	      isFirst=true;
	      pr.second=i+2;
	      tppairs.push_back(pr);
	    }
	   }
	}
	for(int id=0;id<tppairs.size();id++){
	  if(tppairs[id].second-tppairs[id].first>2){
	    point p1 = tanlist[tppairs[id].first];point p2 = tanlist[tppairs[id].second];
	    int dif = tppairs[id].second-tppairs[id].first-1;
	    double ds = 1.0/(double(dif)-1.0);
	    for(int i=0;i<dif;i++){
	      point pnew = p1 + (p2-p1)*ds*double(i);
	      tanlist[tppairs[id].first + i] = pnew;
	     }
	  };
	};
	reconstructCurve(list);
	// redo the turn pt fining
	/*std::ofstream ofile;
	ofile.open("testCurve.dat");
	for(int i=0;i<list.size();i++){
	  ofile<<list[i].getX()<<" "<<list[i].getY()<<" "<<list[i].getZ()<<"\n";
	}
	ofile.close();*/
	dlist.clear();
	tanlist.clear();
	tderiv.clear();
	densityList.clear();
	turnlist.clear();
	turnPts.clear();
	 arcLengthLst.clear();
	 derivList(list);
	 FillDensityList();
	 getListTurnPts(tanlist);
	 // first add the main quaratures due to the indiviudal secttions separted by turning points
	for(int i=0;i<turnlist.size()-1;i=i+2){
		double dummy = writhepl;
		SimpsonsRule(turnlist[i],turnlist[i+1]);
	}
	double actualWrithepl = writhepl*(1.0/(2.0*PI));
	// add any required interpolation arouns the turning points;
	for(int j=1;j<turnlist.size()-1;j=j+2){
		addInterpolatedTan(tanlist[turnlist[j]],densityList[turnlist[j]],tanlist[turnlist[j+1]],densityList[turnlist[j+1]],writhepl);
	}
	writhepl = (1.0/(2.0*PI))*writhepl;
	// now get the integer winding
	std::vector<std::pair<int,int> > pairedTurnList = makeTurnForNonLocal();
	mutualWind2 mw;
	std::vector<double> startAngs;
	std::vector<double> endAngs;
	std::vector<int> integerWindings =mw.getIntegerWindingsStarGenStandard(list,tanlist,pairedTurnList,startAngs,endAngs);
	clearVecs();
	double writhepnl=0.0;
	for(int i=0;i<integerWindings.size();i++){
	  writhepnl = writhepnl + 2.0*double(integerWindings[i]); 
	}
	// add the start angles
	double wrNoAngles = writhepl+writhepnl;
	for(int i=0;i<startAngs.size();i++){
	  //std::cout<<"start angles "<<startAngs[i]<<"\n";
	  writhepnl = writhepnl + (1.0/(PI))*startAngs[i]; 
	}
	// add the end angles
	for(int i=0;i<endAngs.size();i++){
	  writhepnl = writhepnl + (1.0/(PI))*endAngs[i]; 
	}
	std::vector<double> output;
	output.push_back(actualWrithepl);
	output.push_back(writhepnl + (writhepl-actualWrithepl));
	output.push_back(writhepl+writhepnl);
	for(int i=0;i<integerWindings.size();i++){
	  output.push_back(integerWindings[i]); 
	}
	return output;
}


double localWrithe::getEndAngles(std::vector<point> &curvelist1,std::vector<point> &curvelist2,std::vector<int> & turnList1,std::vector<int> & turnList2){
  //
  int curve1Up = curvelist1.size()-1;
  int curve2Up = curvelist2.size()-1;
  double ang2,ang1,angll,anguu,anglu,angul;
  double outsum = 0.0;
  double tanMin1 = curvelist1[1].getZ()-curvelist1[0].getZ();
  double tanMin2 = curvelist2[1].getZ()-curvelist2[0].getZ();
  double tanMax1 = curvelist1[curve1Up].getZ()-curvelist1[curve1Up-1].getZ();
  double tanMax2 = curvelist2[curve2Up].getZ()-curvelist2[curve2Up-1].getZ();
  binaryFind BS;
  double sigma1;
  //check if the start point of curve 1 is part of a downward moving section or
  if(tanMin1>=0 ){
    sigma1=1.0;
  }else{
    sigma1=-1.0;
  }
  double zv = curvelist1[0].getZ();
  for(int k=0;k<turnList2.size()-1;k++){
    double zmin = curvelist2[turnList2[k]].getZ();
    double zmax = curvelist2[turnList2[k+1]].getZ();
    //check if there is an angle
    if(zv >= zmin && zv <=zmax){
      // find the point
      std::vector<point> sublist(curvelist2.begin()+turnList2[k],curvelist2.begin()+turnList2[k+1]);
      int stIndex= 0;
      int endIndex = sublist.size()-1;
      std::pair<int,int> pair  = BS.getContainingPair(sublist,stIndex,endIndex,zv);
      // get the interpolated crossing
      point interpPoint =interp(sublist[pair.first],sublist[pair.second],zv);
      // get the angle
      double xdif = interpPoint.getX()-curvelist1[0].getX();
      double ydif = interpPoint.getY()-curvelist1[0].getY();
      double ang = angle(xdif,ydif);
      //std::cout<<"here angle "<<ang<<"\n";
      // now figure our the indicator function of the second curve sections
      double sigma2;
      if((sublist[1].getZ()-sublist[0].getZ())>=0.0){
	sigma2 = 1.0;
      }else{
	sigma2 = -1.0;
      }
      //check if this angle is the bottom part of a calculation or the top
      if(sigma1>0){
	//here it is the bottom
	outsum=outsum+(-1.0)*sigma1*sigma2*ang;
      }else{
	outsum = outsum +sigma1*sigma2*ang;
      }
    }
  }
  // now do start point 2
  //check if the start point of curve 1 is part of a downward moving section or
  if(tanMin2>=0 ){
    sigma1=1.0;
  }else{
    sigma1=-1.0;
  }
  zv = curvelist2[0].getZ();
  for(int k=0;k<turnList1.size()-1;k++){
    double zmin = curvelist1[turnList1[k]].getZ();
    double zmax = curvelist1[turnList1[k+1]].getZ();
    //check if there is an angle
    if(zv >= zmin && zv <=zmax){
      // find the point
      std::vector<point> sublist(curvelist1.begin()+turnList1[k],curvelist1.begin()+turnList1[k+1]);
      int stIndex= 0;
      int endIndex = sublist.size()-1;
      std::pair<int,int> pair  = BS.getContainingPair(sublist,stIndex,endIndex,zv);
      // get the interpolated crossing
      point interpPoint =interp(sublist[pair.first],sublist[pair.second],zv);
      // get the angle
      double xdif = interpPoint.getX()-curvelist2[0].getX();
      double ydif = interpPoint.getY()-curvelist2[0].getY();
      double ang = angle(xdif,ydif);
      //std::cout<<"here angle 2 "<<ang<<"\n";
      // now figure our the indicator function of the second curve sections
      double sigma2;
      if((sublist[1].getZ()-sublist[0].getZ())>=0.0){
	sigma2 = 1.0;
      }else{
	sigma2 = -1.0;
      }
      //check if this angle is the bottom part of a calculation or the top
      if(sigma1>0){
	//here it is the bottom
	outsum=outsum+(-1.0)*sigma1*sigma2*ang;
      }else{
	outsum = outsum +sigma1*sigma2*ang;
      }
    }
  }
  // now for the end angles
   // and the end angles
  int endCIndex1 = curvelist1.size()-1;
  if(curvelist1[endCIndex1].getZ() -curvelist1[endCIndex1-1].getZ()>=0 ){
    sigma1=1.0;
  }else{
    sigma1=-1.0;
  }
  zv = curvelist1[endCIndex1].getZ();
  for(int k=0;k<turnList2.size()-1;k++){
    double zmin = curvelist2[turnList2[k]].getZ();
    double zmax = curvelist2[turnList2[k+1]].getZ();
    if(zv >= zmin && zv <=zmax){
      // find the point
       std::vector<point> sublist(curvelist2.begin()+turnList2[k],curvelist2.begin()+turnList2[k+1]);
      int stIndex= 0;
      int endIndex = sublist.size()-1;
      std::pair<int,int> pair  = BS.getContainingPair(sublist,stIndex,endIndex,zv);
      // get the interpolated crossing
      point interpPoint =interp(sublist[pair.first],sublist[pair.second],zv);
      // get the angle
      double xdif = interpPoint.getX()-curvelist1[endCIndex1].getX();
      double ydif = interpPoint.getY()-curvelist1[endCIndex1].getY();
      double ang = angle(xdif,ydif);
      //std::cout<<"here 3 "<<ang<<"\n";
      double sigma2;
      if(sigma1>=0.0){
	sigma2 = 1.0;
      }else{
	sigma2 = -1.0;
      }
      //check if this angle is the bottom part of a calculation or the top
      if(sigma1<0){
	// here it is the bottom
	outsum = outsum + (-1.0)*sigma1*sigma2*ang;
      }else{
	outsum = outsum +sigma1*sigma2*ang;
      }
    }
  }
   int endCIndex2 = curvelist2.size()-1;
  if(curvelist2[endCIndex2].getZ() -curvelist2[endCIndex2-1].getZ()>=0 ){
    sigma1=1.0;
  }else{
    sigma1=-1.0;
  }
  for(int k=0;k<turnList1.size()-1;k++){
    double zmin = curvelist1[turnList1[k]].getZ();
    double zmax = curvelist1[turnList1[k+1]].getZ();
    if(curvelist2[endCIndex2].getZ() >= zmin && curvelist2[endCIndex2].getZ() <=zmax){
      // find the point
       std::vector<point> sublist(curvelist1.begin()+turnList1[k],curvelist1.begin()+turnList1[k+1]);
      int stIndex= 0;
      int endIndex = sublist.size()-1;
      double zv = curvelist2[endCIndex2].getZ();
      std::pair<int,int> pair  = BS.getContainingPair(sublist,stIndex,endIndex,zv);
      // get the interpolated crossing
      point interpPoint =interp(sublist[pair.first],sublist[pair.second],zv);
      // get the angle
      double xdif = interpPoint.getX()-curvelist2[endCIndex2].getX();
      double ydif = interpPoint.getY()-curvelist2[endCIndex2].getY();
      double ang = angle(xdif,ydif);
      //std::cout<<"here 4 "<<ang<<"\n";
      double sigma2;
      if(sigma1>=0.0){
	sigma2 = 1.0;
      }else{
	sigma2 = -1.0;
      }
      //check if this angle is the bottom part of a calculation or the top
      if(sigma1<0){
	// here it is the bottom
	outsum = outsum + (-1.0)*sigma1*sigma2*ang;
      }else{
	outsum = outsum +sigma1*sigma2*ang;
      }
    }
  }
return outsum*(1.0/(2.0*PI));
}

/*double localWrithe::getWinding(std::vector<point>& list1,std::vector<point>& list2){
        std::vector<point> tangent1 = derivListReturn(list1);
	std::vector<int> turnList1 = getListTurnPtsReturn(tangent1);
	std::vector<point> tangent2 = derivListReturn(list2);
	std::vector<int> turnList2 = getListTurnPtsReturn(tangent2);
	std::vector<std::pair<int,int> > pairedTurnList1 = makeTurnForNonLocal(turnList1);
	std::vector<std::pair<int,int> > pairedTurnList2 = makeTurnForNonLocal(turnList2);
	// first calculate the integer part
	std::cout<<pairedTurnList1.size()<<"\n";
	std::cout<<pairedTurnList2.size()<<"\n";
	for(int i=0;i<pairedTurnList1.size();i++){
	  std::cout<<pairedTurnList1[i].first<<" "<<pairedTurnList1[i].second<<"\n";
	}
	for(int i=0;i<pairedTurnList2.size();i++){
	  std::cout<<pairedTurnList2[i].first<<" "<<pairedTurnList2[i].second<<"\n";
	}
	mutualWind2 mw;
	int integerWindings=0;
	//integerWindings = mw.getWindingInteger(list1,list2,tangent1,tangent2,pairedTurnList1,pairedTurnList2);
	// next the end angle contributions, several cases
	double endSum = getEndAngles(list1,list2,tangent1,tangent2);
	clearVecs();
	return endSum + double(integerWindings); 
	}*/

double localWrithe::getWinding(std::vector<point>& list1,std::vector<point>& list2){
  std::vector<int> turnList1 = getListTurnPtsReturn(list1);
  std::vector<int> turnList2 = getListTurnPtsReturn(list2);
  // first calculate the integer part
  mutualWind2 mw;
  int integerWindings=0;
  integerWindings = mw.getWindingInteger(list1,list2,turnList1,turnList2);
  // next the end angle contributions, several cases
  double endSum = getEndAngles(list1,list2,turnList1,turnList2);
  clearVecs();
  std::cout<<endSum<<" "<<double(integerWindings)<<"\n";
  return endSum + double(integerWindings); 
}




point localWrithe::dbold(std::vector<point>& pointList,int size,int i,int j){
  i = i%(size-1);
  j = j%(size-1);
  return pointList[i+1].dif(pointList[j+1]);
}

point localWrithe::dboldLk(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i,int j){
  i = i%(size-1);
  j = j%(size2-1);
  return pointList[i+1].dif(pointList2[j+1]);
}

point localWrithe::te(std::vector<point>& pointList,int size,int i){
  point temp = dbold(pointList,size,i,i-1);
  temp.normalise();
  return temp;
}



double localWrithe::acosC(double temp){
  if(temp < -1.0){
    return PI;
  }else if(temp > 1.0){
    return 0;
  }else{
    return acos(temp);
  }
}

double localWrithe::mu(std::vector<point>& pointList,int size,int i,int k,int m, int j){
  point ti= te(pointList,size,i);
  point tj=te(pointList,size,j);
  point dkm = dbold(pointList,size,k,m);
  point un = ti.cross(dkm);
  point deux = dkm.cross(tj);
  if( !un.isNonzero() || !deux.isNonzero()){
    return 0;
  }else{
    un.normalise();
    deux.normalise();
    double temp = un.dotprod(deux);
    point cp = ti.cross(tj);
    double signProd = dkm.dotprod(cp); 
    if(std::abs(signProd)<0.000000001){
      return 0.0;
    }else{
      if(signProd>= 0){
	return acosC(temp);
      }else{
	return (-1.0)*acosC(temp);
      }
    }
   }
}

double localWrithe::muLink(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i, int j){
  int iplus = i+1;
  int iminus = i-1;
  if(i==0){
    iminus=size-1;
  }
  if(i==size-1){
    iplus=0;
  }
  int jplus = j+1;
  int jminus = j-1;
  if(j==0){
    jminus=size2-1;
  }
  if(j==size2-1){
    jplus=0;
  }
  point tiplus= pointList[iplus].dif(pointList[i]);
  point timinus= pointList[i].dif(pointList[iminus]);
  point tjplus= pointList2[jplus].dif(pointList2[j]);
  point tjminus= pointList2[j].dif(pointList2[jminus]);
  double ds1 =0.5*(tiplus.length()+timinus.length());
  double ds2 =0.5*(tjplus.length()+tjminus.length());
  timinus.normalise();
  tjminus.normalise();
  tiplus.normalise();
  tjplus.normalise();
  point ti = (timinus+tiplus)*0.5;
  point tj = (tjminus+tjplus)*0.5;
  point dkm =pointList[i].dif(pointList2[j]);
  double dkmlen = dkm.length();
  double dkmCubed =dkmlen*dkmlen*dkmlen;
  point un = ti.cross(tj);
  return(ds1*ds2*un.dotprod(dkm)/dkmCubed);
}




double localWrithe::wij(std::vector<point>& pointList,int size,int i,int j){
  double temp =0;
  if(j>i){
    double t1 = mu(pointList,size,i,i-1,j-1,j);
    double t2 = mu(pointList,size,i,i,j-1,j);
    double t3 = mu(pointList,size,i,i-1,j,j);
    double t4 = mu(pointList,size,i,i,j,j);
    temp = t1-t2-t3+t4;
  }
  return temp;
}


double localWrithe::DI(std::vector<point>& pointList){
  //pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +wij(pointList,listSize,i,j);
     }
  }
  return sigsum/(2*PI);
}

double localWrithe::DIAbs(std::vector<point>& pointList){
  //pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +std::abs(wij(pointList,listSize,i,j));
     }
  }
  return sigsum/(2*PI);
}

/*
std::vector<double> localWrithe::DI(std::vector<point>& pointList){
  //pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  std::vector<double> writhes;
  double sigsum=0.0;
  double sigsumLoc =0.0;
  double sigsumNonLoc =0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +wij(pointList,listSize,i,j);
     }
  }
  sigsum = sigsum/(2*PI);
  ...
  writhes.push_back(sigsum);
  /// 
  return writhes;
}
*/

double localWrithe::DIClosed(std::vector<point>& pointList){
  pointList.push_back(pointList[0]);
  int listSize = pointList.size();
  double sigsum=0.0;
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize-1;j++){
      sigsum = sigsum +wij(pointList,listSize,i,j);
     }
  }
  return sigsum/(2*PI);
}





double localWrithe::DIClosedLk(std::vector<point>& pointList,std::vector<point>& pointList2){
  int listSize = pointList.size();
  int listSize2 = pointList2.size();
  double sigsum=0.0;
  for(int i=0;i<listSize;i++){
    for(int j=0;j<listSize2;j++){
      sigsum = sigsum +muLink(pointList,pointList2,listSize,listSize2,i,j);
     }
  }
  return sigsum/(4*PI);
}


std::vector<std::pair<std::pair<int,int>,double> > localWrithe::DIDownSample(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DI(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         std::pair<std::pair<int,int>,double> fppoint;
         fppoint.first=pr;fppoint.second=wrval;
         fingerPrintList.push_back(fppoint);
       }
    }   
  }
  return fingerPrintList;
}


std::vector<std::pair<std::pair<int,int>,double> > localWrithe::DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DIAbs(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         std::pair<std::pair<int,int>,double> fppoint;
         fppoint.first=pr;fppoint.second=wrval;
         fingerPrintList.push_back(fppoint);
       }
    }   
  }
  return fingerPrintList;
}



void  localWrithe::DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  std::ofstream fingerFile;
  fingerFile.open(filename);
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DI(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         fingerFile<<k<<" "<<k+l<<" "<<wrval<<"\n";
       }
    }   
  }
  fingerFile.close();
}

void  localWrithe::DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  std::vector<std::pair<std::pair<int,int>,double> > fingerPrintList;
  std::ofstream fingerFile;
  fingerFile.open(filename);
  for(int l=5;l<mol.size();l++){
     for(int k=0;k<mol.size()-5;k++){
       if(k+l<mol.size()){
	 // grab subsection of curve;
	 std::vector<point>::const_iterator first = mol.begin()+k;
         std::vector<point>::const_iterator last = mol.begin()+k+l;
	 std::vector<point> subVec(first, last);
	 double wrval = DIAbs(subVec);
         std::pair<int,int> pr;
         pr.first=k;pr.second=k+l;     
         fingerFile<<k<<" "<<k+l<<" "<<wrval<<"\n";
       }
    }   
  }
  fingerFile.close();
}
