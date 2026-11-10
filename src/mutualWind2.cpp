#include "mutualWind2.h"

windingStorage::windingStorage(std::vector<point> &sublist1In,std::vector<point> &sublist2In,bool &getWindIn,int &sigma1In,int &sigma2In){
  sublist1 = sublist1In;
  sublist2 = sublist2In;
  getWind = getWindIn;
  sigma1 = sigma1In;
  sigma2 = sigma2In;
}

std::vector<point> windingStorage::getSublist1(){
  return sublist1;
};

std::vector<point> windingStorage::getSublist2(){
  return sublist2;
}

bool windingStorage::getIsWind(){
  return getWind;
};

int windingStorage::getSigma1(){
  return sigma1;
}

int windingStorage::getSigma2(){
  return sigma2;
}



mutualWind2::mutualWind2(){
}

point mutualWind2::interp(point& lower,point& upper,double zval){
  double t;
  if(std::abs(upper.getZ() -lower.getZ()) > 0.0001){
  t = (zval -lower.getZ())/(upper.getZ() -lower.getZ());
  }else{
    t = 0.5;
  }
  if(std::abs(t)>1.0){
    std::cout<<"here t is "<<t<<"\n";
    std::cout<<lower.getZ()<<" "<<upper.getZ()<<"  "<<zval<<"\n";
  }
  double xval = lower.getX() + t*(upper.getX()-lower.getX());
  double yval = lower.getY() + t*(upper.getY()-lower.getY());
  point p(xval,yval,zval);
  return p;
}

point mutualWind2::getInterpolatedTop(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList){
  double tanz1 = tanlist[turnMiddle.first].getZ();
  double tanz2 = tanlist[turnMiddle.second].getZ();
  double mids;
  if(std::abs(tanz2-tanz1)>0.0001){
    mids = -tanz1/(tanz2-tanz1);
  }else{
    mids=0.5;
  }
  //std::cout<<"interpolation bad ? "<<mids<<" "<<tanz1<<" "<<tanz2<<"\n";
  double x = curveList[turnMiddle.first].getX() + mids*(curveList[turnMiddle.second].getX()-curveList[turnMiddle.first].getX());
  double y = curveList[turnMiddle.first].getY() + mids*(curveList[turnMiddle.second].getY()-curveList[turnMiddle.first].getY());
  double z = curveList[turnMiddle.first].getZ() + mids*(curveList[turnMiddle.second].getZ()-curveList[turnMiddle.first].getZ());
  point out(x,y,z);
  return out;
}

std::pair<double,double> mutualWind2::getExtremaZ(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList){
  double tanz1 = tanlist[turnMiddle.first].getZ();
  double tanz2 = tanlist[turnMiddle.second].getZ();
  double mids;
  if(std::abs(tanz2-tanz1)>0.0001){
    mids = -tanz1/(tanz2-tanz1);
  }else{
    mids=0.5;
  }
  double z = curveList[turnMiddle.first].getZ() + mids*(curveList[turnMiddle.second].getZ()-curveList[turnMiddle.first].getZ());
  std::pair<double,double> extrema;
  extrema.first = std::max(std::max(curveList[turnMiddle.first].getZ(),curveList[turnMiddle.second].getZ()),z);
  extrema.second = std::min(std::min(curveList[turnMiddle.first].getZ(),curveList[turnMiddle.second].getZ()),z);
  return extrema;
}


int mutualWind2::quadrant(point &pointdif){
	double xdif = pointdif.getX();
	double ydif = pointdif.getY();
	if(xdif< 0){
		if(ydif<0){
			return 3;
		}else{
			return 2;
		}	
	}else{
		if(ydif<0){
			return 4;
		}else{
			return 1;
		}
	}
}


point  mutualWind2::getTurningAngleQuad(std::pair<int,int> &turnMiddle,std::vector<point> &tanlist,std::vector<point> &curveList){
  double tanz1 = tanlist[turnMiddle.first].getZ();
  double tanz2 = tanlist[turnMiddle.second].getZ();
  double mids;
  if(std::abs(tanz2-tanz1)>0.0001){
    mids = -tanz1/(tanz2-tanz1);
  }else{
    mids=0.5;
  }
  
  double x = tanlist[turnMiddle.first].getX() + mids*(tanlist[turnMiddle.second].getX()-tanlist[turnMiddle.first].getX());
  double y = tanlist[turnMiddle.first].getY() + mids*(tanlist[turnMiddle.second].getY()-tanlist[turnMiddle.first].getY());
  double z = tanlist[turnMiddle.first].getZ() + mids*(tanlist[turnMiddle.second].getZ()-tanlist[turnMiddle.first].getZ());
  point out(x,y,z);
  return out;
}


void mutualWind2::branchTracker(int& currBranch,int& prevBranch,int& windingSum,double &rotDirec){
  if(std::abs(prevBranch-currBranch)>1){
    /*If we are here two quadrants have been jumped, we must check teh rotation direction*/
    if(prevBranch==2  && currBranch>2){
      if(rotDirec< 0){
	windingSum = windingSum +1;
      }
    }else if(prevBranch==3 && currBranch<3){
      if(rotDirec> 0){
	windingSum = windingSum -1;
      }
    }// not too forgte 1 to 3 
    else if(prevBranch==1  && currBranch==3){
      if(rotDirec< 0){
	windingSum = windingSum +1;
      }
    }// or 4 to 2
    else if(prevBranch==4  && currBranch==2){
      if(rotDirec> 0){
	windingSum = windingSum -1;
      }
    }
  }else{
    //std::cout<<"in brach tracker"<<prevBranch<<" "<<currBranch<<"\n";
    if(prevBranch ==2 && currBranch == 3){
      windingSum = windingSum +1; 
    }else if(prevBranch ==3 && currBranch == 2){
      windingSum = windingSum -1; 
    }
  }
}

double mutualWind2::angle(double xdif,double ydif){
	if(xdif==0 && ydif ==0){
	return 0;	
	}
	std::complex<double> imag(xdif,ydif);
	double angv = std::arg(imag);
	return angv;
}

/*The following function is used to check the sections share exaclty the same mutal z range*/

bool mutualWind2::checkAlligned(std::vector<point> &sublist1,std::vector<point> &sublist2){
  if(std::abs(sublist1[0].getZ()-sublist2[0].getZ())<0.00000000001 && std::abs(sublist1[sublist1.size()-1].getZ()-sublist2[sublist2.size()-1].getZ())<0.0000000001){
    return true;
  }else{
    return false;
  }
}

/*getNext Index find the next winding vector (possibly using interpolation) and gets its quadrant. It also update the positions of the interator on each curve*/

std::pair<int,point> mutualWind2:: getNextIndex(std::vector<point> &sublist1,std::vector<point> &sublist2,std::pair<int,int> &indexPair,int &prevQuad){
  /*First check if either point is the end of one of the sections, we have to treat these specially */
  std::pair<int,point> newQuadJoinvec;
  bool isUpper1= false;
  if(indexPair.first == sublist1.size()-1){
    isUpper1 =true;
  }
  bool isUpper2 = false;
  if(indexPair.second == sublist2.size()-1){
    isUpper2 =true;
  }
  //std::cout<<"met either end "<<isUpper1<<" "<<isUpper2<<"\n";
  int newQuad;
  point newPoint;
  point joinVec;
  if(isUpper1 ==true && isUpper2 == true){
    newQuadJoinvec.first = prevQuad;
    std::cout<<"Bug: both indicies of getNextIndex are at their lists ends but still getNextINdex has been called\n";
  }
  else if(isUpper1==true && isUpper2 ==false){
    //std::cout<<"should be here "<<"\n";
    if(indexPair.second ==  sublist2.size()-2){
      indexPair.second = indexPair.second+1;
      joinVec = sublist2[indexPair.second].dif(sublist1[indexPair.first]); 
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec;
    }
    else{
      indexPair.second = indexPair.second+1;
      newPoint = interp(sublist1[indexPair.first-1],sublist1[indexPair.first],sublist2[indexPair.second].getZ());
      joinVec = sublist2[indexPair.second].dif(newPoint);
      newQuadJoinvec.first = quadrant(joinVec);
       newQuadJoinvec.second = joinVec;
    }
  }
  else if(isUpper1== false && isUpper2 ==true){
    if(indexPair.first ==  sublist1.size()-2){
      indexPair.first = indexPair.first+1;
      joinVec = sublist2[indexPair.second].dif(sublist1[indexPair.first]); 
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec;
    }
    else{
      indexPair.first = indexPair.first+1;
      newPoint = interp(sublist2[indexPair.second-1],sublist2[indexPair.second],sublist1[indexPair.first].getZ());
      joinVec = newPoint.dif(sublist1[indexPair.first]);
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec;
    }
  }
  else{
    /*here neither indices are at the end of either curve need to know which point is above the other or if they are equal*/
    //std::cout<<sublist1[indexPair.first+1].getZ()<<" "<<sublist2[indexPair.second+1].getZ()<<" in routine\n";
    if(sublist1[indexPair.first+1].getZ()<sublist2[indexPair.second+1].getZ()){
      // std::cout<<"here1"<<"\n";
      /*the point on sublist2 is above that of 1, interpolate on sublist 2*/
      indexPair.first = indexPair.first+1;
      newPoint = interp(sublist2[indexPair.second],sublist2[indexPair.second+1],sublist1[indexPair.first].getZ());
      joinVec = newPoint.dif(sublist1[indexPair.first]);
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec;
    }else if(sublist1[indexPair.first+1].getZ()>sublist2[indexPair.second+1].getZ()){
      //std::cout<<"here2"<<"\n";
       /*the point on sublist1 is above that of 2, interpolate on sublist 1*/
      indexPair.second = indexPair.second+1;
      newPoint = interp(sublist1[indexPair.first],sublist1[indexPair.first+1],sublist2[indexPair.second].getZ());
      joinVec = sublist2[indexPair.second].dif(newPoint);
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec;
    }else{
      // std::cout<<"here3"<<"\n";
      /*Here the points happend to be at the same height*/
      indexPair.first = indexPair.first+1;
      indexPair.second = indexPair.second+1;
      joinVec = sublist2[indexPair.second].dif(sublist1[indexPair.first]); 
      newQuadJoinvec.first = quadrant(joinVec);
      newQuadJoinvec.second = joinVec; 
      /*std::cout<<indexPair.first<<"\n";
       std::cout<<indexPair.second<<"\n";
       std::cout<<joinVec.getX()<<" "<<joinVec.getY()<<" "<<joinVec.getZ()<<"\n";*/
    }
  }
  return  newQuadJoinvec;
}

/*The following function takes two lists which share a mutual point an extends them so that their end points are given by the
same tangent interpolation we use to idenify the actual turning point of the curves. This case is special as two poinst on thw two curves are actaull the same, so we have a local minimum or maximum*/

void mutualWind2::AppendExtremumJoin(std::vector<point> &sublist1,std::vector<point> &sublist2,point &newS1minExtremum,point &newS2maxExtremum,point &newJointExtremum){
  /*first determine when we have a local minimum or maximum*/
  if(sublist1[0].getZ()< sublist1[sublist1.size()-1].getZ()){
    /*we have a local maximum*/
    if(newJointExtremum.getZ()>sublist1[sublist1.size()-1].getZ()){
      sublist1.push_back(newJointExtremum);
    }else{
      sublist1[sublist1.size()-1] = newJointExtremum;
    }
    if(newJointExtremum.getZ()>sublist2[0].getZ()){
      sublist2.insert(sublist2.begin(),newJointExtremum);
    }else{
      sublist2[0] = newJointExtremum;
    }
    // also add the non turning point extrema
    if(newS1minExtremum.getZ()<sublist1[0].getZ()){
      sublist1.insert(sublist1.begin(),newS1minExtremum);
    }else{
      sublist1[0] = newS1minExtremum;
    }
    if(newS2maxExtremum.getZ()<sublist2[sublist2.size()-1].getZ()){
      sublist2.push_back(newS2maxExtremum);
    }else{
      sublist2[sublist2.size()-1] = newS2maxExtremum;
    }
  }else{
    /*we have a local minimum*/
    if(newJointExtremum.getZ()<sublist1[sublist1.size()-1].getZ()){
      sublist1.push_back(newJointExtremum);
    }else{
      sublist1[sublist1.size()-1] = newJointExtremum;
    }
    if(newJointExtremum.getZ()<sublist2[0].getZ()){
      sublist2.insert(sublist2.begin(),newJointExtremum);
    }else{
      sublist2[0] = newJointExtremum;
    }
    // also add the non turning point extrema
    if(newS1minExtremum.getZ()>sublist1[0].getZ()){
      sublist1.insert(sublist1.begin(),newS1minExtremum);
    }else{
      sublist1[0] = newS1minExtremum;
    }
    if(newS2maxExtremum.getZ()>sublist2[sublist2.size()-1].getZ()){
      sublist2.push_back(newS2maxExtremum);
    }else{
      sublist2[sublist2.size()-1] = newS2maxExtremum;
    }
  }
}



void mutualWind2::AppendExtremumJoinClosed(std::vector<point> &sublist1,std::vector<point> &sublist2,point &newS1minExtremum,point &newS2maxExtremum,point &newJointExtremum){
  /*first determine when we have a local minimum or maximum*/
    /*we have a local maximum*/
  if(newJointExtremum.getZ()>sublist1[0].getZ()){
    sublist1.insert(sublist1.begin(),newJointExtremum);
  }else{
    sublist1[0] = newJointExtremum;
  }
  if(newJointExtremum.getZ()>sublist2[sublist2.size()-1].getZ()){
    sublist2[sublist2.size()-1] = newJointExtremum;
  }
  // also add the non turning point extrema
  if(newS1minExtremum.getZ()<sublist1[0].getZ()){
    sublist1.insert(sublist1.begin(),newS1minExtremum);
  }else{
    sublist1[0] = newS1minExtremum;
  }
  if(newS2maxExtremum.getZ()<sublist2[sublist2.size()-1].getZ()){
    sublist2[sublist2.size()-1] = newS2maxExtremum;
  }
}




/*The following function takes two lists extends them so that their end points are given by the
same tangent interpolation we use to idenify the actual turning point of the curves.*/

void mutualWind2::checkAppendExtrema(std::vector<point> &sublist,point &minExtremum,point &maxExtremum){
  /*First check if the curve is moving up or down*/
  if(sublist[0].getZ()<sublist[sublist.size()-1].getZ()){
    /*Upward moving section*/
    if(minExtremum.getZ()<sublist[0].getZ()){
      sublist.insert(sublist.begin(),minExtremum);
    }else{
      sublist[0] =minExtremum; 
    }
    if(maxExtremum.getZ()>sublist[sublist.size()-1].getZ()){
      sublist.push_back(maxExtremum);
    }else{
      sublist[sublist.size()-1]= maxExtremum;
    }
  }else{
    if(minExtremum.getZ()>sublist[0].getZ()){
      sublist.insert(sublist.begin(),minExtremum);
    }else{
      sublist[0] =minExtremum; 
    }
    if(maxExtremum.getZ()<sublist[sublist.size()-1].getZ()){
      sublist.push_back(maxExtremum);
    }else{
      sublist[sublist.size()-1]= maxExtremum;
    }
  }
}


std::vector<std::vector<point> > mutualWind2::interpolateForWinding(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > turningList){
  point newS1minExtremum;
  point newS1maxExtremum;
  std::vector<std::vector<point> > sectionList;
  for(int i=0;i<turningList.size()-1;i++){
    if(i==0){
      // find interpolated extrema
      if(turningList[i].first != turningList[i].second){
	newS1minExtremum = getInterpolatedTop(turningList[i],tanList,curveList);
      }else{
	newS1minExtremum = curveList[turningList[i].first];
      }
      if(turningList[i+1].first != turningList[i+1].second){
	newS1maxExtremum = getInterpolatedTop(turningList[i+1],tanList,curveList);
      }else{
	newS1maxExtremum = curveList[turningList[i+1].second];
      }
      std::vector<point> sublist1(curveList.begin()+turningList[i].second,curveList.begin()+turningList[i+1].first +1);
      checkAppendExtrema(sublist1,newS1minExtremum,newS1maxExtremum);
      sectionList.push_back(sublist1);
      //std::cout<<"pre here "<<sublist1.front().getZ()<<" "<<sublist1.back().getZ()<<"\n";
    }else{
      if(turningList[i+1].first != turningList[i+1].second){
	newS1maxExtremum = getInterpolatedTop(turningList[i+1],tanList,curveList);
      }else{
	newS1maxExtremum = curveList[turningList[i+1].second];
      }
      std::vector<point> sublist1(curveList.begin()+turningList[i].second,curveList.begin()+turningList[i+1].first +1);
      checkAppendExtrema(sublist1,sectionList[i-1].back(),newS1maxExtremum);
      sectionList.push_back(sublist1);
      //std::cout<<"here "<<sublist1.front().getZ()<<" "<<sublist1.back().getZ()<<"\n";
    }
  }
  return sectionList;
}


std::vector<windingStorage>  mutualWind2::getMutualSections(std::vector<std::vector<point> > &superSublist1,std::vector<std::vector<point> > &superSublist2,int &firstIndex){
  std::vector<point> sublist2Prev;
  std::vector<point> sublist1;
  std::vector<point> sublist2;
  std::vector<windingStorage> windingInfo;
  int sigma1,sigma2;
  for(int j=0;j<superSublist2.size();j++){
    bool isTurningPoint = false;
    if(j>0){
      sublist2Prev = sublist2;
    }
    sublist1 = superSublist1[firstIndex];
    sublist2 = superSublist2[j];
    sigma1 =1;
    sigma2 =1;
    int s1Up = sublist1.size()-1;
    int s2Up = sublist2.size()-1;
    // we want both lists to be ordered so them move up in z, check if they need to an reverse if necessar
    if(sublist1[s1Up].getZ()<sublist1[0].getZ()){
      std::reverse(sublist1.begin(),sublist1.end());
      sigma1 = -1;
    }
    if(sublist2[s2Up].getZ()<sublist2[0].getZ()){
      std::reverse(sublist2.begin(),sublist2.end());
      sigma2 = -1;
    }
    double z1min = sublist1[0].getZ();
    double z2min = sublist2[0].getZ();
    double z1max = sublist1[s1Up].getZ();
    double z2max = sublist2[s2Up].getZ();
    //std::cout<<"z values post chop non joined "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n"; 
    //std::cout<<j<<" sigma 2 is "<<sigma2<<"\n";
    bool getWind;
    binaryFind BS;
    /*Now chop the domains so that they mutually ovelrap*/
    if(z1max<= z2min || z2max<= z1min){
      // no overlap of the domains
      getWind=false;
      std::cout<<"fail ?"<<"\n";
    }else if(z1min<=z2min && z1max>= z2max){
      // [z2min,z2max] \subset [z1min,z1max]
      int low = 0;
      std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
      std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
      point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
      point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
      std::vector<point> sublistFin1(sublist1.begin()+pair1.second,sublist1.begin()+ pair2.first+1);
      sublistFin1.insert(sublistFin1.begin(),newLower);
      sublistFin1.push_back(newUpper);
      sublist1 = sublistFin1;
      getWind=true;
      //std::cout<<"here 1\n";
      isTurningPoint = true;
    }else if(z2min<=z1min && z2max>= z1max){
      // [z1min,z1max] \subset [z2min,z2max]
      int low = 0;
      std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
      std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
      point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
      point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
      std::vector<point> sublistFin2(sublist2.begin()+ pair1.second,sublist2.begin() +pair2.first+1);
      sublistFin2.insert(sublistFin2.begin(),newLower);
      sublistFin2.push_back(newUpper);
      sublist2 = sublistFin2;
      getWind=true;
    }
    else if(z2min<=z1max && z2min>= z1min){
      // neither section contians the other but z1max is below z2max
      int low = 0;
      std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
      std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
      point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
      point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
      std::vector<point> sublistFin1(sublist1.begin()+ pair1.second, sublist1.begin() + s1Up+1);
      std::vector<point> sublistFin2(sublist2.begin(),sublist2.begin()+pair2.first+1);
      sublistFin1.insert(sublistFin1.begin(),newLower);
      sublistFin2.push_back(newUpper);
      sublist2 = sublistFin2;
      sublist1= sublistFin1;
      getWind=true;
      if(sigma2 == 1){
	//std::cout<<"here 2\n";
	isTurningPoint = true;
      }
    }
    else if(z1min<=z2max && z1min>= z2min){
      // neither section contians the other but z1max is below z2max
      int low = 0;
      std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
      std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
      point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
      point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
      std::vector<point> sublistFin2(sublist2.begin()+pair1.second,sublist2.begin()+s2Up+1);
      std::vector<point> sublistFin1(sublist1.begin(),sublist1.begin()+pair2.first+1);
      sublistFin2.insert(sublistFin2.begin(),newLower);
      sublistFin1.push_back(newUpper);
      sublist2 = sublistFin2;
      sublist1= sublistFin1;
      getWind=true;
      if(sigma2 == -1){
	//std::cout<<"here 3\n";
	isTurningPoint = true;
      }
    }
    // now check if we need to do any cleaning
    if(j>0){
      if(sigma2==-1){
	if(isTurningPoint){
	  sublist2[sublist2.size()-1] = sublist2Prev.back();
	}
      }else{
	if(isTurningPoint){
	  sublist2[0] = sublist2Prev.front();
	}
      }
    }
    //std::cout<<"re check chop  "<<sublist1.front().getZ()<<" "<<sublist2.front().getZ()<<" "<<sublist1.back().getZ()<<" "<<sublist2.back().getZ()<<"\n";
    windingStorage w(sublist1,sublist2,getWind,sigma1,sigma2);
    windingInfo.push_back(w);
  }
  return windingInfo;
}

/*now the actual routine to calculate the integer windings of a pair of secions sharinga  local extremum*/


int mutualWind2::getIntegerWindingOfPairJoined(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp){
  point newS1minExtremum;
  point newJointExtremum;
  point newS2maxExtremum;
  std::pair<int,int> turnPair;
  /*find the interpolated extrema*/
  //std::cout<<s1minInp.first<<" "<<s1minInp.second<<" "<<s2minInp.first<<" "<<s2minInp.second<<"\n";
  //std::cout<<curveList[s1minInp.first].getZ()<<" "<<curveList[s1minInp.second].getZ()<<" "<<curveList[s2minInp.first].getZ()<<" "<<curveList[s2minInp.second].getZ()<<"\n";
  if(s1minInp.first != s1minInp.second){
    newS1minExtremum = getInterpolatedTop(s1minInp,tanList,curveList);
  }else{
    newS1minExtremum = curveList[s1minInp.first];
  }
  if(s2maxInp.first != s2maxInp.second){
    newS2maxExtremum = getInterpolatedTop(s2maxInp,tanList,curveList);
  }else{
    newS2maxExtremum = curveList[s2maxInp.first];
  }
  /*The turning point*/
  if(s1maxInp.first != s1maxInp.second){
    newJointExtremum = getInterpolatedTop(s1maxInp,tanList,curveList);
  }else{
    newJointExtremum =  curveList[s1maxInp.first];
  }
  //get the relevant parts of the list
  std::vector<point> sublist1(curveList.begin()+s1minInp.second,curveList.begin()+s1maxInp.first +1);
  std::vector<point> sublist2(curveList.begin()+s2minInp.second,curveList.begin()+s2maxInp.first +1);
  //alter to account for our interpolated extrema
  AppendExtremumJoin(sublist1,sublist2,newS1minExtremum,newS2maxExtremum,newJointExtremum);
  int s1Up = sublist1.size()-1;
  int s2Up = sublist2.size()-1;
  // we want both lists to be ordered so them move up in z, check if they need to an reverse if necessary
  bool isRev1 =false;
  int  sigma1=1;
  //check the sections have a sufficinet size
  if(std::abs(sublist1[s1Up].getZ()-sublist1[0].getZ())>0.0000001 && std::abs(sublist2[s2Up].getZ()-sublist2[0].getZ())>0.0000001){
  if(sublist1[s1Up].getZ()<sublist1[0].getZ()){
    std::reverse(sublist1.begin(),sublist1.end());
    isRev1 = true;
    sigma1 = -1;
  }
  bool isRev2 =false;
  int  sigma2=1;
  if(sublist2[s2Up].getZ()<sublist2[0].getZ()){
    std::reverse(sublist2.begin(),sublist2.end());
    isRev2 = true;
    sigma2 = -1;
  }
  double z1min = sublist1[0].getZ();
  double z2min = sublist2[0].getZ();
  double z1max = sublist1[s1Up].getZ();
  double z2max = sublist2[s2Up].getZ();
  //std::cout<<"z values post chop joined "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n";
  //check for an infecltion point
  //std::cout<<std::abs(z1max-z2min)<<" "<<std::abs(z1min-z2max)<<"\n";
  if(std::abs(z1max-z2min)>0.000000001 && std::abs(z1min-z2max)>0.000000001){
    /*Bug check, in the joined case one of these lists MUST have been reversed*/
    if(isRev1 == true && isRev2== true || isRev1 == false && isRev2== false){
      std::cout<<"Bug, for a pair of sections sharing a mutual point both lists have been reversed, this shjould not happen for continous curve and indicates a bug, please email christopher.prior@durham.ac.uk and send the data file you were using\n";
    }
    //std::cout<<"z values post chop joined "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n"; 
    /*need to mark if this si a local minimum of local maximum*/
    bool joinUp;
    bool algoBug=false;
    if(z1min ==z2min){
      joinUp = false;
    }else{
      if(z1max == z2max){
	joinUp = true;
      }else{
      algoBug= true;
      }
    }
    binaryFind BS;
    // now clip the lists to have the same z domain
   int integerwind=0;
   //std::cout<<"here joinUp "<<joinUp<<" "<<"algoBug "<<algoBug<<"\n";
   if(algoBug==false){
     if(joinUp == true){
       // local maximum, clipping on the minimum values
       if(z1min < z2min){
	 int low =0;
	 std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min);
	 point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
	 std::vector<point> sublistFin1(sublist1.begin()+ pair1.second,sublist1.begin()+s1Up+1);
	 sublistFin1.insert(sublistFin1.begin(),newLower);
	 sublist1 = sublistFin1;
       }else{
	 int low =0;
	 std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
	 point newLower = interp(sublist2[pair2.first],sublist2[pair2.second],z1min);
	 std::vector<point> sublistFin2(sublist2.begin() +pair2.second,sublist2.begin()+s2Up+1);
	 sublistFin2.insert(sublistFin2.begin(),newLower);
	 sublist2 = sublistFin2;
       }
     }else{
       if(z1max < z2max){
	 int low =0;
	 std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1max);
	 point newUpper = interp(sublist2[pair1.first],sublist2[pair1.second],z1max);
	 std::vector<point> sublistFin2(sublist2.begin(),sublist2.begin() +pair1.first+1);
	 sublistFin2.push_back(newUpper);
	 sublist2 = sublistFin2;
       }else{
	 int low =0;
	 std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max); 
	 point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
	 std::vector<point> sublistFin1(sublist1.begin(),sublist1.begin()+pair2.first+1);
	 sublistFin1.push_back(newUpper);
	 sublist1 = sublistFin1;
       }
     }
     /*std::cout<<sublist1[0].getX()<<" "<<sublist1[0].getY()<<" "<<sublist1[0].getZ()<<"\n";
	std::cout<<sublist2[0].getX()<<" "<<sublist2[0].getY()<<" "<<sublist2[0].getZ()<<"\n";
	std::cout<<"joined "<<sublist1[0].getZ()<<" "<<sublist1[sublist1.size()-1].getZ()<<" "<<sublist2[0].getZ()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";
	std::cout<<"is max ?"<<joinUp<<"\n";
	std::cout<<"size1 "<<sublist1.size()<<"\n";
	std::cout<<"size2 "<<sublist2.size()<<"\n";*/
     /*First check if the two sections now have a mutual z range, if not thorw a bug*/
     if(checkAlligned(sublist1,sublist2)==false){
       std::cout<<"A bug has occured and two sections of the curve involved in the non-local algorithm do not sharea  mutual range, pleas email christopher.prior@durham.ac.uk with the data file for which this occurred\n";
     }else{
       //begin searching for the mutual winding, two different rouitnes depending on a local min or local maximum
       int currQuad;
       int newQuad;
       point joinVec;
       point oldVec;
       double rotDirec;
       point zhat(0.0,0.0,1.0);
       std::pair<int,point> newQuadAndJoin;
       if(joinUp==true){
	 joinVec = sublist2[0].dif(sublist1[0]);
	 currQuad = quadrant(joinVec);
       }else{
	 joinVec = getTurningAngleQuad(s2minInp,tanList,curveList);
	 currQuad = quadrant(joinVec);
       }
       oldVec = joinVec;
       rotDirec = zhat.scalarTriple(zhat,joinVec,oldVec);
       std::pair<int,int> indicies(0,0);
       s1Up = sublist1.size()-1;
       s2Up = sublist2.size()-1;
       int maxIterations =s2Up +s1Up+10;
       int checker =1;
       if(joinUp ==false){
	 while(indicies.first<= s1Up && indicies.second<= s2Up && checker<maxIterations+10){
	   /*I believe that we should never have the case where on index is at its other limit and the other not, we check for this*/
	   if(indicies.first == s1Up && indicies.second < s2Up-1 ||indicies.first <s1Up-1 && indicies.second == s2Up){
	     std::cout<<"assumption error: the end of one list has been reached while the index of the other is not second to last. Please contact christopher.prior@durhma.ac.uk to report this bug\n";
	   }
	   if(indicies.first == s1Up && indicies.second == s2Up){
	     // the end of the lists, do nothing and push the indices to the end of the loop
	     indicies.first = indicies.first + 2*maxIterations;
	     indicies.second = indicies.second + 2*maxIterations;
	   }else{
	     newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	     branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	     rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	     currQuad =newQuadAndJoin.first;
	     //std::cout<<currQuad<<"\n";
	     oldVec = newQuadAndJoin.second;
	     checker++;
	   }
	 }
       }else{
	 /*If we are here we have a local maximum and we must iterate unitl both points are below the end points*/
	 while(indicies.first< s1Up && indicies.second< s2Up && checker<maxIterations+10){
	   /*I believer we should alwyas get to the point where both lists are one below thier ends for a local maximum, chekc this is true*/
	   if(indicies.first == s1Up || indicies.second == s2Up){
	     std::cout<<"assumption error, IN the non-local joined section routine, for a local maximum one list has reach the end whilst the other has not, this should not happen, please report this bug to christopher.prior@duham.ac.uk\n";
	   }
	   if(indicies.first==s1Up-1 &&indicies.second==s2Up-1){
	     joinVec = getTurningAngleQuad(s2minInp,tanList,curveList);
	     newQuad = quadrant(joinVec);
	     branchTracker(newQuad,currQuad,integerwind,rotDirec);
	     // and push the indices out of bounds
	     indicies.first = indicies.first + 2*maxIterations;
	     indicies.second = indicies.second + 2*maxIterations;
	   }else{
	     newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	     branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	     currQuad =newQuadAndJoin.first;
	     //std::cout<<currQuad<<"\n";
	     rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	     oldVec = newQuadAndJoin.second;
	     checker++;
	   }
      }
       }
     }
     
   }else{
     std::cout<<"bug: there is a bug in the non-local routine for joined sections. The mutual turning point has faieled to be recognised PLease email christopher.prior@durham.ac.uk\n";
     return 0.0;
   }
   //std::cout<<"joined algo"<<" "<<sigma1*sigma2*integerwind<<"\n";
   return sigma1*sigma2*integerwind;
  }else{
    //here we met at an inflection point
    //std::cout<<"hit the inflection tracker\n";
    return 0.0;
  }
  }else{
    return 0.0;
  }
}


int mutualWind2::getIntegerWindingOfPairJoinedClosedEnd(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp){
  point newS1minExtremum;
  point newJointExtremum;
  point newS2maxExtremum;
  int currQuad;
  std::pair<int,int> turnPair;
  /*find the interpolated extrema*/
  if(s1minInp.first != s1minInp.second){
    newS1minExtremum = getInterpolatedTop(s1minInp,tanList,curveList);
  }else{
    newS1minExtremum = curveList[s1minInp.first];
  }
  if(s2maxInp.first != s2maxInp.second){
    newS2maxExtremum = getInterpolatedTop(s2maxInp,tanList,curveList);
  }else{
    newS2maxExtremum = curveList[s2maxInp.first];
  }
  /*The turning point here this is an interpolation between the first and last points*/
  point tanz1 = tanList[0];
  point tanz2 = tanList[tanList.size()-1];
  double tanzv1 = tanz1.getZ();
  double tanzv2 = tanz2.getZ();
  double mids;
  if(std::abs(tanzv2-tanzv1)>0.01){
    mids = -tanzv1/(tanzv2-tanzv1);
  }else{
    mids=0.5;
  }
  /*double x =curveList[0].getX()+ tanz1.getX()*(mids-0.001) + (mids-0.001)*(mids-0.001)*(tanz2.getX()-tanz1.getX());
  double y = curveList[0].getY()+ tanz1.getY()*(mids-0.001) + (mids-0.001)*(mids-0.001)*(tanz2.getY()-tanz1.getY());
  double z = curveList[0].getZ()+ tanz1.getZ()*(mids-0.001) + (mids-0.001)*(mids-0.001)*(tanz2.getZ()-tanz2.getZ());
  double x2 =curveList[0].getX()+ tanz1.getX()*(mids+0.001) + (mids+0.001)*(mids+0.001)*(tanz2.getX()-tanz1.getX());
  double y2 = curveList[0].getY()+ tanz1.getY()*(mids+0.001) + (mids+0.001)*(mids+0.001)*(tanz2.getY()-tanz1.getY());
  double z2 = curveList[0].getZ()+ tanz1.getZ()*(mids+0.001) + (mids+0.001)*(mids+0.001)*(tanz2.getZ()-tanz2.getZ());*/
  double x =curveList[0].getX()+ tanz1.getX()*(mids) + (mids)*(mids)*(tanz2.getX()-tanz1.getX());
  double y = curveList[0].getY()+ tanz1.getY()*(mids) + (mids)*(mids)*(tanz2.getY()-tanz1.getY());
  double z = curveList[0].getZ()+ tanz1.getZ()*(mids) + (mids)*(mids)*(tanz2.getZ()-tanz2.getZ());
  point out(x,y,z);
  //point out2(x2,y2,z2);
  /*point toptan = tanz1.getX()*(mids) + *(tanz2.getX()-tanz1.getX())*tanz2;
    int endQuad = quadrant(toptan);*/
  newJointExtremum =  out;
  //get the relevant parts of the list
  std::vector<point> sublist1(curveList.begin()+s1minInp.second,curveList.begin()+s1maxInp.first +1);
  std::vector<point> sublist2(curveList.begin()+s2minInp.second,curveList.begin()+s2maxInp.first +1);
  AppendExtremumJoin(sublist1,sublist2,newS1minExtremum,newS2maxExtremum,newJointExtremum);
  //alter to account for our interpolated extrem
  /*std::cout<<"curve1"<<"\n";
  for(int i=0;i<4;i++){
    sublist1[i].printPoint();
  }
  std::cout<<"curve2"<<"\n";
  for(int i=sublist2.size()-3;i<sublist2.size();i++){
    sublist2[i].printPoint();
  }*/
  int s1Up = sublist1.size()-1;
  int s2Up = sublist2.size()-1;
  // we want both lists to be ordered so them move up in z, check if they need to an reverse if necessary
  bool isRev1 =false;
  int  sigma1=1;
  //check the sections have a sufficinet size
  if(std::abs(sublist1[s1Up].getZ()-sublist1[0].getZ())>0.0000001 && std::abs(sublist2[s2Up].getZ()-sublist2[0].getZ())>0.0000001){
  if(sublist1[s1Up].getZ()<sublist1[0].getZ()){
    std::reverse(sublist1.begin(),sublist1.end());
    isRev1 = true;
    sigma1 = -1;
  }
  //std::cout<<isRev1<<"\n";
  bool isRev2 =false;
  int  sigma2=1;
  if(sublist2[s2Up].getZ()<sublist2[0].getZ()){
    std::reverse(sublist2.begin(),sublist2.end());
    isRev2 = true;
    sigma2 = -1;
  }
  double z1min = sublist1[0].getZ();
  double z2min = sublist2[0].getZ();
  double z1max = sublist1[s1Up].getZ();
  double z2max = sublist2[s2Up].getZ();
  //check for an infecltion point
  if(std::abs(z1max-z2min)>0.000000001 && std::abs(z1min-z2max)>0.000000001){
    /*Bug check, in the joined case one of these lists MUST have been reversed*/
    if(isRev1 == true && isRev2== true || isRev1 == false && isRev2== false){
      std::cout<<"Bug, for a pair of sections sharing a mutual point both lists have been reversed, this shjould not happen for continous curve and indicates a bug, please email christopher.prior@durham.ac.uk and send the data file you were using\n";
    }
    //std::cout<<"z values post chop joined "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n"; 
    /*need to mark if this si a local minimum of local maximum*/
    bool joinUp=true;
    bool algoBug=false;
    binaryFind BS;
    // now clip the lists to have the same z domain
   int integerwind=0;
   //std::cout<<"here joinUp "<<joinUp<<"\n";
   if(algoBug==false){
       // local maximum, clipping on the minimum values
       if(z1min < z2min){
	 int low =0;
	 std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min);
	 point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
	 std::vector<point> sublistFin1(sublist1.begin()+ pair1.second,sublist1.begin()+s1Up+1);
	 sublistFin1.insert(sublistFin1.begin(),newLower);
	 sublist1 = sublistFin1;
       }else{
	 int low =0;
	 std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
	 point newLower = interp(sublist2[pair2.first],sublist2[pair2.second],z1min);
	 std::vector<point> sublistFin2(sublist2.begin() +pair2.second,sublist2.begin()+s2Up+1);
	 sublistFin2.insert(sublistFin2.begin(),newLower);
	 sublist2 = sublistFin2;
       }
       if(z1max < z2max){
	 int low =0;
	 std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1max);
	 point newUpper = interp(sublist2[pair1.first],sublist2[pair1.second],z1max);
	 std::vector<point> sublistFin2(sublist2.begin(),sublist2.begin() +pair1.first+1);
	 sublistFin2.push_back(newUpper);
	 sublist2 = sublistFin2;
       }else{
	 int low =0;
	 std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max); 
	 point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
	 std::vector<point> sublistFin1(sublist1.begin(),sublist1.begin()+pair2.first+1);
	 sublistFin1.push_back(newUpper);
	 sublist1 = sublistFin1;
       }
       /* std::cout<<sublist1[sublist1.size()-1].getX()<<" "<<sublist1[sublist1.size()-1].getY()<<" "<<sublist1[sublist1.size()-1].getZ()<<"\n";
	std::cout<<sublist2[sublist2.size()-1].getX()<<" "<<sublist2[sublist2.size()-1].getY()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";
	std::cout<<"joined "<<sublist1[0].getZ()<<" "<<sublist1[sublist1.size()-1].getZ()<<" "<<sublist2[0].getZ()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";
	point enddif = sublist1[sublist1.size()-1]-sublist2[sublist2.size()-1];
	std::cout<<"here the end angle is "<<angle(enddif.getX(),enddif.getY())<<"\n";
	std::cout<<"is max ?"<<joinUp<<"\n";
	std::cout<<"size1 "<<sublist1.size()<<"\n";
	std::cout<<"size2 "<<sublist2.size()<<"\n";*/
     /*First check if the two sections now have a mutual z range, if not thorw a bug*/
     if(false){
       std::cout<<"Closed Case: A bug has occured and two sections of the curve involved in the non-local algorithm do not sharea  mutual range, pleas email christopher.prior@durham.ac.uk with the data file for which this occurred\n";
     }else{
       //begin searching for the mutual winding, two different rouitnes depending on a local min or local maximum
       //int currQuad;
       int newQuad;
       point joinVec;
       point oldVec;
       double rotDirec;
       point zhat(0.0,0.0,1.0);
       std::pair<int,point> newQuadAndJoin;
       //std::cout<<"classed as local min or max "<<joinUp<<"\n";
       if(joinUp==true){
	 joinVec = sublist2[0].dif(sublist1[0]);
	 currQuad = quadrant(joinVec);
       }else{
	 joinVec = getTurningAngleQuad(s2minInp,tanList,curveList);
	 currQuad = quadrant(joinVec);
       }
       oldVec = joinVec;
       rotDirec = zhat.scalarTriple(zhat,joinVec,oldVec);
       std::pair<int,int> indicies(0,0);
       s1Up = sublist1.size()-1;
       s2Up = sublist2.size()-1;
       int maxIterations =s2Up +s1Up+10;
       int checker =1;
       if(joinUp ==false){
	 while(indicies.first<= s1Up && indicies.second<= s2Up && checker<maxIterations+10){
	   /*I believe that we should never have the case where on index is at its other limit and the other not, we check for this*/
	   if(indicies.first == s1Up && indicies.second < s2Up-1 ||indicies.first <s1Up-1 && indicies.second == s2Up){
	     std::cout<<"assumption error: the end of one list has been reached while the index of the other is not second to last. Please contact christopher.prior@durhma.ac.uk to report this bug\n";
	   }
	   if(indicies.first == s1Up && indicies.second == s2Up){
	     // the end of the lists, do nothing and push the indices to the end of the loop
	     indicies.first = indicies.first + 2*maxIterations;
	     indicies.second = indicies.second + 2*maxIterations;
	   }else{
	     newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	     branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	     rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	     currQuad =newQuadAndJoin.first;
	     //std::cout<<currQuad<<"\n";
	     oldVec = newQuadAndJoin.second;
	     checker++;
	   }
	 }
       }else{
	 /*If we are here we have a local maximum and we must iterate unitl both points are below the end points*/
	 while(indicies.first< s1Up && indicies.second< s2Up && checker<maxIterations+10){
	   /*I believer we should alwyas get to the point where both lists are one below thier ends for a local maximum, chekc this is true*/
	   if(indicies.first == s1Up || indicies.second == s2Up){
	     std::cout<<"assumption error, IN the non-local joined section routine, for a local maximum one list has reach the end whilst the other has not, this should not happen, please report this bug to christopher.prior@duham.ac.uk\n";
	   }
	   if(indicies.first==s1Up-1 &&indicies.second==s2Up-1){
	     joinVec = getTurningAngleQuad(s2minInp,tanList,curveList);
	     newQuad = quadrant(joinVec);
	     branchTracker(newQuad,currQuad,integerwind,rotDirec);
	     // and push the indices out of bounds
	     indicies.first = indicies.first + 2*maxIterations;
	     indicies.second = indicies.second + 2*maxIterations;
	   }else{
	     newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	     branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	     currQuad =newQuadAndJoin.first;
	     //std::cout<<currQuad<<"\n";
	     rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	     oldVec = newQuadAndJoin.second;
	     checker++;
	   }
	 }
       }
     }
   }else{
     std::cout<<"bug: there is a bug in the non-local routine for joined sections. The mutual turning point has faieled to be recognised PLease email christopher.prior@durham.ac.uk\n";
     return 0.0;
   }
   //for this closed curve end join we chekc for the last and final change
   /*std::cout<<"curr quad "<<currQuad<<"\n";
   std::cout<<"final quad "<<endQuad<<"\n";
   if(endQuad<3  && currQuad>2){
	integerwind =integerwind +1;
   }
   if(endQuad>2  && currQuad<3){
     integerwind =integerwind -1;
     }*/
   //std::cout<<"joined algo"<<" "<<sigma1*sigma2*integerwind<<"\n";
   return sigma1*sigma2*integerwind;
  }else{
    //here we met at an inflection point
    //std::cout<<"hit the inflection tracker\n";
    return 0.0;
  }
  }else{
    return 0.0;
  }
}

int mutualWind2::getIntegerWindingOfPair(std::vector<point> &curveList,std::vector<point> &tanList,std::pair<int,int> &s1minInp,std::pair<int,int> &s1maxInp,std::pair<int,int> &s2minInp,std::pair<int,int> &s2maxInp){
  //first reshape he lists to have there maxima and minima as the tangent interpolated ends
  point newS1minExtremum;
  point newS1maxExtremum;
  point newS2minExtremum;
  point newS2maxExtremum;
  std::pair<int,int> turnPair;
  //std::cout<<"mindexes "<<s1minInp.first<<" "<<s1minInp.second<<" "<<s1maxInp.first<<" "<<s1maxInp.second<<"\n";
  if(s1minInp.first != s1minInp.second){
    newS1minExtremum = getInterpolatedTop(s1minInp,tanList,curveList);
  }else{
    newS1minExtremum = curveList[s1minInp.first];
  }
  if(s1maxInp.first != s1maxInp.second){
    newS1maxExtremum = getInterpolatedTop(s1maxInp,tanList,curveList);
  }else{
    newS1maxExtremum = curveList[s1maxInp.first];
  }
  if(s2minInp.first != s2minInp.second){
    newS2minExtremum = getInterpolatedTop(s2minInp,tanList,curveList);
  }else{
    newS2minExtremum =  curveList[s2minInp.first];
  }
  if(s2maxInp.first != s2maxInp.second){
    newS2maxExtremum = getInterpolatedTop(s2maxInp,tanList,curveList);
  }else{
    newS2maxExtremum =  curveList[s2maxInp.first];
  }
  std::vector<point> sublist1(curveList.begin()+s1minInp.second,curveList.begin()+s1maxInp.first+1);
  std::vector<point> sublist2(curveList.begin()+s2minInp.second,curveList.begin()+s2maxInp.first+1);
  //std::cout<<"pre chop sublist non joined "<<sublist1[0].getZ()<<" "<<sublist1[sublist1.size()-1].getZ()<<" "<<sublist2[0].getZ()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";
  //alter to account for our interpolated extrema
  checkAppendExtrema(sublist1,newS1minExtremum,newS1maxExtremum);
  checkAppendExtrema(sublist2,newS2minExtremum,newS2maxExtremum);
  /*std::cout<<"extrema non-joined "<<newS1minExtremum.getZ()<<" "<<newS1maxExtremum.getZ()<<" "<<newS2minExtremum.getZ()<<" "<<newS2maxExtremum.getZ()<<"\n";
  std::cout<<"just after chop sublist non joined "<<sublist1[0].getZ()<<" "<<sublist1[sublist1.size()-1].getZ()<<" "<<sublist2[0].getZ()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";*/
  int s1Up = sublist1.size()-1;
  int s2Up = sublist2.size()-1;
  // reverse the lists if necessary
  bool isRev1 =false;
  int  sigma1=1;
  if(sublist1[s1Up].getZ()<sublist1[0].getZ()){
    std::reverse(sublist1.begin(),sublist1.end());
    isRev1 = true;
    sigma1 = -1;
  }
  bool isRev2 =false;
  int  sigma2=1;
  if(sublist2[s2Up].getZ()<sublist2[0].getZ()){
    std::reverse(sublist2.begin(),sublist2.end());
    isRev2 = true;
    sigma2 = -1;
  }
  double z1min = sublist1[0].getZ();
  double z2min = sublist2[0].getZ();
  double z1max = sublist1[s1Up].getZ();
  double z2max = sublist2[s2Up].getZ();
  //std::cout<<"pre-chop "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n";
  /*for(int i=0;i<sublist1.size();i++){
      std::cout<<"pre weird insertion "<<sublist1[i].getZ()<<"\n";
    }
    for(int i=0;i<sublist2.size();i++){
      std::cout<<"pre weird insertion2 "<<sublist2[i].getZ()<<"\n";
    }
  */
  int integerWind;
  bool getWind;
  binaryFind BS;
  /*Now chop the domains so that they mutually ovelrap*/
  if(z1max<= z2min || z2max<= z1min){
    // no overlap of the domains
    integerWind =0;
    getWind=false;
  }else if(z1min<=z2min && z1max>= z2max){
    // [z2min,z2max] \subset [z1min,z1max]
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
    point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
    point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
    std::vector<point> sublistFin1(sublist1.begin()+pair1.second,sublist1.begin()+ pair2.first+1);
    sublistFin1.insert(sublistFin1.begin(),newLower);
    sublistFin1.push_back(newUpper);
    sublist1 = sublistFin1;
    getWind=true;
  }else if(z2min<=z1min && z2max>= z1max){
    // [z1min,z1max] \subset [z2min,z2max]
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
    point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
    point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
    std::vector<point> sublistFin2(sublist2.begin()+ pair1.second,sublist2.begin() +pair2.first+1);
    sublistFin2.insert(sublistFin2.begin(),newLower);
    sublistFin2.push_back(newUpper);
    sublist2 = sublistFin2;
    getWind=true;
  }
  else if(z2min<=z1max && z2min>= z1min){
    // neither section contians the other but z1max is below z2max
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
    point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
    point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
    std::vector<point> sublistFin1(sublist1.begin()+ pair1.second, sublist1.begin() + s1Up+1);
    std::vector<point> sublistFin2(sublist2.begin(),sublist2.begin()+pair2.first+1);
    sublistFin1.insert(sublistFin1.begin(),newLower);
    sublistFin2.push_back(newUpper);
    sublist2 = sublistFin2;
    sublist1= sublistFin1;
    getWind=true;
  }
  else if(z1min<=z2max && z1min>= z2min){
    // neither section contians the other but z2max is below z1max
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
    point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
    point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
    std::vector<point> sublistFin2(sublist2.begin()+pair1.second,sublist2.begin()+s2Up+1);
    std::vector<point> sublistFin1(sublist1.begin(),sublist1.begin()+pair2.first+1);
    sublistFin2.insert(sublistFin2.begin(),newLower);
    sublistFin1.push_back(newUpper);
    sublist2 = sublistFin2;
    sublist1= sublistFin1;
    getWind=true;
  }
  //std::cout<<"non joined "<<sublist1[0].getZ()<<" "<<sublist1[sublist1.size()-1].getZ()<<" "<<sublist2[0].getZ()<<" "<<sublist2[sublist2.size()-1].getZ()<<"\n";
  int integerwind =0;
  if(getWind){
    /*for(int i=0;i<sublist1.size();i++){
      std::cout<<"weird insertion "<<sublist1[i].getZ()<<"\n";
    }
    for(int i=0;i<sublist2.size();i++){
      std::cout<<"weird insertion2 "<<sublist2[i].getZ()<<"\n";
      }*/
    // here we have overlap and we search for an integer winding
     if(checkAlligned(sublist1,sublist2)==false){
       std::cout<<"A bug has occured and two sections of the curve involved in the non-local algorithm do not sharea  mutual range, pleas email christopher.prior@durham.ac.uk with the data file for which this occurred\n";
      }else{
        int currQuad;
	int newQuad;
	point joinVec;
	point oldVec;
	double rotDirec;
	point zhat(0.0,0.0,1.0);
	std::pair<int,point> newQuadAndJoin;
	joinVec = sublist2[0].dif(sublist1[0]);
	currQuad = quadrant(joinVec);
	oldVec = joinVec;
	rotDirec = 0.0;
	std::pair<int,int> indicies(0,0);
	s1Up = sublist1.size()-1;
	s2Up = sublist2.size()-1;
	int maxIterations =s2Up +s1Up+10;
	int checker =1;
	while(indicies.first<= s1Up && indicies.second<= s2Up && checker<maxIterations+10){
	/*I believe that we should never have the case where on index is at its other limit and the other not, we check for this*/
	  //std::cout<<"indcheck "<<indicies.first<<" "<<indicies.second<<" "<<s1Up<<" "<<s2Up<<"\n";
	  if(indicies.first == s1Up && indicies.second < s2Up-1 ||indicies.first <s1Up-1 && indicies.second == s2Up){
	  std::cout<<"assumption error: the end of one list has been reached while the index of the other is not second to last. Please contact christopher.prior@durhma.ac.uk to report this bug\n";
	  }
	  //std::cout<<sublist1[indicies.first].getZ()<<" "<<sublist2[indicies.second].getZ()<<"\n"; 
	  if(indicies.first == s1Up && indicies.second == s2Up){
	  // the end of the lists, do nothing and push the indices to the end of the loop
	    indicies.first = indicies.first + 2*maxIterations;
	    indicies.second = indicies.second + 2*maxIterations;
	  }else{
	    newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	    branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	    rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	    currQuad =newQuadAndJoin.first;
	    oldVec = newQuadAndJoin.second;
	    checker++;
	  }
	}
     }
  }
  //std::cout<<"nonjoined algo"<<" "<<sigma1*sigma2*integerwind<<"\n";
  return sigma1*sigma2*integerwind;
}


int mutualWind2::getIntegerWindingOfPairLink(std::vector<point> &curveList1,std::vector<point> &curveList2,int &s1minInp,int &s1maxInp,int &s2minInp,int &s2maxInp){
  //for the linking number no interpolation is needed
  point newS1minExtremum;
  point newS1maxExtremum;
  point newS2minExtremum;
  point newS2maxExtremum;
  std::vector<point> sublist1(curveList1.begin()+s1minInp,curveList1.begin()+s1maxInp+1);
  std::vector<point> sublist2(curveList2.begin()+s2minInp,curveList2.begin()+s2maxInp+1);
  int s1Up = sublist1.size()-1;
  int s2Up = sublist2.size()-1;
  // reverse the lists if necessary
  bool isRev1 =false;
  int  sigma1=1;
  if(sublist1[s1Up].getZ()<sublist1[0].getZ()){
    std::reverse(sublist1.begin(),sublist1.end());
    isRev1 = true;
    sigma1 = -1;
  }
  bool isRev2 =false;
  int  sigma2=1;
  if(sublist2[s2Up].getZ()<sublist2[0].getZ()){
    std::reverse(sublist2.begin(),sublist2.end());
    isRev2 = true;
    sigma2 = -1;
  }
  double z1min = sublist1[0].getZ();
  double z2min = sublist2[0].getZ();
  double z1max = sublist1[s1Up].getZ();
  double z2max = sublist2[s2Up].getZ();
  //std::cout<<"z values post chop non joined "<<z1min<<" "<<z1max<<" "<<z2min<<" "<<z2max<<"\n"; 
  int integerWind;
  bool getWind;
  binaryFind BS;
  /*Now chop the domains so that they mutually ovelrap*/
  if(z1max<= z2min || z2max<= z1min){
    // no overlap of the domains
    integerWind =0;
    getWind=false;
  }else if(z1min<=z2min && z1max>= z2max){
    // [z2min,z2max] \subset [z1min,z1max]
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
    point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
    point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
    std::vector<point> sublistFin1(sublist1.begin()+pair1.second,sublist1.begin()+ pair2.first+1);
    sublistFin1.insert(sublistFin1.begin(),newLower);
    sublistFin1.push_back(newUpper);
    sublist1 = sublistFin1;
    getWind=true;
  }else if(z2min<=z1min && z2max>= z1max){
    // [z1min,z1max] \subset [z2min,z2max]
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
    point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
    point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
    std::vector<point> sublistFin2(sublist2.begin()+ pair1.second,sublist2.begin() +pair2.first+1);
    sublistFin2.insert(sublistFin2.begin(),newLower);
    sublistFin2.push_back(newUpper);
    sublist2 = sublistFin2;
    getWind=true;
  }
  else if(z2min<=z1max && z2min>= z1min){
    // neither section contians the other but z2max is below z1max
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist1,low,s1Up,z2min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist2,low,s2Up,z1max);
    point newLower = interp(sublist1[pair1.first],sublist1[pair1.second],z2min);
    point newUpper = interp(sublist2[pair2.first],sublist2[pair2.second],z1max);
    std::vector<point> sublistFin1(sublist1.begin()+ pair1.second, sublist1.begin() + s1Up+1);
    std::vector<point> sublistFin2(sublist2.begin(),sublist2.begin()+pair2.first+1);
    sublistFin1.insert(sublistFin1.begin(),newLower);
    sublistFin2.push_back(newUpper);
    sublist2 = sublistFin2;
    sublist1= sublistFin1;
    getWind=true;
  }
  else if(z1min<=z2max && z1min>= z2min){
    // neither section contians the other but z1max is below z2max
    int low = 0;
    std::pair<int,int> pair1  = BS.getContainingPair(sublist2,low,s2Up,z1min); 
    std::pair<int,int> pair2  = BS.getContainingPair(sublist1,low,s1Up,z2max);
    point newLower = interp(sublist2[pair1.first],sublist2[pair1.second],z1min);
    point newUpper = interp(sublist1[pair2.first],sublist1[pair2.second],z2max);
    std::vector<point> sublistFin2(sublist2.begin()+pair1.second,sublist2.begin()+s2Up+1);
    std::vector<point> sublistFin1(sublist1.begin(),sublist1.begin()+pair2.first+1);
    sublistFin2.insert(sublistFin2.begin(),newLower);
    sublistFin1.push_back(newUpper);
    sublist2 = sublistFin2;
    sublist1= sublistFin1;
    getWind=true;
  }
  
  int integerwind =0;
  if(getWind){
    // here we have overlap and we search for an integer winding
     if(checkAlligned(sublist1,sublist2)==false){
       std::cout<<"A bug has occured and two sections of the curve involved in the non-local algorithm do not sharea  mutual range, pleas email christopher.prior@durham.ac.uk with the data file for which this occurred\n";
      }else{
       int currQuad;
	int newQuad;
	point joinVec;
	point oldVec;
	double rotDirec;
	point zhat(0.0,0.0,1.0);
	std::pair<int,point> newQuadAndJoin;
	joinVec = sublist2[0].dif(sublist1[0]);
	currQuad = quadrant(joinVec);
	oldVec = joinVec;
	rotDirec = zhat.scalarTriple(zhat,joinVec,oldVec);
	std::pair<int,int> indicies(0,0);
	s1Up = sublist1.size()-1;
	s2Up = sublist2.size()-1;
	int maxIterations =s2Up +s1Up+10;
	int checker =1;
	while(indicies.first<= s1Up && indicies.second<= s2Up && checker<maxIterations+10){
	/*I believe that we should never have the case where on index is at its other limit and the other not, we check for this*/
	  if(indicies.first == s1Up && indicies.second < s2Up-1 ||indicies.first <s1Up-1 && indicies.second == s2Up){
	  std::cout<<"assumption error: the end of one list has been reached while the index of the other is not second to last. Please contact christopher.prior@durhma.ac.uk to report this bug\n";
	  }
	  if(indicies.first == s1Up && indicies.second == s2Up){
	  // the end of the lists, do nothing and push the indices to the end of the loop
	    indicies.first = indicies.first + 2*maxIterations;
	    indicies.second = indicies.second + 2*maxIterations;
	  }else{
	    newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	    branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	    rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	    currQuad =newQuadAndJoin.first;
	    oldVec = newQuadAndJoin.second;
	    checker++;
	  }
	}
     }
  }
  //std::cout<<"nonjoined algo"<<" "<<sigma1*sigma2*integerwind<<"\n";
  return sigma1*sigma2*integerwind;
}


int mutualWind2::getIntegerWindingOfPairLinkTesting(windingStorage &w){
  std::vector<point> sublist1 = w.getSublist1();
  std::vector<point> sublist2 = w.getSublist2();
  bool getWind = w.getIsWind();
  int sigma1 = w.getSigma1();
  int sigma2 = w.getSigma2();
  int integerwind =0;
  //sublist1.front().printPoint();
  //sublist1.back().printPoint();
  // sublist2.front().printPoint();
  //sublist2.back().printPoint();
  //std::cout<<" \n";
  if(getWind){
    // here we have overlap and we search for an integer winding
     if(checkAlligned(sublist1,sublist2)==false){
       std::cout<<"A bug has occured and two sections of the curve involved in the non-local algorithm do not share a  utual range, please email christopher.prior@durham.ac.uk with the data file for which this occurred\n";
      }else{
        int currQuad;
	int newQuad;
	point joinVec;
	point oldVec;
	double rotDirec;
	point zhat(0.0,0.0,1.0);
	std::pair<int,point> newQuadAndJoin;
	joinVec = sublist2[0].dif(sublist1[0]);
	currQuad = quadrant(joinVec);
	oldVec = joinVec;
	rotDirec = zhat.scalarTriple(zhat,joinVec,oldVec);
	std::pair<int,int> indicies(0,0);
	int s1Up = sublist1.size()-1;
	int s2Up = sublist2.size()-1;
	int maxIterations =s2Up +s1Up+10;
	int checker =1;
	while(indicies.first<= s1Up && indicies.second<= s2Up && checker<maxIterations+10){
	/*I believe that we should never have the case where on index is at its other limit and the other not, we check for this*/
	  if(indicies.first == s1Up && indicies.second < s2Up-1 ||indicies.first <s1Up-1 && indicies.second == s2Up){
	  std::cout<<"assumption error: the end of one list has been reached while the index of the other is not second to last. Please contact christopher.prior@durhma.ac.uk to report this bug\n";
	  }
	  if(indicies.first == s1Up && indicies.second == s2Up){
	    //std::cout<<"final quad is "<<currQuad<<"\n";
	  // the end of the lists, do nothing and push the indices to the end of the loop
	    indicies.first = indicies.first + 
2*maxIterations;
	    indicies.second = indicies.second + 2*maxIterations;
	  }else{
	    if(checker>1){
	      newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	      branchTracker(newQuadAndJoin.first,currQuad,integerwind,rotDirec);
	      rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	      currQuad =newQuadAndJoin.first;
	      /*std::cout<<"currQuad is "<<currQuad<<"\n";
		std::cout<<"integerWind is "<<integerwind<<"\n";*/
	      oldVec = newQuadAndJoin.second;
	      checker++;
	      //std::cout<<"checker is "<<checker<<" "<<indicies.first<<" "<<indicies.second<<" "<<" "<<s1Up<<" "<<s2Up<<"\n";
	    }
	    else{
	      // if the first point ignore the branch check as we cannot track direction.
	      newQuadAndJoin= getNextIndex(sublist1,sublist2,indicies,currQuad);
	      rotDirec = zhat.scalarTriple(zhat,newQuadAndJoin.second,oldVec);
	      oldVec = newQuadAndJoin.second;
	      checker++;
	    }
	  }
	}
     }
  }
  //std::cout<<"nonjoined algo"<<" "<<sigma1*sigma2*integerwind<<"\n";
  return sigma1*sigma2*integerwind;
}

void mutualWind2::extendAboveSpan(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2){
   std::pair<double,double> extremalPair;
   bool addBeggining2=false;
   bool addEnd2=false;
   bool addBeggining1=false;
   bool addEnd1=false;
   if(turnList2.size()>0){
     extremalPair = getExtremaZ(turnList2[0],tanList2,curveList2);
     int CS2 = curveList2.size()-1;
     double globalMax =extremalPair.first;
     double globalMin =extremalPair.second;
     for(int i=1;i<turnList2.size();i++){
       extremalPair = getExtremaZ(turnList2[i],tanList2,curveList2);
       if(extremalPair.first>globalMax){
	 globalMax =extremalPair.first;
       }
       if(extremalPair.second<globalMin){
	 globalMin =extremalPair.second;
       }
     } 
     /*now check if either of these is bigger/smaller then the curve's endpoints */
     if(curveList2[0].getZ() >globalMax){
       globalMax = curveList2[0].getZ();
     }
     if(curveList2[CS2].getZ()> globalMax){
       globalMax = curveList2[CS2].getZ();
     }
     if(curveList2[0].getZ() <globalMin){
       globalMin = curveList2[0].getZ();
     }
     if(curveList2[CS2].getZ()< globalMin){
       globalMin = curveList2[CS2].getZ();
     }
   /*Now check if we should extend curve1 */
  std::pair<int,int> pst(0,0);
  std::pair<int,int> pe(curveList1.size()-1,curveList1.size()-1);
  if((curveList1[1].getZ()-curveList1[0].getZ())<0.0){
    // here we have added an effective turning point the the curve's begnining (as the close moves n the same way as its tangent
    turnList1.insert(turnList1.begin(),pst);
    addBeggining1=true;
    };
  if((curveList1[curveList1.size()-1].getZ()-curveList1[curveList1.size()-2].getZ())<0.0){
    // here we have added an effective turning point the the curve's begnining (as the close moves n the smae way as its tangent
    turnList1.push_back(pe);
    addEnd1=true;
  };
  int CS = curveList1.size()-1;
   if(curveList1[0].getZ()>globalMin){
     addBeggining2=true;
      double minDif = globalMin - curveList1[0].getZ();
      int newPoints=100;
      double step = 1.5*(minDif-0.001)/double(newPoints);
      double origZ = curveList1[0].getZ();
      for(int i=1;i<= newPoints;i++){
	point p(curveList1[0].getX(),curveList1[0].getY(),origZ+ i*step);
	curveList1.insert(curveList1.begin(),p);
	point pt(0.0,0.0,1.0);
	tanList1.insert(tanList1.begin(),pt);
      }
   }
    // now check if the global maximum is gretaer than the last point, if so extend the curve down
   CS = curveList1.size()-1;
    if(curveList1[CS].getZ()<globalMax){
      double maxDif = globalMax - curveList1[CS].getZ();
      int newPoints=100;
      double step = 1.5*(maxDif+0.001)/double(newPoints);
      for(int i=1;i<= newPoints;i++){
	point p(curveList1[CS].getX(),curveList1[CS].getY(),curveList1[CS].getZ()+i*step);
	curveList1.push_back(p);
	point pt(0.0,0.0,1.0);
	tanList1.push_back(pt);
      }
      addEnd2=true;
    }
    if(addBeggining2){
       for(int i=0;i<turnList1.size();i++){
	 turnList1[i].first = turnList1[i].first + 100;
	 turnList1[i].second = turnList1[i].second + 100;
       }
    }
    /*Finally depending on the initial tangent it may be we have added a turning point.*/
  
   }
    //if we have added points to the beginnig of the list we must push back the turning point
   if(addBeggining2==true||addBeggining1==false){
      std::pair<int,int> pstart(0,0);
      turnList1.insert(turnList1.begin(),pstart);
    }
    if(addEnd2==true||addEnd1==false){
      std::pair<int,int> pend(curveList1.size()-1,curveList1.size()-1);
      turnList1.push_back(pend);
    }
}

void mutualWind2::extendBelowSpan(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2){
  bool addBeggining2=false;
   bool addEnd2=false;
   bool addBeggining1=false;
   bool addEnd1=false;
   if(turnList2.size()>0){
   std::pair<double,double> extremalPair;
   extremalPair = getExtremaZ(turnList2[0],tanList2,curveList2);
   double globalMax =extremalPair.first;
   double globalMin =extremalPair.second;
   for(int i=1;i<turnList2.size();i++){
      extremalPair = getExtremaZ(turnList2[i],tanList2,curveList2);
     if(extremalPair.first>globalMax){
	  globalMax =extremalPair.first;
       }
     if(extremalPair.second<globalMin){
	  globalMin =extremalPair.second;
       }
    } 
   /*now check if either of these is bigger/smaller then the curve's endpoints */
   int CS2 = curveList2.size()-1;
   if(curveList2[0].getZ() >globalMax){
     globalMax = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()> globalMax){
     globalMax = curveList2[CS2].getZ();
   }
   if(curveList2[0].getZ() <globalMin){
     globalMin = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()< globalMin){
     globalMin = curveList2[CS2].getZ();
   }
   std::pair<int,int> pst(0,0);
   std::pair<int,int> pe(curveList1.size()-1,curveList1.size()-1);
   if((curveList1[1].getZ()-curveList1[0].getZ())>0.0){
     // here we have added an effective turning point the the curve's begnining (as the close moves n the smae way as its tangent
     turnList1.insert(turnList1.begin(),pst);
     addBeggining1=true;
   };
   if((curveList1[curveList1.size()-1].getZ()-curveList1[curveList1.size()-2].getZ())>0.0){
     // here we have added an effective turning point the the curve's begnining (as the close moves n the smae way as its tangent
     turnList1.push_back(pe);
     addEnd1=true;
   };
   /*Now check if we should extend curve1 */
  int CS = curveList1.size()-1;
   if(curveList1[0].getZ() < globalMax){
      double maxDif = globalMax - curveList1[0].getZ();
      int newPoints=100;
      double step = 1.5*(maxDif+0.001)/double(newPoints);
      double origZ = curveList1[0].getZ();
      for(int i=1;i<=newPoints;i++){
	point p(curveList1[0].getX(),curveList1[0].getY(),origZ+ i*step);
	curveList1.insert(curveList1.begin(),p);
	point pt(0.0,0.0,1.0);
	tanList1.insert(tanList1.begin(),pt);
	addBeggining2=true;
      }
   }
   CS = curveList1.size()-1;
   // now check if the global maximum is gretaer than the last point, if so extend the curve down
   if(curveList1[CS].getZ() > globalMin){
     double minDif = globalMin - curveList1[CS].getZ();
      int newPoints=100;
      double step = 1.5*(minDif-0.001)/double(newPoints);
      for(int i=1;i<= newPoints;i++){
	point p(curveList1[CS].getX(),curveList1[CS].getY(),curveList1[CS].getZ()+i*step);
	curveList1.push_back(p);
	point pt(0.0,0.0,1.0);
	tanList1.push_back(pt);
      }
      addEnd2=true;
   }
   if(addBeggining2){
     for(int i=0;i<turnList1.size();i++){
       turnList1[i].first = turnList1[i].first + 100;
       turnList1[i].second = turnList1[i].second + 100;
     }
   }
   //if we have added points to the beginnig of the list
   if(addBeggining2==true||addBeggining1==false){
     std::pair<int,int> pstart(0,0);
     turnList1.insert(turnList1.begin(),pstart);
   }
   if(addEnd2==true||addEnd1==false){
     std::pair<int,int> pend(curveList1.size()-1,curveList1.size()-1);
      turnList1.push_back(pend);
   }
   }
}



void mutualWind2::extendBothPosNeg(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2){
   std::pair<double,double> extremalPair;
   extremalPair = getExtremaZ(turnList2[0],tanList2,curveList2);
   double globalMax =extremalPair.first;
   double globalMin =extremalPair.second;
   for(int i=1;i<turnList2.size();i++){
      extremalPair = getExtremaZ(turnList2[i],tanList2,curveList2);
     if(extremalPair.first>globalMax){
	  globalMax =extremalPair.first;
       }
     if(extremalPair.second<globalMin){
	  globalMin =extremalPair.second;
       }
    } 
   /*now check if either of these is bigger/smaller then the curve's endpoints */
   int CS2 = curveList2.size()-1;
   if(curveList2[0].getZ() >globalMax){
     globalMax = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()> globalMax){
     globalMax = curveList2[CS2].getZ();
   }
   if(curveList2[0].getZ() <globalMin){
     globalMin = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()< globalMin){
     globalMin = curveList2[CS2].getZ();
   }
   /*Now check if we should extend curve1 */
   std::vector<point> subcurve;
   bool addBeggining=false;
   if(curveList1[0].getZ() < globalMax){
      double maxDif = globalMax - curveList1[0].getZ();
      int newPoints=100;
      double step = 1.5*(maxDif+0.001)/100.0;
      double origZ = curveList1[0].getZ();
      for(int i=1;i<= 100;i++){
	point p(curveList1[0].getX(),curveList1[0].getY(),origZ+ i*step);
	subcurve.push_back(p);
	point pt(0.0,0.0,1.0);
	tanList1.insert(tanList1.begin(),pt);
	addBeggining=true;
      }
      std::reverse(subcurve.begin(),subcurve.end()); 
      curveList1.insert(curveList1.begin(),subcurve.begin(),subcurve.end());
   }
    int CS = curveList1.size()-1;
    // now check if the global maximum is gretaer than the last point, if so extend the curve down
    if(curveList1[CS].getZ() > globalMin){
      double minDif = globalMin - curveList1[CS].getZ();
      int newPoints=100;
      double step = 1.5*(minDif-0.001)/100.0;
      for(int i=1;i<= 100;i++){
	point p(curveList1[CS].getX(),curveList1[CS].getY(),curveList1[CS].getZ()+i*step);
	curveList1.push_back(p);
	point pt(0.0,0.0,-1.0);
	tanList1.push_back(pt);
      }
    }
    if(addBeggining){
       for(int i=0;i<turnList1.size();i++){
	 turnList1[i].first = turnList1[i].first + 100;
	 turnList1[i].second = turnList1[i].second + 100;
       }
     }
     //if we have added points to the beginnig of the list we must push back the turning points
    std::pair<int,int> pstart(0,0);
    std::pair<int,int> pend(curveList1.size()-1,curveList1.size()-1);
    turnList1.insert(turnList1.begin(),pstart);
    turnList1.push_back(pend);
}

void mutualWind2::extendBothNegPos(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<point> &tanList1,std::vector<point> &tanList2,std::vector<std::pair<int,int> > &turnList1,std::vector<std::pair<int,int> > &turnList2){
   std::pair<double,double> extremalPair;
   extremalPair = getExtremaZ(turnList2[0],tanList2,curveList2);
   double globalMax =extremalPair.first;
   double globalMin =extremalPair.second;
   for(int i=1;i<turnList2.size();i++){
      extremalPair = getExtremaZ(turnList2[i],tanList2,curveList2);
     if(extremalPair.first>globalMax){
	  globalMax =extremalPair.first;
       }
     if(extremalPair.second<globalMin){
	  globalMin =extremalPair.second;
       }
    } 
   /*now check if either of these is bigger/smaller then the curve's endpoints */
   int CS2 = curveList2.size()-1;
   if(curveList2[0].getZ() >globalMax){
     globalMax = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()> globalMax){
     globalMax = curveList2[CS2].getZ();
   }
   if(curveList2[0].getZ() <globalMin){
     globalMin = curveList2[0].getZ();
   }
   if(curveList2[CS2].getZ()< globalMin){
     globalMin = curveList2[CS2].getZ();
   }
   /*Now check if we should extend curve1 */
  bool addBeggining=false;
  std::vector<point> subcurve;
   if(curveList1[0].getZ() < globalMin){
      double maxDif = globalMin - curveList1[0].getZ();
      int newPoints=100;
      double step = 1.5*(0.001-maxDif)/100.0;
      double origZ = curveList1[0].getZ();
      for(int i=1;i<= 100;i++){
	point p(curveList1[0].getX(),curveList1[0].getY(),origZ+ i*step);
	subcurve.push_back(p);
	point pt(0.0,0.0,-1.0);
	tanList1.insert(tanList1.begin(),pt);
	addBeggining=true;
      }
      std::reverse(subcurve.begin(),subcurve.end()); 
      curveList1.insert(curveList1.begin(),subcurve.begin(),subcurve.end());
      
   }
    int CS = curveList1.size()-1;
    // now check if the global maximum is gretaer than the last point, if so extend the curve down
    if(curveList1[CS].getZ() < globalMax){
      double maxDif = globalMax - curveList1[CS].getZ();
      int newPoints=100;
      double step = 1.5*(maxDif+0.001)/100.0;
      for(int i=1;i<= 100;i++){
	point p(curveList1[CS].getX(),curveList1[CS].getY(),curveList1[CS].getZ()+i*step);
	curveList1.push_back(p);
	point pt(0.0,0.0,1.0);
	tanList1.push_back(pt);
      }
    }
    if(addBeggining){
       for(int i=0;i<turnList1.size();i++){
	 turnList1[i].first = turnList1[i].first + 100;
	 turnList1[i].second = turnList1[i].second + 100;
       }
     }
     //if we have added points to the beginnig of the list we must push back the turning points
    std::pair<int,int> pstart(0,0);
    std::pair<int,int> pend(curveList1.size()-1,curveList1.size()-1);
    turnList1.insert(turnList1.begin(),pstart);
    turnList1.push_back(pend);
}

std::vector<int> mutualWind2::getIntegerWindingsStarGen(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList){
    // find the maximum and minimum values amongst the turning points, we need to chekc if they are above he ends
    int curveSize = curveList.size()-1;
    int tanSize=tanList.size()-1;
    /*if(tanList[0].getZ()>= 0 && tanList[tanSize].getZ()>= 0){
      extendAboveSpan(curveList,curveList,tanList,tanList,turnList,turnList);
    }else if(tanList[0].getZ()< 0 && tanList[tanSize].getZ()<0){
      //std::cout<<"here2 \n";
      extendBelowSpan(curveList,curveList,tanList,tanList,turnList,turnList);
    }else if(tanList[0].getZ()>= 0 && tanList[tanSize].getZ()<0){
      //std::cout<<"here3 \n";
      extendBothPosNeg(curveList,curveList,tanList,tanList,turnList,turnList);
    }else if(tanList[0].getZ()< 0 && tanList[tanSize].getZ()>= 0){
      //std::cout<<"here4 \n";
      extendBothNegPos(curveList,curveList,tanList,tanList,turnList,turnList);
    }*/
    //pre compute the without extesion end non-local contirbutions
    if(curveList[0].getZ() < curveList[curveSize].getZ()){
      extendAboveSpan(curveList,curveList,tanList,tanList,turnList,turnList);
    }else{
      extendBelowSpan(curveList,curveList,tanList,tanList,turnList,turnList);
    }
    /*for(int i=turnList[0].first;i<=turnList[1].first;i++){
      curveList[i].printPoint();
     }
    for(int i=turnList[1].first;i<=turnList[2].first;i++){
      curveList[i].printPoint();
      }*/
    /*We are now ready to calculate the integer widnings*/
    //std::cout<<starContributionsIndex.size()<<" "<<turnList.size()-oldTurnList.size()<<"\n";
    std::vector<int> output;
    if(turnList.size() >2){
      for(int k=0;k<=turnList.size()-3;k++){
	      for(int l=k+1;l<=turnList.size()-2;l++){
	  //std::cout<<"indicies "<<k<<" "<<k+1<<" "<<l<<" "<<l+1<<"\n";
	  //std::cout<<turnList[k].first<<" "<<turnList[k+1].first<<" "<<turnList[l].first<<" "<<turnList[l+1].first<<"\n";
	    if(l==k+1){
	      int currInt = getIntegerWindingOfPairJoined(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	    output.push_back(currInt);
	  }else{
	    if((turnList[k+1].first-turnList[k].first)>1 && turnList[l+1].first-turnList[l].first>1){
	      int currInt = getIntegerWindingOfPair(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	      output.push_back(currInt);
	    }else{
	      // section too small
	      int currInt = 0;
	      output.push_back(currInt);
	    }
	  }
	  //std::cout<<integerWind<<"\n";
	}
      }
    }    //integerWind = integerWind + getIntegerWindingOfPairJoined(curveList,tanList,turnList[3],turnList[4],turnList[4],turnList[5]);
    return output;
  }

std::vector<int> mutualWind2::getIntegerWindingsClosed(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList){
    // find the maximum and minimum values amongst the turning points, we need to chekc if they are above he ends
    std::vector<int> output;
    std::pair<int,int> pstart(0,0);
    std::pair<int,int> pend(curveList.size()-1,curveList.size()-1);
    turnList.insert(turnList.begin(),pstart);
    turnList.push_back(pend);
    //std::cout<<"curve lentgth "<<curveList.size()<<"\n";
    //std::cout<<"tan length "<<tanList.size()<<"\n";
    //std::cout<<turnList.size()<<"\n";
    if(turnList.size() >2){
      for(int k=0;k<=turnList.size()-3;k++){
	    for(int l=k+1;l<=turnList.size()-2;l++){
	  //if(l==k+1 || (k==0 && l==(turnList.size()-2))){
    if(l==k+1){
	    int currInt = getIntegerWindingOfPairJoined(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	    output.push_back(currInt);
	  }else{
     if((k==0 && l==(turnList.size()-2))){
       int currInt = getIntegerWindingOfPairJoinedClosedEnd(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
       //std::cout<<"closure integer "<<currInt<<"\n";
	    output.push_back(currInt);
     }else{
	   if((turnList[k+1].first-turnList[k].first)>1 && turnList[l+1].first-turnList[l].first>1){
	      int currInt = getIntegerWindingOfPair(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	      output.push_back(currInt);
	    }else{
	      // section too small
	      int currInt = 0;
	      output.push_back(currInt);
	    }
     } 
	  }
	  //std::cout<<integerWind<<"\n";
	}
      }
    }
    return output;
}

std::vector<int> mutualWind2::getIntegerWindingsStarGenStandard(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList,std::vector<double> &startAngles,std::vector<double> &endAngles){
    // find the maximum and minimum values amongst the turning points, we need to chekc if they are above he ends
    std::vector<int> output;
     std::pair<int,int> pstart(0,0);
    std::pair<int,int> pend(curveList.size()-1,curveList.size()-1);
    turnList.insert(turnList.begin(),pstart);
    turnList.push_back(pend);
    if(turnList.size() >2){
      for(int k=0;k<=turnList.size()-3;k++){
	for(int l=k+1;l<=turnList.size()-2;l++){
	  if(l==k+1){
	    int currInt = getIntegerWindingOfPairJoined(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	    output.push_back(currInt);
	  }else{
	   if((turnList[k+1].first-turnList[k].first)>1 && turnList[l+1].first-turnList[l].first>1){
	      int currInt = getIntegerWindingOfPair(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	      output.push_back(currInt);
	    }else{
	      // section too small
	      int currInt = 0;
	      output.push_back(currInt);
	    }
	  }
	  //std::cout<<integerWind<<"\n";
	}
      }
    }
    // add the start/end pounts
    // now add the start angles
    binaryFind BS;
    double sigma1;
    if(tanList[0].getZ()>=0 ){
      sigma1=1.0;
    }else{
      sigma1=-1.0;
    }
    if(turnList.size()>2){
      for(int k=1;k<=turnList.size()-2;k++){
	// potentially make new extremum of the section
	point newS1minExtremum;point newS1maxExtremum;
	if(turnList[k].first != turnList[k].second){
	  newS1minExtremum = getInterpolatedTop(turnList[k],tanList,curveList);
	}else{
	  newS1minExtremum = curveList[turnList[k].first];
	}
	if(turnList[k+1].first != turnList[k+1].second){
	  newS1maxExtremum = getInterpolatedTop(turnList[k+1],tanList,curveList);
	}else{
	  newS1maxExtremum = curveList[turnList[k+1].first];
	}
	// get the relevant subsection
	std::vector<point> sublist1(curveList.begin()+turnList[k].second,curveList.begin()+turnList[k+1].first +1);
	// add the extremal points
	checkAppendExtrema(sublist1,newS1minExtremum,newS1maxExtremum);
	double zmin = std::min(sublist1[0].getZ(),sublist1[sublist1.size()-1].getZ());
	double zmax = std::max(sublist1[0].getZ(),sublist1[sublist1.size()-1].getZ());
	//check if there is an angle
	//std::cout<<zmin<<" "<<zmax<<" "<<curveList[0].getZ()<<"\n";
	if( curveList[0].getZ() >= zmin && curveList[0].getZ() <=zmax){
	  // find the point
	  int stIndex= 0;
	  int endIndex = sublist1.size()-1;
	  double zv = curveList[0].getZ();
	  if(sublist1[0].getZ()>sublist1[endIndex].getZ()){
	    std::reverse(sublist1.begin(),sublist1.end());
	  }
	  std::pair<int,int> pair  = BS.getContainingPair(sublist1,stIndex,endIndex,zv);
	  // get the interpolated crossing
	  point interpPoint =interp(sublist1[pair.first],sublist1[pair.second],curveList[0].getZ());
	  // get the angle
	  double xdif = interpPoint.getX()-curveList[0].getX();
	  double ydif = interpPoint.getY()-curveList[0].getY();
	  double ang = angle(xdif,ydif);
	  //std::cout<<"here angle "<<ang<<"\n";
	  double sigma2;
	  if(tanList[turnList[k].second+1].getZ()>=0.0){
	    sigma2 = 1.0;
	  }else{
	    sigma2 = -1.0;
	  }
	  //check if this angle is the bottom part of a calculation or the top
	  if(sigma1>0){
	    //here it is the bottom
	    startAngles.push_back((-1.0)*sigma1*sigma2*ang);
	  }else{
	    startAngles.push_back(sigma1*sigma2*ang);
	  }
	}else{
	  startAngles.push_back(0.0);
	}
      }
      // and the end angles
      int endCIndex = curveList.size()-1;
      int endTIndex = tanList.size()-1;
      if(tanList[endTIndex].getZ()>=0 ){
	sigma1=1.0;
      }else{
	sigma1=-1.0;
      }
      for(int k=0;k<=turnList.size()-3;k++){
	// potentially make new extremum of the section
	point newS1minExtremum;point newS1maxExtremum;
	if(turnList[k].first != turnList[k].second){
	  newS1minExtremum = getInterpolatedTop(turnList[k],tanList,curveList);
	}else{
	  newS1minExtremum = curveList[turnList[k].first];
	}
	if(turnList[k+1].first != turnList[k+1].second){
	  newS1maxExtremum = getInterpolatedTop(turnList[k+1],tanList,curveList);
	}else{
	  newS1maxExtremum = curveList[turnList[k+1].first];
	}
	// get the relevant subsection
	std::vector<point> sublist1(curveList.begin()+turnList[k].second,curveList.begin()+turnList[k+1].first +1);
	// add the extremal points
	checkAppendExtrema(sublist1,newS1minExtremum,newS1maxExtremum);
	double zmin = std::min(sublist1[0].getZ(),sublist1[sublist1.size()-1].getZ());
	double zmax = std::max(sublist1[0].getZ(),sublist1[sublist1.size()-1].getZ());
	//check if there is an angle
	//std::cout<<zmin<<" "<<zmax<<" "<<curveList[endCIndex].getZ()<<"\n";
	if( curveList[endCIndex].getZ() >= zmin && curveList[endCIndex].getZ() <=zmax){
	  // find the point
	  int stIndex= 0;
	  int endIndex = sublist1.size()-1;
	  double zv = curveList[endCIndex].getZ();
	  //reverse if needs be for the inary find algo
	  if(sublist1[0].getZ()>sublist1[endIndex].getZ()){
	    std::reverse(sublist1.begin(),sublist1.end());
	  }
	  std::pair<int,int> pair  = BS.getContainingPair(sublist1,stIndex,endIndex,zv);
	  // get the interpolated crossing
	  point interpPoint =interp(sublist1[pair.first],sublist1[pair.second],zv);
	  // get the angle
	  double xdif = curveList[endCIndex].getX()-interpPoint.getX();
	  double ydif = curveList[endCIndex].getY()-interpPoint.getY();
	  double ang = angle(xdif,ydif);
	  double sigma2;
	  if(tanList[turnList[k].second+1].getZ()>=0.0){
	    sigma2 = 1.0;
	  }else{
	    sigma2 = -1.0;
	  }
	  //check if this angle is the bottom part of a calculation or the top
	  if(sigma1<0){
	    // here it is the bottom
	    endAngles.push_back((-1.0)*sigma1*sigma2*ang);
	  }else{
	    endAngles.push_back(sigma1*sigma2*ang);
	  }
	}else{
	  endAngles.push_back(0.0);
	}
      }
    }
    return output;
  }


int mutualWind2::getIntegerWindingsStar(std::vector<point> &curveList,std::vector<point> &tanList,std::vector<std::pair<int,int> > &turnList){
  if(turnList.size()>0){
    // find the maximum and minimum values amongst the turning points, we need to chekc if they are above he ends
    std::pair<double,double> extremalPair;
    extremalPair = getExtremaZ(turnList[0],tanList,curveList);
    double globalMax =extremalPair.first;
    double globalMin =extremalPair.second;
    for(int i=1;i<turnList.size();i++){
      extremalPair = getExtremaZ(turnList[i],tanList,curveList);
     if(extremalPair.first>globalMax){
	  globalMax =extremalPair.first;
       }
     if(extremalPair.second<globalMin){
	  globalMin =extremalPair.second;
       }
    }
    // now check if the global maximum is gretaer than the last point, if so extend the curve up
    int CS= curveList.size()-1;
    if(curveList[CS].getZ()<globalMax){
      double maxDif = globalMax - curveList[CS].getZ();
      int newPoints=100;
      double step = 1.5*maxDif/100.0;
      for(int i=1;i<= 100;i++){
	point p(curveList[CS].getX(),curveList[CS].getY(),curveList[CS].getZ()+i*step);
	curveList.push_back(p);
	point pt(0.0,0.0,1.0);
	tanList.push_back(pt);
      }
    }
    bool addBeggining=false;
    // now check if the global maximum is gretaer than the last point, if so extend the curve down
    if(curveList[0].getZ()>globalMin){
      double minDif = globalMin - curveList[0].getZ();
      int newPoints=100;
      double step = 1.5*minDif/100.0;
      double origZ = curveList[0].getZ();
      for(int i=1;i<= 100;i++){
	point p(curveList[0].getX(),curveList[0].getY(),origZ+ i*step);
	curveList.insert(curveList.begin(),p);
	point pt(0.0,0.0,1.0);
	tanList.insert(tanList.begin(),pt);
	addBeggining=true;
      }
    }
     //if we have added points to the beginnig of the list we must push back the turning points
     if(addBeggining){
       for(int i=0;i<turnList.size();i++){
	 turnList[i].first = turnList[i].first + 100;
	 turnList[i].second = turnList[i].second + 100;
       }
     }
  }
  // we now add the two end points to the tunring point list
  std::pair<int,int> pstart(0,0);
  std::pair<int,int> pend(curveList.size()-1,curveList.size()-1);
  turnList.insert(turnList.begin(),pstart);
  turnList.push_back(pend);
  /*We are now ready to calculate the integer widnings*/
  int integerWind=0;
  /*for(int i=0;i<turnList.size();i++){
    std::cout<<turnList[i].first<<" "<<turnList[i].second<<"\n";
    std::cout<<curveList[turnList[i].first].getZ()<<" "<<curveList[turnList[i].second].getZ()<<"\n";
    }*/
  if(turnList.size()>2){
    for(int k=0;k<=turnList.size()-3;k++){
      for(int l=k+1;l<=turnList.size()-2;l++){
	if(l==k+1){
	  integerWind = integerWind + getIntegerWindingOfPairJoined(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	}else{
	  integerWind = integerWind + getIntegerWindingOfPair(curveList,tanList,turnList[k],turnList[k+1],turnList[l],turnList[l+1]);
	}
      }
    }
  }
  else{
    integerWind=0;
  }
  //integerWind = integerWind + getIntegerWindingOfPairJoined(curveList,tanList,turnList[3],turnList[4],turnList[4],turnList[5]);
   return integerWind;
}


/*The following routine is for calculating the winding number*/

int mutualWind2::getWindingInteger(std::vector<point> &curveList1,std::vector<point> &curveList2,std::vector<int> &turnList1,std::vector<int> &turnList2){
    // find the maximum and minimum values amongst the turning points, we need to chekc if they are above he ends
    /*int tanSize1 = tanList1.size()-1;
    int tanSize2 = tanList2.size()-1;
    if(tanList1[0].getZ()>= 0 && tanList1[tanSize1].getZ()>= 0){
      extendAboveSpan(curveList1,curveList2,tanList1,tanList2,turnList1,turnList2);
    }else if(tanList1[0].getZ()< 0 && tanList1[tanSize1].getZ()<0){
      extendBelowSpan(curveList1,curveList2,tanList1,tanList2,turnList1,turnList2);
    }else if(tanList1[0].getZ()>= 0 && tanList1[tanSize1].getZ()<0){
      extendBothNegPos(curveList1,curveList2,tanList1,tanList2,turnList1,turnList2);
    }else if(tanList1[0].getZ()< 0 && tanList1[tanSize1].getZ()>= 0){
      extendBothPosNeg(curveList1,curveList2,tanList1,tanList2,turnList1,turnList2);
    }
    if(tanList2[0].getZ()>= 0 && tanList2[tanSize2].getZ()>= 0){
      extendAboveSpan(curveList2,curveList1,tanList2,tanList1,turnList2,turnList1);
    }else if(tanList2[0].getZ()< 0 && tanList2[tanSize2].getZ()<0){
      extendBelowSpan(curveList2,curveList1,tanList2,tanList1,turnList2,turnList1);
    }else if(tanList2[0].getZ()>= 0 && tanList2[tanSize2].getZ()<0){
      extendBothNegPos(curveList2,curveList1,tanList2,tanList1,turnList2,turnList1);
    }else if(tanList2[0].getZ()< 0 && tanList2[tanSize2].getZ()>= 0){
      extendBothPosNeg(curveList2,curveList1,tanList2,tanList1,turnList2,turnList1);
      }*/
  /*We are now ready to calculate the integer widnings*/
    int integerWind=0;
    for(int k=0;k<turnList1.size()-1;k++){
      for(int l=0;l<turnList2.size()-1;l++){
	integerWind = integerWind + getIntegerWindingOfPairLink(curveList1,curveList2,turnList1[k],turnList1[k+1],turnList2[l],turnList2[l+1]);
      }
    } 
    return integerWind;
  }


