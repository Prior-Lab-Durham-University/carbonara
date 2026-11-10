# include "binaryFind.h"

binaryFind::binaryFind(){
	 isIn =false;
}

void binaryFind::checkInSec(double lowerSec, double upperSec,double& zval){
	if(zval <= upperSec && zval >= lowerSec){
		isIn =true;
	}
}

bool binaryFind::checkInSecInitial(double lowerSec, double upperSec,double& zval){
	bool isInC;
	if(zval <= upperSec && zval >= lowerSec){
		isInC =true;
	}
	return isInC;
}

std::pair<int,int>  binaryFind::getContainingPair(std::vector<point>& pointList,int& lowerIn,int& upperIn,double& zv){
	lowerIndex = lowerIn;
	upperIndex = upperIn;
	zval = zv;
	int indexDif;
	bool switchPoinstAtEnd = false;
	if(upperIndex>(lowerIndex+1)){
	  if(lowerIndex>upperIndex){
	    switchPoinstAtEnd = true; 
	    int lowerIndexTemp = lowerIndex;
	    lowerIndex = upperIndex;
	    upperIndex =lowerIndexTemp;
	  }
	  indexDif= upperIndex-lowerIndex;
	  int midIndex;
	  if(indexDif==2){
	    midIndex = lowerIndex+1;
	  }
	  while(indexDif>2){
	    if(indexDif%2 ==0){
	      midIndex = lowerIndex + indexDif/2;
	    }else{
	      midIndex = lowerIndex+ ceil(short(indexDif)/2.0);
	    }
	    //std::cout<<pointList[midIndex].getZ()<<" "<<pointList[upperIndex].getZ()<<"\n";
	    if(switchPoinstAtEnd){
	      checkInSec(pointList[upperIndex].getZ(),pointList[midIndex].getZ(),zval);
	    }else{
	      checkInSec(pointList[midIndex].getZ(),pointList[upperIndex].getZ(),zval);
	    }
	    if(isIn){
	      lowerIndex = midIndex;	
	    }else{
	      upperIndex = midIndex;
	    }
	    isIn=false;
	    indexDif = upperIndex-lowerIndex;
	  }
	  //after this look we should have narrowed to 3 indicies, sometimes the mid index will be equal to one or the other
	  if(midIndex == lowerIndex){
	    midIndex++; 
	  }
	  if(midIndex == upperIndex){
	    midIndex--; 
	  }
	  if(switchPoinstAtEnd){
	    checkInSec(pointList[upperIndex].getZ(),pointList[midIndex].getZ(),zval);
	    if(isIn){
	      boundingIndicies.first =upperIndex;
	      boundingIndicies.second = midIndex;
	    }else{
	      boundingIndicies.first = midIndex;
	      boundingIndicies.second =lowerIndex;
	    }	
	  }else{
	    checkInSec(pointList[midIndex].getZ(),pointList[upperIndex].getZ(),zval);
	    if(isIn){
	      boundingIndicies.first = midIndex;
	      boundingIndicies.second = upperIndex;
	    }else{
	      boundingIndicies.first = lowerIndex;
	      boundingIndicies.second =midIndex;
	    }	
	  }
	}else{
	  boundingIndicies.first=lowerIndex;
	  boundingIndicies.second=upperIndex;
	}
	// reset the isIn checker
	isIn=false;
	return boundingIndicies;
	}
	
