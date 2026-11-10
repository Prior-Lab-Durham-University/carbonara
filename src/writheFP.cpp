#include "writheFP.h"

const double PI  =3.141592653589793238463;

writheFP::writheFP(){};

point writheFP::dbold(std::vector<point>& pointList,int size,int i,int j){
  i = i%(size-1);
  j = j%(size-1);
  return pointList[i+1].dif(pointList[j+1]);
}

point writheFP::dboldLk(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i,int j){
  i = i%(size-1);
  j = j%(size2-1);
  return pointList[i+1].dif(pointList2[j+1]);
}

point writheFP::te(std::vector<point>& pointList,int size,int i){
  point temp = dbold(pointList,size,i,i-1);
  temp.normalise();
  return temp;
}



double writheFP::acosC(double temp){
  if(temp < -1.0){
    return PI;
  }else if(temp > 1.0){
    return 0;
  }else{
    return acos(temp);
  }
}

double writheFP::mu(std::vector<point>& pointList,int size,int i,int k,int m, int j){
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

double writheFP::muLink(std::vector<point>& pointList,std::vector<point>& pointList2,int size,int size2,int i, int j){
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




double writheFP::wij(std::vector<point>& pointList,int size,int i,int j){
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


double writheFP::DI(std::vector<point>& pointList){
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

double writheFP::DIAbs(std::vector<point>& pointList){
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

double writheFP::DIClosed(std::vector<point>& pointList){
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





double writheFP::DIClosedLk(std::vector<point>& pointList,std::vector<point>& pointList2){
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


std::vector<std::pair<std::pair<int,int>,double> > writheFP::DIDownSample(std::vector<std::vector<point> >& pointListIn){
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




std::vector<std::pair<std::pair<int,int>,double> > writheFP::DIDownSampleAbs(std::vector<std::vector<point> >& pointListIn){
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

double  writheFP::DIDownSampleSingle(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  double wrval = DI(mol);
  return wrval;
}

double writheFP::DIDownSampleAbsSingle(std::vector<std::vector<point> >& pointListIn){
  // first downsample
  std::vector<point> mol;
  for(int i=0;i<pointListIn.size();i++){
    mol.push_back(pointListIn[i][0]);
  }
  std::vector<point> lastSec = pointListIn[pointListIn.size()-1];
  mol.push_back(lastSec[lastSec.size()-1]);
  double wrval = DIAbs(mol);
  return wrval;
}

double writheFP::DIDownSampleAbsSingleSKMT(ktlMolecule&molin){
  // first downsample
  skmt s;
  std::vector<point> mol = s.getSKMTCurve(molin);
  double wrval = DIAbs(mol);
  return wrval;
}

std::pair<double,double> writheFP::sumWr(int& iv,int& jv,std::vector<std::vector<double> >& wijVals){
  double sigsum=0.0;
  double sigsumAbs=0.0;
  for(int i=iv;i<jv-1;i++){
    for(int j=iv+1;j<jv;j++){
      sigsum = sigsum +wijVals[i][j];
      sigsumAbs = sigsumAbs + abs(wijVals[i][j]);
     }
  }
  std::pair<double,double> p(sigsum/(2.0*PI),sigsumAbs/(2.0*PI));
  return p;
}

std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > >  writheFP::getWritheFingerprints(std::vector<point>& pointList){
  int listSize = pointList.size();
  std::vector<std::vector<double> > wijVals(listSize,std::vector<double>(listSize,0.0));
  // calculate all the winding densities
  for(int i=0;i<listSize-1;i++){
    for(int j=i+1;j< listSize;j++){
      wijVals[i][j] =wij(pointList,listSize,i,j);
     }
  }
  std::vector<std::vector<point> > wr;
  std::vector<std::vector<point> > absWr;
  std::vector<point> wrRow;
  std::vector<point> absWrRow;
  for(int i=0;i<listSize;i++){
    for(int j=i+1;j<listSize;j++){
      std::pair<double,double> wrvals =sumWr(i,j,wijVals);
      point p(i,j,wrvals.first);
      point pabs(i,j,wrvals.second);
      wrRow.push_back(p);
      absWrRow.push_back(pabs);
    }
    wr.push_back(wrRow);
    absWr.push_back(absWrRow);
    wrRow.clear();
    absWrRow.clear();
  }
  std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > > outpair;
  outpair.first = wr;
  outpair.second = absWr;
  return outpair;
}



void writheFP::DIDownSampleWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
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

void writheFP::DIDownSampleAbsWrite(std::vector<std::vector<point> >& pointListIn,const char* filename){
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

double writheFP::compare(std::vector<point>& sec1,std::vector<point>& sec2,int& maxSz){
  double sum=0.0;
  for(int i=0;i<sec1.size();i++){
    double denomW = 0.2*(sec1[i].getY()-sec1[i].getX());
    sum = sum + (std::abs(sec1[i].getZ()-sec2[i].getZ()))/denomW;
  }
  return sum/sec1.size();
}



bool compareOnSecond(std::pair<std::pair<int,int>,std::pair<int,int> > p1,std::pair<std::pair<int,int>,std::pair<int,int> > p2)
{
    return (p1.second.first < p2.second.first);
}



bool compareOnFirst(std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > p1,std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > >  p2)
{
    return (p1[0].first.first < p2[0].first.first);
}

bool compareOnFirst2(std::pair<std::pair<int,int>,std::pair<int,int> >  p1,std::pair<std::pair<int,int>,std::pair<int,int> >   p2)
{
    return (p1.first.first < p2.first.first);
}

bool writheFP::checkOverlap(std::pair<std::pair<int,int>,std::pair<int,int> >& p1,std::pair<std::pair<int,int>,std::pair<int,int> >& p2){
  if(((p1.first.first>=p2.first.first && p1.first.second<=p2.first.second) ||(p2.first.first>=p1.first.first && p2.first.second<=p1.first.second)||(p1.first.first>=p2.first.first && p1.first.second<=p2.first.second))&&((p1.second.first>=p2.second.first && p1.second.second<=p2.second.second) ||(p2.second.first>=p1.second.first && p2.second.second<=p1.second.second)||(p1.second.first>=p2.second.first && p1.second.second<=p2.second.second)) ){
    return true;
  }else{
    return false;
  }
}

bool writheFP::checkOverlap(std::pair<int,int>& p1,std::pair<int,int>& p2){
  std::pair<bool,int> out;
  if((p1.first>=p2.first && p1.second<=p2.second) ||(p2.first>=p1.first && p2.second<=p1.second) ){
     return true;
  }else{
    return false;
  }
}


std::pair<bool,double>  writheFP::checkOverlapLargest(std::pair<int,int>& p1,std::pair<int,int>& p2){
  std::pair<bool,int> out;
  if((p1.first<=p2.first && p1.second >= p2.first)||(p1.first>=p2.first && p1.first <=p2.second)|| (p1.first>=p2.first && p1.second<=p2.second) ||(p2.first>=p1.first && p2.second<=p1.second) ){
    // choose the largest one
    out.first = true;
    if((p1.second-p1.first) >(p2.second-p2.first)){ 
      out.second=p1.second-p1.first;
      return out;
    }else{
      out.second=p2.second-p2.first;
      return out;
    }
  }else{
    out.first=false;
    out.second=p2.second-p2.first+p1.second-p1.first;
    return out;
  }
}


std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > writheFP::compareFingerPrints(std::vector<std::vector<point> >& fp1In,int& len1In,std::vector<std::vector<point> >& fp2In,int& len2In,double& cutOff){
  int maxZ = 8;
  std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > allPairs;
  std::vector<std::vector<point> > fp1;std::vector<std::vector<point> > fp2;
  int len1,len2;
  if(len2In>len1In){
    fp1=fp1In;
    fp2=fp2In;
    len1 =len1In;
    len2 = len2In;
  }else{
    fp1 = fp2In;
    fp2 = fp1In;
    len2 =len1In;
    len1 = len2In;
  }
    // loop over all sections of a length from 8 to the length of the smaller first molecule
  std::pair<std::pair<int,int>,std::pair<int,int> > secPr;
  std::pair<int,int> sec1;
  std::pair<int,int> sec2;
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > pairSet;
  for(int i=10;i<len1;i++){
    //loop over all sections of molecule one of the given length
    int lenDif = len1-i;
    for(int j=0;j<lenDif;j++){
      // select the required section of molecule 1
      std::vector<point> sub1 = {fp1[j].begin(),fp1[j].begin()+i};
      int lenDif2 =len2-i;
      for(int k=0;k<lenDif2;k++){
	std::vector<point> sub2 = {fp2[k].begin(),fp2[k].begin()+i};
	double compVal = compare(sub1,sub2,maxZ);
	if(compVal<cutOff){
	  sec1.first=j;sec1.second=j+i;
	  sec2.first=k;sec2.second=k+i;
	  secPr.first = sec1;
	  secPr.second = sec2;
	  pairSet.push_back(secPr);
	  //std::cout<<"whu ? "<<j<<" "<<j+i<<" "<<k<<" "<<k+i<<" "<<compVal<<"\n";
	}
      }
    }
    if(pairSet.size()>0){
      allPairs.push_back(pairSet);
      pairSet.clear();
    }
  }
  // now find the largest overlap and check which sets are disjoint with it, this is what allPairsLeftOver has
  if(allPairs.size()>0){
    std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > allPairsFinal;
    std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > allPairsLeftOver;
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > pairSetBiggest = allPairs[allPairs.size()-1];
    allPairsFinal.push_back(pairSetBiggest);
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > newLayer;
    for(int i=allPairs.size()-2;i>-1;i--){
      for(int j=0;j<allPairs[i].size();j++){
	bool pairCheck=false;
	for(int k=0;k<pairSetBiggest.size();k++){
	  if(checkOverlap(allPairs[i][j].first,pairSetBiggest[k].first)==false){
	    pairCheck=true;
	  }
	}
	if(pairCheck){
	  newLayer.push_back(allPairs[i][j]);
	}
      }
      if(newLayer.size()>0){
	allPairsLeftOver.push_back(newLayer);
	newLayer.clear();
      }
    }

    std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > allPairsLeftOverNew;
    int maxIterations = allPairsLeftOver.size()+2;
    // check for non disjoint sets with this largest element
    int l=0;
    while(allPairsLeftOver.size()>0 && l<maxIterations){
      l++;
      pairSetBiggest = allPairsLeftOver[0];
      allPairsFinal.push_back(allPairsLeftOver[0]);
      if(allPairsLeftOver.size()>1){
	//std::cout<<" how many left over ?"<<allPairsLeftOver.size()<<"\n";
	for(int i=1;i<allPairsLeftOver.size()-1;i++){
	  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > newLayer;
	  for(int j=0;j<allPairsLeftOver[i].size();j++){
	    bool isOkay=false;
	    for(int k=0;k<pairSetBiggest.size();k++){
	      if((allPairsLeftOver[i][j].first.first<=pairSetBiggest[k].first.first &&allPairsLeftOver[i][j].first.second<=pairSetBiggest[k].first.first) ||(allPairsLeftOver[i][j].first.first>=pairSetBiggest[k].first.second && allPairsLeftOver[i][j].first.second>=pairSetBiggest[k].first.second)){
		isOkay=true;
	      }
	    }
	    if(isOkay){
	      newLayer.push_back(allPairsLeftOver[i][j]);
	    }
	  }
	  if(newLayer.size()>0){
	    allPairsLeftOverNew.push_back(newLayer);
	  }
	}
        allPairsLeftOver=allPairsLeftOverNew;
	allPairsLeftOverNew.clear();
      }else{
	l=maxIterations+10;
      }
    }
    // now sort on the second molecule
    // first we sort the sets of the first molecule as we want an ordered comparison
    std::vector< std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > nonRepeats;
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > >  nonRepeatSlice;
    if(allPairsFinal.size()>0){
      std::sort(allPairsFinal.begin(),allPairsFinal.end(),compareOnFirst);
      for(int i=0;i<allPairsFinal.size();i++){
	bool hasPushedBack=false;
	for(int j=0;j< allPairsFinal[i].size()-1;j++){
	  bool isCool =true;
	  for(int k=j+1;k<allPairsFinal[i].size();k++){
	    if(checkOverlap(allPairsFinal[i][j],allPairsFinal[i][k])==true){
	      isCool= false;
	    }
	  }
	  if(isCool==true){
	    nonRepeatSlice.push_back(allPairsFinal[i][j]);
	    if(j==allPairsFinal[i].size()-2&&allPairsFinal[i].size()>1){
	      nonRepeatSlice.push_back(allPairsFinal[i][allPairsFinal[i].size()-1]);
	    }
	    hasPushedBack=true;
	  }
	}
	if(hasPushedBack==false){
	  nonRepeatSlice.push_back(allPairsFinal[i][0]);
	}
	nonRepeats.push_back(nonRepeatSlice);
	nonRepeatSlice.clear();
      }
    }
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > fSet = finalSort(nonRepeats);
    // ensure the output ordering matches the users order
    if(len2In>len1In){
      return fSet;
    }else{
      //swap back
      for(int i=0;i<fSet.size();i++){
	std::swap(fSet[i].first,fSet[i].second);
      }
      return fSet;
    }
  }else{
    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > out;
    std::pair<int,int> p1(0,0);
    std::pair<std::pair<int,int>,std::pair<int,int> > pr(p1,p1);
    out.push_back(pr);
    return out;
  }
}
  

std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > writheFP::finalSort(std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > & prSet){
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > >  bestCombo;
  if(prSet.size()>0){
    int noSlices=prSet.size();
    std::vector<int> permSizes;
    int permtot=1.0;
    for(int i=0;i<noSlices;i++){
      permSizes.push_back(prSet[i].size());
      permtot=permtot*prSet[i].size();
    }
    if(noSlices>1 &&permtot<500){
      std::vector<int> permSizesSummed;
      int runningTot=0;
      for(int i=0;i<permSizes.size();i++){
	if(i==0){
	  runningTot=1;
	}else{
	  runningTot=runningTot*permSizes[i-1];
	}
	permSizesSummed.push_back(runningTot);
      }
      std::reverse(permSizesSummed.begin(),permSizesSummed.end());
      std::vector<int> perm (noSlices,0);
      std::vector<int> permTest (noSlices,0);
      std::vector<std::vector<int> > permutations;
      for(int i=1;i<=runningTot*permSizes[permSizes.size()-1];i++){
	for(int j=0;j<permSizesSummed.size();j++){
	  if(i%permSizesSummed[j]==0){
	    permTest[j]=permTest[j]+1;
	    perm[j] = permTest[j]%permSizes[j];
	  }
	}
	permutations.push_back(perm);
      }
      //now test the permutations for the longest set;
      int Maxlen =0;
      int bestIndex=0;
      std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > bestSet;
      for(int i=0;i<permutations.size();i++){
	std::pair<int,std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > comb =findBiggestOverlap(prSet,permutations[i]);
	if(comb.first>Maxlen){
	  Maxlen =comb.first;
	  bestSet = comb.second;
	}
      }
      //now grab the best combination
      bestCombo.insert( bestCombo.end(), bestSet.begin(), bestSet.end() );
    }else{
      //just choose one form the list
      bestCombo.push_back(prSet[0][0]);
    }
  }
  return bestCombo;
}




std::pair<int,std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > writheFP::findBiggestOverlap(std::vector<std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > >& separatedCombo,std::vector<int>& perm){
  int outLen=0;
  //generate possible combinations
  std::vector<int> slices;
  for(int i=0;i<perm.size();i++){
    slices.push_back(i);
  }
  // find all apirs which cannot overlap
  std::vector<std::pair<int,int> > forbiddenPairs;
  for(int i=0;i<perm.size()-1;i++){
    for(int j=i+1;j<perm.size();j++){
      std::pair<bool,double> chk1 = checkOverlapLargest(separatedCombo[i][perm[i]].second,separatedCombo[j][perm[j]].second);
      std::pair<bool,double> chk2 = checkOverlapLargest(separatedCombo[i][perm[i]].first,separatedCombo[j][perm[j]].first);
      if(chk1.first==true || chk2.first==true){
	std::pair<int,int> pr(i,j);
	forbiddenPairs.push_back(pr);
      }
    }
  }
  // generate all possible combinations
  std::vector<std::vector<int> > sets = getAllSubsets(slices);

  // isolate the combinations which are valid
  std::vector<std::vector<int> > validSets;
  bool isSetOkay;
  for(int i = 0;i<sets.size();i++){
    isSetOkay=true;
    if(sets[i].size()>1){
      for(int k=0;k<forbiddenPairs.size();k++){
	if((std::find(sets[i].begin(),sets[i].end(),forbiddenPairs[k].first)!=sets[i].end())&&(std::find(sets[i].begin(),sets[i].end(),forbiddenPairs[k].second)!=sets[i].end())){
	      isSetOkay=false;
	}
      }
    }
    if(isSetOkay){
      validSets.push_back(sets[i]);
     }
  }
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > bestCombo;
  std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > tempCombo;
  int maxLen=0;
  int bestSet=0;
  for(int i=0;i<validSets.size();i++){
    int currLen = 0;
    for(int j=0;j<validSets[i].size();j++){
      currLen = currLen +separatedCombo[validSets[i][j]][perm[validSets[i][j]]].second.second- separatedCombo[validSets[i][j]][perm[validSets[i][j]]].second.first;
      tempCombo.push_back(separatedCombo[validSets[i][j]][perm[validSets[i][j]]]);
    }
    if(i==0){
      bestCombo=tempCombo;
    }
    if(currLen>maxLen){
      maxLen= currLen;
      bestSet =i;
      bestCombo=tempCombo;
    }
    tempCombo.clear();
  }
  std::pair<int,std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > > output;
  output.first =maxLen;
  output.second =bestCombo;
  return output;
}

std::vector< std::vector<int> > writheFP::getAllSubsets(std::vector<int>& set)
{
  std::vector< std::vector<int> > subset;
    std::vector<int> empty;
    subset.push_back( empty );

    for (int i = 0; i < set.size(); i++)
    {
        std::vector< std::vector<int> > subsetTemp = subset;  //making a copy of given 2-d vector.

        for (int j = 0; j < subsetTemp.size(); j++)
            subsetTemp[j].push_back( set[i] );   // adding set[i] element to each subset of subsetTemp. like adding {2}(in 2nd iteration  to {{},{1}} which gives {{2},{1,2}}.

        for (int j = 0; j < subsetTemp.size(); j++)
            subset.push_back( subsetTemp[j] );  //now adding modified subsetTemp to original subset (before{{},{1}} , after{{},{1},{2},{1,2}}) 
    }
    return subset;
}