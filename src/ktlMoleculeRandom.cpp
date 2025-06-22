#include "ktlMoleculeRandom.h"
#include "point.h"

ktlMolecule::ktlMolecule(){
  kapvallink = 0.36932;
  kapvalbeta = 0.343782;
  kapvalalpha =0.380928;
  tauvallink = 0.13439;
  tauvalbeta = 0.123222;
  tauvalalpha = 0.145831;
  alvallink = 4.22128;
  alvalbeta = 4.13944;
  alvalalpha = 4.2329;
};

void ktlMolecule::setParams(double &rminIn,double &rmaxIn,double &lminIn){
 rmg.setParams(rminIn,rmaxIn,lminIn);
}

std::vector<int> ktlMolecule::getUnitNos(){
  return noPts;
}

std::vector<std::vector<point> > ktlMolecule::getTangents(){
  return tanlist;
};

std::vector<std::vector<point> > ktlMolecule::getNormals(){
  return normlist;
};

std::vector<std::vector<point> > ktlMolecule::getBinormals(){
  return binormlist;
};

std::vector<std::vector<point> > ktlMolecule::getCoordinates(){
  return coords;
};

std::vector<point> ktlMolecule::getCoordinatesSection(int i){
  return coords[i];
}

std::vector<std::vector<std::string> > ktlMolecule::getAminoList(){
  return aminoList;
}

int ktlMolecule::getSubsecSize(int sec){
  // the minus 1 is becasue the labellingnumbers are 1,2,3 e.t.c but for chainList 1 is at index0
  return chainList[sec-1].second-chainList[sec-1].first+1;
}

 std::vector<std::vector<point> > ktlMolecule::getSubsecCoordinates(int &sec){
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  return subcoords;
}

std::vector<std::pair<std::string,int> > ktlMolecule::getNameSizeListOfSection(int &sec){
  std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[sec].first;
  std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[sec].second+1;
  std::vector<std::pair<std::string,int> > subns(first,second);
  return subns;
}

double ktlMolecule::maxNeighbourDistSec(int &sec){
  // the minus 1 is becasue the labellingnumbers are 1,2,3 e.t.c but for chainList 1 is at index0
  double dmax=0.0;
  for(int i = chainList[sec-1].first;i<= chainList[sec-1].second;i++){
    for(int j=0;j<coords[i].size();j++){      double d;
      if(j==(coords[i].size()-1) && i<chainList[sec-1].second){
	d = coords[i][j].eDist(coords[i+1][0]);
      }else if(j==(coords[i].size()-1) && i==chainList[sec-1].second){
	d=0.0;
      }else{
	d = coords[i][j].eDist(coords[i][j+1]);
      }
      if(d>dmax){
	dmax=d;
	if(d>4.0){
	  //std::cout<<i<<" "<<j<<"\n";
	}
      }
    }
  }
  return dmax;
}

/*std::vector<double> ktlMolecule::getDistChanges(){
  return distChanges;
  }*/


double ktlMolecule::getCurvatureJoined(int index){
  double kapval;
  int sz =coords[index].size();
  if(sz<=2){
    kapval = kapvallink;
  }else if(3 < sz && sz < 7){
    kapval = kapvalbeta;
  }else{
    kapval = kapvalalpha;
  }
  return kapval;
}


double ktlMolecule::getTorsionJoined(int index){
  double tauval;
  int sz =coords[index].size();
  if(sz<=2){
    tauval = tauvallink;
  }else if(3 < sz && sz < 7){
    tauval = tauvalbeta;
  }else{
    tauval = tauvalalpha;
  }
  return tauval;
}


 int ktlMolecule::getUnitNo(int index){
   return coords[index].size();
}

point ktlMolecule::getTangent(int mindex,int subindex){
  return tanlist[mindex][subindex];
};

point ktlMolecule::getNormal(int mindex,int subindex){
  return normlist[mindex][subindex];
};

point ktlMolecule::getBinormal(int mindex,int subindex){
  return binormlist[mindex][subindex];
};

point ktlMolecule::getCoordinate(int mindex,int subindex){
  if(subindex==-1){
    return coords[mindex-1][coords[mindex-1].size()-1];
  }else{
    return coords[mindex][subindex];
  }
};


double ktlMolecule::getAlbeadJoined(int index){
  double length;
  int sz =coords[index].size();
  if(sz<=2){
    length = alvallink;
  }else if(3 < sz && sz < 7){
    length = alvalbeta;
  }else{
    length = alvalalpha;
  }
  return length;
}

std::string ktlMolecule::getType(int &index){
  return nameSizeList[index].first;
}

std::string ktlMolecule::getType(int &chainNo,int &index){
  int fullIndex = chainList[chainNo-1].first;
  return nameSizeList[fullIndex+index].first;
}


/*double ktlMolecule::getDistChange(int index){
  return distChanges[index];
  }*/

double ktlMolecule::getMaxDistChange(){
  return maxDistChange;
}

/*double ktlMolecule::getMaxScatLength(){
  double maxLen =0.0;
  double len=0.0;
  for(int i=0;i<coords.size();i++){
    for(int j=i+1;j<coords.size();j++){
      len = coords[i].eDist(coords[j]);
      if(len>maxLen){
        maxLen = len;
      }
    }
  }
  return 2.0*maxLen;
  }*/

int ktlMolecule::noSecSize(){
  return coords.size();
}

int ktlMolecule::noChains(){
  return chainList.size();
}

int ktlMolecule::molSize(){
  int sz = 0;
  for(int i=0;i<coords.size();i++){
    if(nameSizeList[i].first=="Helix" || nameSizeList[i].first=="Strand"){
      sz = sz + 1;
    }else{
      sz = sz + coords[i].size();
    }
  }
  return sz;
}

int ktlMolecule::getNoAminos(){
  int sz = 0;
  for(int i=0;i<coords.size();i++){
    sz=sz+coords[i].size();
  }
  return sz;
}

std::pair<double,double> ktlMolecule::getMaxPossibleLength(){
  double maxL=0.0;
  double totLength=0.0;
  double dst;
  for(int i=0;i<coords.size()-1;i++){
    for(int j=i+1;j<coords.size();j++){
      dst =coords[i][0].eDist(coords[j][0]);
      if(j==i+1){
	totLength= totLength+dst;
      }
      if(dst>maxL){
	maxL = dst;
      }
    }
  }
  maxL = maxL*2.0;
  std::pair<double,double> mxpr;
  mxpr.first = totLength;
  mxpr.second = maxL;
  return mxpr;
}








void ktlMolecule::readInSequence(const char* filename,double &rmin,double &rmax,double &lmin){
  int npts;
  std::ifstream myfile;
  myfile.open(filename);
  std::string output;
  double val,prevval,X,Y,Z;
  int power;
  int n =0;
  bool isFirst=true;
  int nunits;
  int nosecs=0;
  std::string sequence;
  std::string predictions;
  std::string empty;
  std::string chainNo;
  std::string type,prevType;
  std::vector<std::string> aminoType;
  int noChains;
  if(myfile.is_open()) {
    // read in the first line which tells us howmany chains there are
    std::getline(myfile,chainNo);
    std::stringstream ss(chainNo);
    ss>>noChains;
    std::cout<<"number of chains "<<noChains<<"\n";
    distSetsSecs.resize(noChains);
    for(int i=1;i<=noChains;i++){
      std::pair<int,int> p;
      p.first=nameSizeList.size();
      // grab the amino sequence (currently not used)
      std::getline(myfile,empty);
      std::getline(myfile,sequence);
      //std::cout<<sequence<<"\n";
      // grab the predictions
      std::getline(myfile,empty);
      std::getline(myfile,predictions);
      //std::cout<<predictions<<"\n";
      // input should read as a letter (amino type, currently ignore) and then the type H (helix) S (strand) or C (loop)
      n=0;
      for(int i=0;i<predictions.size();i++){
	if(i==0){
	  if(predictions.compare(i,1,"-")==0){
	    type = "loop";
	  }else if(predictions.compare(i,1,"E")==0 || predictions.compare(i,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "Helix";
	  }
	  n=1;
	  prevType=type;
	  aminoType.push_back(sequence.substr(i,1));
	}else{
	  if(predictions.compare(i,1,"-")==0){
	    type = "loop";
	  }else if(predictions.compare(i,1,"E")==0 || predictions.compare(i,1,"S")==0){
	    type = "Strand";
	  }else{
	    type = "Helix";
	  }
	  if(prevType==type){
	    // keep increasing the section size
	    n++;
	    aminoType.push_back(sequence.substr(i,1));
	  }else{
	    //store the previous type and length, if it is a helix check for breaks
	      //strand or linker, no need to separate
	      std::pair<std::string,int> stpr;
	      stpr.first = prevType;stpr.second = n;
	      nameSizeList.push_back(stpr);
	      distChanges.push_back(0.0);
	      prevType = type;
	      aminoList.push_back(aminoType);
	      aminoType.clear();
	      // empty the aminoVector
	      aminoType.push_back(sequence.substr(i,1));
	      n=1;
	  }
	}
      }
      std::pair<std::string,int> stpr;
      stpr.first = prevType;stpr.second = n;
      nameSizeList.push_back(stpr);
      distChanges.push_back(0.0);
      aminoList.push_back(aminoType);
      aminoType.clear();
      p.second = nameSizeList.size()-1;
      chainList.push_back(p);
      /*std::cout<<chainList.size()<<"\n";
      for(int iv=chainList[chainList.size()-2].second+1;iv<nameSizeList.size();iv++){
	std::cout<<nameSizeList[iv].first<<" "<<nameSizeList[iv].second<<"\n";
	}*/
    }
  //don't forget to set the random moelcule parameets
  rmg.setParams(rmin,rmax,lmin);
}else{
    // here no sequence secondary predictions provided, will just assume one large linker
    std::cout<<"no sequence provided for chain number "<<chainNo<<"\n";
  }
}

void ktlMolecule::readInCoordinates(const char* filename){
  std::ifstream coordFile(filename);
  std::string line;
  //std::cout<<"here?\n";
  if(coordFile.is_open()){
    std::vector<point> secondarySec;
    for(int iv=0;iv<nameSizeList.size();iv++){
      //keep reading until we get a coordinate
      bool waitToFill=true;
      int noInSection = nameSizeList[iv].second;
      int currNoMols=0;
      //std::cout<<"section number "<<iv<<" "<<nameSizeList[iv].first<<" "<<noInSection<<" "<<nameSizeList.size()<<"\n";
      while((coordFile.eof()==false)&&waitToFill==true){
	std::getline(coordFile,line);
	if(line.find("chain")==std::string::npos){
	  if(line.size()>1){
	    point p(line);
	    currNoMols++;
	    secondarySec.push_back(p);
	    if(currNoMols==noInSection){
	      // std::cout<<" this big pushed back"<<secondarySec.size()<<"\n";
	      coords.push_back(secondarySec);
	      secondarySec.clear();
	      waitToFill=false;
	      //std::cout<<"here?\n";
	    }
	  }
	}
      }
      //std::cout<<coords.size()<<" "<<coords[iv].size()<<" "<<nameSizeList.size()<<" "<<nameSizeList[iv].second<<"\n";
    }
  }else{
    std::cout<<"there is no valid coordinate file supplied\n";
  }
}


std::vector<double> ktlMolecule::getHydrophobicDistance(std::vector<std::vector<point> > &solventList,double &maxSolDist){
  // maxSolDist is the size of the spehere I check for solvents in
  std::vector<double> meanHydroAminoDists;
  for(int i=0;i<hydroPhobicList.size();i++){
    std::pair<int,int> pr = hydroPhobicList[i];
    point coord = coords[pr.first][pr.second];
    double solDist = 0.0;
    int noDistances=0;
    // loop over all solvents, I don't care where they come from.
    for(int k=0;k<solventList.size();k++){
      for(int l=0;l<solventList[k].size();l++){
	double sddist = coord.eDist(solventList[k][l]);
	if(sddist<maxSolDist){
	  solDist = solDist + sddist;
	  noDistances++;
	}
      }
    }
    if(noDistances>0){
      solDist = solDist/double(noDistances);
    }
    //std::cout<<solDist<<" "<<noDistances<<"\n";
    meanHydroAminoDists.push_back(double(noDistances)/double(hydroPhobicList.size()));
  }
  return meanHydroAminoDists;
}


void ktlMolecule::getHydrophobicResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="A"||aminoList[i][j]=="I"||aminoList[i][j]=="L"||aminoList[i][j]=="M"||aminoList[i][j]=="F"||aminoList[i][j]=="V"||aminoList[i][j]=="P"||aminoList[i][j]=="G"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	hydroPhobicList.push_back(pr);
      }
    }
  }
}

void ktlMolecule::getCoiledCoilResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="L"||aminoList[i][j]=="I"||aminoList[i][j]=="G"){
	//check if the amino belongs to a helix
	if(nameSizeList[i].first=="Helix"){
	  std::pair<int,int> pr;
	  pr.first = i;pr.second = j;
	  coiledCoilList.push_back(pr);
	}
      }
    }
  }
  /* std::cout<<"in coiled coild routine "<<coiledCoilList.size()<<"\n";
  for(int i=0;i<coiledCoilList.size();i++){
    std::pair<int,int> pr = coiledCoilList[i];
    std::cout<<"in section "<<pr.first<<" of type "<<nameSizeList[pr.first].first<<" at residue "<<pr.second<<"\n";
    }*/
}

double ktlMolecule::coiledCoilPotential(){
  double potentialOut=0.0;
  int noCalcs=0;
  for(int i=0;i<coiledCoilList.size()-1;i++){
    for(int j=i+1;j<coiledCoilList.size();j++){
      if(coiledCoilList[i].first!=coiledCoilList[j].first){
	point p1=coords[coiledCoilList[i].first][coiledCoilList[i].second];
	point p2=coords[coiledCoilList[j].first][coiledCoilList[j].second];
	double d = p1.eDist(p2);d=d-7.3;
  double dsix = d*d;
	potentialOut = potentialOut + dsix;
       noCalcs++;
      }
    }
  }
  if(noCalcs>1){
  potentialOut= potentialOut/double(noCalcs);
  }
  return potentialOut;
}


double ktlMolecule::coiledCoilPotentialBetween(int &secNo){
  double potentialOut=0.0;
  int noCalcs=0;
  int firstIndex =chainList[secNo].first;
  int secondIndex = chainList[secNo].second;
  // get the sublist of hydrophobic residues in this section
  std::vector<std::pair<int,int> > coiledCoilSubList;
  for(int i=0;i<coiledCoilList.size();i++){
      if((coiledCoilList[i].first>=firstIndex)&&(coiledCoilList[i].first<=secondIndex)){
        coiledCoilSubList.push_back(coiledCoilList[i]);
      }
  }
  std::vector<std::pair<int,int> > coiledCoilNotSubList;
  for(int i=0;i<coiledCoilList.size();i++){
      if((coiledCoilList[i].first<firstIndex)||(coiledCoilList[i].first>secondIndex)){
        coiledCoilNotSubList.push_back(coiledCoilList[i]);
      }
  }
  for(int i=0;i<coiledCoilSubList.size();i++){
    for(int j=0;j<coiledCoilNotSubList.size();j++){
      point p1=coords[coiledCoilList[i].first][coiledCoilSubList[i].second];
      point p2=coords[coiledCoilNotSubList[j].first][coiledCoilNotSubList[j].second];
      double d = p1.eDist(p2);d=d-7.3;
      double dsix = d*d;
      potentialOut = potentialOut + dsix;
      noCalcs++;
    }
  }
  if(noCalcs>1){
    potentialOut= potentialOut/double(noCalcs);
  }
  return potentialOut;
}

double ktlMolecule::coiledCoilPotentialBetween(){
  double coiledCoilPotential=0.0;
  int noChains = chainList.size();
  //std::cout<<"no chains are ? "<<noChains<<"\n";
  for(int i=0;i<noChains;i++){
    coiledCoilPotential = coiledCoilPotential + 0.5*coiledCoilPotentialBetween(i);
  }
  return coiledCoilPotential;
}


void ktlMolecule::getPolarResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="Q"||aminoList[i][j]=="N"||aminoList[i][j]=="H"||aminoList[i][j]=="S"||aminoList[i][j]=="T"||aminoList[i][j]=="Y"||aminoList[i][j]=="C"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	polarList.push_back(pr);
      }
    }
  }
}

void ktlMolecule::getPositiveResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="K"||aminoList[i][j]=="Y"||aminoList[i][j]=="R"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	posChargeList.push_back(pr);
      }
    }
  }
}


void ktlMolecule::getNegativeResidues(){
  for(int i=0;i<aminoList.size();i++){
    for(int j=0;j<aminoList[i].size();j++){
      if(aminoList[i][j]=="D"||aminoList[i][j]=="E"){
	std::pair<int,int> pr;
	pr.first = i;pr.second = j;
	negChargeList.push_back(pr);
      }
    }
  }
}



point ktlMolecule::getCentreOfMass(std::vector<std::vector<point> > &cdSet){
  point mean(0.0,0.0,0.0);
  int noMols=0;
  for(int i=0;i<cdSet.size();i++){
    for(int j=0;j<cdSet[i].size();j++){
    mean = mean + cdSet[i][j];
    noMols=noMols+1;
    }
  }
  return mean*(1.0/double(noMols));
}

std::vector<int> ktlMolecule::checkOverlap(std::vector<std::vector<point> > &cdsIN){
  std::vector<int> overlappedSecs;
  for(int i=0;i<cdsIN.size()-1;i++){
    bool triggered =false;
    for(int j=i+1;j<cdsIN.size();j++){
      for(int k=0;k<cdsIN[i].size();k++){
	for(int l=0;l<cdsIN[j].size();l++){
	  double dist= cdsIN[i][k].eDist(cdsIN[j][l]);
	  if(k==(cdsIN[i].size()-1)&&l==0 && j==i+1){
	    dist=10.0;
	  }
	  if(dist<5.0){
	    triggered = true;
	  }
	}
      }
    }
    if(triggered==true){
      overlappedSecs.push_back(i);
    }
  }
  return overlappedSecs;
}




std::vector<double> ktlMolecule::checkOverlapWithRad(double &wRad,int &sec){
  std::vector<double> overlappedSecs;
  distSetsSecs[sec-1].clear();
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec-1].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec-1].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  for(int i=0;i< subcoords.size()-1;i++){
    bool triggered =false;
    for(int j=i+1;j<subcoords.size();j++){
      for(int k=0;k<subcoords[i].size();k++){
	for(int l=0;l<subcoords[j].size();l++){
	  double dist= subcoords[i][k].eDist(subcoords[j][l]);
          distSetsSecs[sec-1].push_back(dist);
          if(k==(subcoords[i].size()-1)&&l==0 && j==i+1){
	    dist=10.0;
          }
	  if(dist<wRad){
	    triggered = true;
	    overlappedSecs.push_back(dist);
	  }
	}
      }
    }
    //if(triggered==true){
    // overlappedSecs.push_back(i);
    //}
  }
  return overlappedSecs;
}


std::vector<double> ktlMolecule::checkOverlapWithRad(double &wRad){
  std::vector<double> overlappedSecs;
  distSets.clear();
  minimumintraMolecularDistances.clear();
  minimumintraMolecularDistances.resize(chainList.size(),std::vector<double>(chainList.size()));
  for(int i=0;i<minimumintraMolecularDistances.size();i++){
    for(int j=0;j<minimumintraMolecularDistances.size();j++){
      minimumintraMolecularDistances[i][j]=10000.0;
    }
  }
  minimumintraMolecularDistanceMean  =0;
  minimumintraMoleculaeDistancePerChain.clear();
  int chainIndex1=0;
  int chainIndex2=0;
  for(int i=0;i<coords.size();i++){
    bool triggered =false;
    if(i>chainList[chainIndex1].second){
      // if this is an nmer we check which section we are in so as to monitor intra-moelcular distances
      chainIndex1++;
    }
    chainIndex2=chainIndex1;
    for(int j=i;j<coords.size();j++){
      if(j>chainList[chainIndex2].second){
      // if this is an nmer we check which section we are in so as to monitor intra-moelcular distances
      chainIndex2++;
      }
      for(int k=0;k<coords[i].size();k++){
	for(int l=0;l<coords[j].size();l++){
	  double dist= coords[i][k].eDist(coords[j][l]);
	  if(i==j){
	    if(l>k){
	      distSets.push_back(dist);
	    }
	  }else{
	    distSets.push_back(dist);
	  }
	  if(chainIndex1!=chainIndex2){
	    if(dist<minimumintraMolecularDistances[chainIndex1][chainIndex2]){
	      minimumintraMolecularDistances[chainIndex1][chainIndex2]=dist;
	      minimumintraMolecularDistances[chainIndex2][chainIndex1]=dist;
	    }
	  }
          if((k==(coords[i].size()-1)&&l==0 && j==i+1)||j==i){
	    dist=10.0;
          }
	  if(dist<wRad){
	    triggered = true;
	    overlappedSecs.push_back(dist);
	  }
	}
      }
    }
    //if(triggered==true){
    // overlappedSecs.push_back(i);
    //}
  }
  // add last minimum diatnce
  for(int i=0;i<minimumintraMolecularDistances.size();i++){
    /*for(int k=0;k<minimumintraMolecularDistances[i].size();k++){
      std::cout<<" check "<<minimumintraMolecularDistances[i][k]<<"\n";
      }*/
    auto ptrMinElement = std::min_element(minimumintraMolecularDistances[i].begin(),minimumintraMolecularDistances[i].end());
    auto min = *ptrMinElement;
    if(min >  minimumintraMolecularDistanceMean && min<9999.9){
      minimumintraMolecularDistanceMean = min;
    }
    //std::cout<<"min dist for chain "<<i<<" is "<<min<<"\n";
    minimumintraMoleculaeDistancePerChain.push_back(min);
  }
  return overlappedSecs;
}

double ktlMolecule::getMininumIntraMoelcularDistance(){
  return minimumintraMolecularDistanceMean;
}

std::vector<double> ktlMolecule::getMinimumintraMoleculaeDistancePerChain(){
  return minimumintraMoleculaeDistancePerChain;
}

std::vector<std::vector<double> >  ktlMolecule::getMinimumintraMolecularDistances(){
  return minimumintraMolecularDistances;
}

std::vector<double> ktlMolecule::getDistSet(){
  return distSets;
};

double ktlMolecule::compareDistances(std::vector<std::vector<point> > &coords2){
  std::vector<int> overlappedSecs;
  double distDifTot =0.0;
  double n=0;;
  for(int i=0;i<coords.size()-1;i++){
    for(int j=i+1;j<coords.size();j++){
      for(int k=0;k<coords[i].size();k++){
	for(int l=0;l<coords[j].size();l++){
	  double dist= coords[i][k].eDist(coords[j][l]);
	  double dist2= coords2[i][k].eDist(coords2[j][l]);
	  distDifTot = distDifTot+ std::abs(dist2-dist);
	  n=n+1.0;
	}
      }
    }
  }
  return distDifTot/n;
}




bool ktlMolecule::checkCalphas(std::vector<std::vector<point> > &coordsIn){
  bool violated=false;
  for(int i=0;i<coordsIn.size();i++){
    for(int j=0;j<coordsIn[i].size();j++){
      double dist =0.0;
      if((j==coordsIn[i].size()-1) && (i<(coordsIn.size()-1))){
	dist = coordsIn[i][j].eDist(coordsIn[i+1][0]);
      }else if((j==coordsIn[i].size()-1) && (i==coordsIn.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coordsIn[i][j].eDist(coordsIn[i][j+1]);
      }
      std::cout<<i<<" "<<j<<" violated "<<dist<<"\n";
      if((dist>4.0)|| (dist<3.3)){
	violated =true;
	std::cout<<i<<" "<<j<<" violated "<<dist<<"\n";
      }
    }
  }
  return violated;
}

bool ktlMolecule::checkCalphas(){
  bool violated=false;
  for(int i=0;i<coords.size();i++){
    for(int j=0;j<coords[i].size();j++){
      double dist =0.0;
      if((j==coords[i].size()-1) && (i<(coords.size()-1))){
	dist = coords[i][j].eDist(coords[i+1][0]);
      }else if((j==coords[i].size()-1) && (i==coords.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coords[i][j].eDist(coords[i][j+1]);
      }
      if((dist>4.0)|| (dist<3.3)){
	violated =true;
      }
    }
  }
  return violated;
}


bool ktlMolecule::checkCalphas(int &secNo){
  int sec = secNo-1;
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > coordsSec(firstc,secondc);
  bool violated=false;
  for(int i=0;i<coordsSec.size();i++){
    //std::cout<<coordsSec[i].size()<<"\n";
    for(int j=0;j<coordsSec[i].size();j++){
      double dist =0.0;
      if((j==coordsSec[i].size()-1) && (i<(coordsSec.size()-1))){
	dist = coordsSec[i][j].eDist(coordsSec[i+1][0]);
      }else if((j==coordsSec[i].size()-1) && (i==coordsSec.size()-1)){
	dist = 3.7;
      }
      else{
	dist = coordsSec[i][j].eDist(coordsSec[i][j+1]);
      }
      if((dist>4.6)|| (dist<2.6)){
	//std::cout<<i<<" "<<j<<" "<<dist<<"\n";
	violated =true;
      }
    }
  }
  return violated;
}

bool ktlMolecule::checkCalphas(int &secNo,ktlMolecule &ktl){
  int sec = secNo-1;
  std::vector<std::vector<point> > coordsOld = ktl.getCoordinates();
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> >::const_iterator firstcOld =coordsOld.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondcOld =coordsOld.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > coordsSec(firstc,secondc);
  std::vector<std::vector<point> > coordsSecOld(firstcOld,secondcOld);
  bool violated=false;
  for(int i=0;i<coordsSec.size();i++){
    //std::cout<<coordsSec[i].size()<<"\n";
    for(int j=0;j<coordsSec[i].size();j++){
      double dist =0.0;double distOld=0.0;
      if((j==coordsSec[i].size()-1) && (i<(coordsSec.size()-1))){
	dist = coordsSec[i][j].eDist(coordsSec[i+1][0]);
	distOld = coordsSecOld[i][j].eDist(coordsSecOld[i+1][0]);
      }else if((j==coordsSec[i].size()-1) && (i==coordsSec.size()-1)){
	dist = 3.7;
	distOld = 3.7;
      }
      else{
	dist = coordsSec[i][j].eDist(coordsSec[i][j+1]);
	distOld = coordsSecOld[i][j].eDist(coordsSecOld[i][j+1]);
      }
      if(std::abs(dist-distOld)>0.9){
	violated =true;
      }
    }
  }
  return violated;
}


void ktlMolecule::changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList){
  point sp(0.0,0.0,0.0);
  bool suc=true;
  maxDistChange =rmg.reshapeMol(nameSizeSubList,cdsIn,index,sp,suc);
}

int ktlMolecule::getRandomMolecule(){
  coords.clear();
  int noOverLaps=0;
  for(int i=0;i<chainList.size();i++){
     std::cout<<"chain "<<i<<"\n";
    std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[i].first;
    std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[i].second+1;
    std::vector<std::pair<std::string,int> > subns(first,second);
    point sp(0.0,0.0,0.0);
    bool suc=true;
    std::vector<std::vector<point> > submol = rmg.makeRandomMolecule(subns,sp,suc);
    std::cout<<"made man\n?";
    std::vector<int> overLapList = checkOverlap(submol);
    std::vector<int> currOverLapList= overLapList;
    std::random_device rdev{};
    std::default_random_engine generator1{rdev()};
    int l=0;
    while(l<10000 && currOverLapList.size()>0){
      std::vector<std::vector<point> > tempCoords = submol;
      //std::cout<<currOverLapList.size()<<" "<<tempCoords.size()<<" "<<l<<"\n";
      if(currOverLapList.size()>l){
	std::uniform_int_distribution<int> distributionR(0,currOverLapList.size()-1);
	int iv= distributionR(generator1);
	if(currOverLapList[iv]==0){
	  changeMoleculeSingle(currOverLapList[1],tempCoords,subns);
	}else{
	  changeMoleculeSingle(currOverLapList[iv],tempCoords,subns);
	}
      }else{
	changeMoleculeSingle(currOverLapList[0],tempCoords,subns);
      }
      overLapList = checkOverlap(tempCoords);
      bool caca = checkCalphas(tempCoords);
      if(overLapList.size()<currOverLapList.size() && caca==false){
	currOverLapList=overLapList;
	submol = tempCoords;
      }
      l++;
    }
    noOverLaps = noOverLaps + currOverLapList.size();
    // now add the section to the whole
    for(int j=0;j<submol.size();j++){
      coords.push_back(submol[j]);
    }
  }
  return noOverLaps;
}





void ktlMolecule::changeMoleculeSingleMulti(int &index,int secIn){
  int sec = secIn-1;
  point sp(0.0,0.0,0.0);
  std::vector<std::pair<std::string,int> >::const_iterator first =nameSizeList.begin()+chainList[sec].first;
  std::vector<std::pair<std::string,int> >::const_iterator second =nameSizeList.begin()+chainList[sec].second+1;
  std::vector<std::pair<std::string,int> > subns(first,second);
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  bool suc=true;
  //calculate the size of the section to check after if we have lost any molecules
  int preSize = 0;
  for(int i=0;i<subcoords.size();i++){
    preSize = preSize+subcoords[i].size();
  }
  if(index<subcoords.size()-2){
    // // check for repeat heliced
    if(subns[index].first=="loop"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
      //change loo then encounter repeat helices
       //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<subcoords.size()&&subns[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolLoopThenHelixSet(subns,subcoords,index,noHelices,sp,suc);
    }else if(subns[index].first=="Helix"&&subns[index+1].first=="Helix"&&subns[index+2].first=="Helix"){
      //change helix and encounter repreat helices
      int noHelices=0;int k=index+1;
      while(k<subcoords.size()&&subns[k].first=="Helix"){
	noHelices++;
	k++;
      }
      maxDistChange =rmg.reshapeMolHelixSet(subns,subcoords,index,noHelices,sp,suc);
    }else{
      maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
    }
  }else{
    maxDistChange =rmg.reshapeMol(subns,subcoords,index,sp,suc);
  }
  // now replace the new sections into the main chain
  //check the size
  int postSize = 0;
  for(int i=0;i<subcoords.size();i++){
    postSize = postSize+subcoords[i].size();
  }
  //std::cout<<preSize<<" "<<postSize<<"\n";
  if(preSize==postSize){
    int d = chainList[sec].second+1-chainList[sec].first;
    std::copy_n(subcoords.begin(),d,&coords[chainList[sec].first]);
  }
}



void ktlMolecule::changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec){
  point sp(0.0,0.0,0.0);
  int sec = secIn-1;
  std::vector<std::vector<point> >::const_iterator firstc =coords.begin()+chainList[sec].first;
  //rotate from here until end of molecule
  std::vector<std::vector<point> >::const_iterator secondc =coords.begin()+chainList[sec].second+1;
  std::vector<std::vector<point> > subcoords(firstc,secondc);
  bool suc=true;
  point com = getCentreOfMass(subcoords);
  rotateSection(subcoords,com,k,angle,transVec);
  // now replace the new sections into the main chain
  int di = chainList[sec].second+1-chainList[sec].first;
  std::copy_n(subcoords.begin(),di,&coords[chainList[sec].first]);
}

void ktlMolecule::replicateMolecule(int &noReplications){
  std::random_device rdev{};
  std::default_random_engine generator1{rdev()};
  std::uniform_real_distribution<> rotAng(0.0,0.5);
  std::uniform_real_distribution<> theAng(0.0,3.14159265359);
  std::uniform_real_distribution<> phiAng(0.0,6.28318530718);
  point com = getCentreOfMass(coords);
  double maxRad=0.0;
  for(int i=0;i<coords.size();i++){
    for(int j=0;j<coords[i].size();j++){
      double rad = com.eDist(coords[i][j]);
      if(rad>maxRad){
        maxRad = rad;
      }
    }
  }
  std::vector<std::vector<point> > replicatedCoords = coords;
  std::vector<std::pair<std::string,int> > replicatedNameSizeList=nameSizeList;
  std::vector<std::pair<int,int> > chainListCopy = chainList;
  std::vector<std::vector<std::string> > aminoListCopy = aminoList;
  for(int i=1;i<=noReplications;i++){
    std::vector<std::vector<point> > newPts = coords;
    double angle = rotAng(generator1);double theta=theAng(generator1);double phi = phiAng(generator1);
    point k(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));
    theta=theAng(generator1);phi = phiAng(generator1);
    point translate(std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta));
    translate.scalarMult(2.2*maxRad);
    rotateSection(newPts,com,k,angle,translate);
    // now append
    replicatedCoords.insert(replicatedCoords.end(),newPts.begin(),newPts.end());
    // now add on to the nameList separator
    std::vector<std::pair<std::string,int> > copyNameList;
    for(int j=0;j<nameSizeList.size();j++){
      std::pair<std::string,int> pr;
      pr.first=nameSizeList[j].first;
      pr.second=nameSizeList[j].second;
      copyNameList.push_back(pr);
    }
    replicatedNameSizeList.insert(replicatedNameSizeList.end(),copyNameList.begin(),copyNameList.end());
    std::pair<int,int> prInts;
    prInts.first=chainListCopy[chainListCopy.size()-1].second+1;
    prInts.second=chainListCopy[chainListCopy.size()-1].second+chainList[chainList.size()-1].second+1;
    chainListCopy.push_back(prInts);
    //finally replicate the aminoList and hydration list
    aminoListCopy.insert(aminoListCopy.begin(),aminoList.begin(),aminoList.end());
  }
  aminoList = aminoListCopy;
  nameSizeList=replicatedNameSizeList;
  chainList=chainListCopy;
  coords=replicatedCoords;
}



void ktlMolecule::rotation3D(point &p,point &centre,point &k,double &cosangle,double &sinangle){
  point out;
  // translate so the centre is at the origin
  out = p-centre;
  point term1 = out*cosangle;
  point term2 = k.cross(out);term2 = term2*sinangle;
  double dotProd = k.dotprod(out);
  point term3 = k*dotProd;
  term3 = term3*(1.0-cosangle);
  out = term1+term2 +term3;
  out = out +centre;
  p=out;
}

void ktlMolecule::rotateSection(std::vector<std::vector<point> >  &section,point &centre,point &k,double &angle,point &transVec){
  double cosangle = std::cos(angle);double sinangle = std::sin(angle);
  for(int i=0;i<section.size();i++){
    for(int j=0;j<section[i].size();j++){
      rotation3D(section[i][j],centre,k,cosangle,sinangle);
      section[i][j] = section[i][j]+transVec;
    }
  }
}

void ktlMolecule::writeMoleculeToFile(const char* filename){
  std::ofstream ofile;
  ofile.open(filename);
  //std::cout<<filename<<"\n";
  if(ofile.is_open()){
     for(int k=0;k<chainList.size();k++){
      for(int i=chainList[k].first;i<=chainList[k].second;i++){
       for(int j=0;j<coords[i].size();j++){
	 //std::cout<<coords[i][j].getX()<<" "<<coords[i][j].getY()<<" "<<coords[i][j].getZ()<<"\n";
	   ofile<<coords[i][j].getX()<<" "<<coords[i][j].getY()<<" "<<coords[i][j].getZ()<<"\n";
       }
       ofile<<"\n";
      }
      ofile<<"End chain "<<k+1<<"\n";
    }
  }else{
    ofile<<"cannot write molecule to file";
  }
  ofile.close();
}

std::vector<std::pair<double,double> > ktlMolecule::getKapTauVals(){
  std::vector<std::pair<double,double> > ktllst;
  for(int k=0;k<chainList.size();k++){
    for(int i=chainList[k].first;i<=chainList[k].second;i++){
      for(int j=0;j<coords[i].size();j++){
         if(i!=chainList[k].second){
          if(j<coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i][j+3]));
          }else if(j==coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i+1][0]));
          }else if(j==coords[i].size()-2){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i+1][0],coords[i+1][1]));
          }else if(j==coords[i].size()-1){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][0],coords[i+1][1],coords[i+1][2]));
          }
        }else{
          if(j<coords[i].size()-3){
            ktllst.push_back(rmg.kapTau(coords[i][j],coords[i][j+1],coords[i][j+2],coords[i][j+3]));
          }
        }
      }
    }
  }
  return ktllst;
}



double ktlMolecule::lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff){
  //for actual Leennard jones
  /*double rat = idealDist/currDist;
  double edist = 1.0;
  double sixpow =rat*rat*rat*rat*rat*rat;
  return (edist*(sixpow*sixpow-2.0*sixpow)+edist)/double(noConPred);
  */
  double dif = (idealDist-currDist)/idealDist;
  return weightCoeff*dif*dif;
}

std::vector<double> ktlMolecule::solMolDists(std::vector<std::vector<point> > &pts1){
  std::vector<double> distSet;
  for(int i=0;i<pts1.size();i++){
    for(int j=0;j<pts1[i].size();j++){
      for(int k=0;k<coords.size();k++){
        for(int l=0;l<coords[k].size();l++){
          double d =pts1[i][j].eDist(coords[k][l]);
          distSet.push_back(d);
        }
      }
    }
  }
  return distSet;
}


void ktlMolecule::loadContactPredictions(const char* contactloc){
  std::ifstream cpfile;
  cpfile.open(contactloc);
  if(cpfile.is_open()){
     std::string output;
     while(!cpfile.eof()){
       int ind1,ind2,ind3,ind4;double distance,percentage;
       std::getline(cpfile,output);
       std::stringstream ss(output);
       if(output.length()>0){
       ss>>ind1;
       ss.ignore();
       ss>>ind2;
       std::pair<int,int> pr1;
       pr1.first=ind1;pr1.second=ind2;
       ss.ignore();
       ss>>ind3;
       ss.ignore();
       ss>>ind4;
       ss.ignore();
       std::pair<int,int> pr2;
       pr2.first=ind3;pr2.second=ind4;
       ss>>distance;
       ss.ignore();
       ss>>percentage;
       std::pair<double,double> pr3;
       pr3.first=distance;pr3.second=percentage;
       std::tuple<std::pair<int,int>,std::pair<int,int>,std::pair<double,double> > tp;
       std::get<0>(tp) = pr1;std::get<1>(tp) = pr2;std::get<2>(tp) = pr3;
       contactPairList.push_back(tp);
       }
     }
     int cpls =contactPairList.size();
  }
}

double ktlMolecule::getLennardJonesContact(){
  double ljval =0.0;
  for(int i=0;i<contactPairList.size();i++){
    std::tuple<std::pair<int,int>,std::pair<int,int>,std::pair<double,double>> tp = contactPairList[i];
    std::pair<int,int> pr1 =  std::get<0>(tp);
    std::pair<int,int> pr2 =  std::get<1>(tp);
    point cd1 = coords[pr1.first][pr1.second];
    point cd2 = coords[pr2.first][pr2.second];
    double prWiseDist = cd1.eDist(cd2);
    //std::cout<<"pair "<<i<<"\n";
    //std::cout<<pr1.first<<" "<<pr1.second<<" "<<pr2.first<<" "<<pr2.second<<"\n";
    //cd1.printPoint();
    //cd2.printPoint();
    //std::cout<<"size1 "<<coords[pr1.first].size()<<"\n";
    //std::cout<<"size2 "<<coords[pr2.first].size()<<"\n";
    //std::cout<<"d comp "<<prWiseDist<<" "<<std::get<2>(tp).first<<"\n";
    double distFrac=(prWiseDist-std::get<2>(tp).first)/(std::get<2>(tp).first);
    //std::cout<<distFrac<<"\n";
    double distFracWeighted = distFrac/(std::get<2>(tp).second);
    //std::cout<<distFracWeighted<<" "<<std::get<2>(tp).second<<"\n";
    //std::cout<<"\n";
    ljval = ljval + 0.0001*distFracWeighted*distFracWeighted*distFracWeighted*distFracWeighted;
    //std::cout<<"pen bounded ?"<<distFracWeighted<<"\n";
    //cd1.printPoint();
    //cd2.printPoint();
    //double drat = prWiseDist/std::get<2>(tp);
    //double indrat  = 1.0/(std::sqrt(2.0)*drat);
    //ljval =ljval+ (4.0*(indrat*indrat*indrat*indrat -indrat*indrat)+1.0);
  }
  /*if(contactPairList.size()>0){
    return ljval;
  }else{
    return ljval;
    }*/
  return ljval;
}
void ktlMolecule::loadFixedSections(const char* fixedsecloc){
  std::ifstream fsfile;
  fsfile.open(fixedsecloc);
  if(fsfile.is_open()){
    std::string output;
     while(!fsfile.eof()){
       std::getline(fsfile,output);
       std::stringstream ss(output);
       int section;
       ss>>section;
       unchangedSections.push_back(section);
     }
  }
}


double ktlMolecule::getBetaSheetProximityReward() {

    std::vector<std::vector<point>> betaSheets;

    // Find all betasheets, add section coords to vec
    for (int sec=0; sec<coords.size(); sec++) {

        if (nameSizeList[sec].first == "Strand") {
            betaSheets.push_back(coords[sec]);
        }
    }

    numBetaSheets = betaSheets.size();

    // If your protein has fewer than 2 sheets, no reward buddy!
    if (betaSheets.size() < 2) {
        return 0.0;
    }

    double totalReward = 0.0;

    // Go through and find min distance between any of the sheets
    // maybe I want a count of number of dist under 5A or w/e
    for (int i=0; i<betaSheets.size()-1; i++) {
        for (int j=i+1; j<betaSheets.size(); j++) {

            double minDistance = 10000; // rand large number

            for (auto& point1 : betaSheets[i]) {
                for (auto& point2 : betaSheets[j]){
                    double distance = point1.eDist(point2);
                    if (distance < minDistance) {
                        minDistance = distance;
                    }
                }
            }

            // check is sheet is (anti) parallel
            bool isAligned = areSheetsPairallelOrAntiparallel(betaSheets[i], betaSheets[j]);

            // some reward metric based on closeness and aligned-ness
            if(isAligned) {

                            double reward = 0.0;

                            // Reward close, aligned sheets. Adjust these parameters as needed.
                            double Dopt = 5.0;  // Optimal distance in Angstroms
                            double Dmax = 10.0;  // Max distance for any reward
                            double Dmin = 3.8; // Min distance - dont want to reward for getting too close!

                            if(minDistance <= Dmin) {

                                reward = 0.0;

                            } else if (minDistance <= Dmax) {

                                reward = 1.0 - std::abs(minDistance - Dopt) / (Dmax - Dmin);

                            }

                            totalReward += reward;

                            }
        }
    }

    return totalReward;

}


bool ktlMolecule::areSheetsPairallelOrAntiparallel(std::vector<point>& sheet1, std::vector<point>& sheet2) {
    // Calculate direction vectors for both sheets
    point dir1 = sheet1.back() - sheet1.front();
    point dir2 = sheet2.back() - sheet2.front();

    // Normalize the vectors
    dir1.normalise();
    dir2.normalise();

    // Calculate (unsigned) angle between the two sheets
    double dotProduct = std::abs(dir1.dotprod(dir2));

    // If dot product is close to 1 or -1, sheets are parallel or antiparallel
    return dotProduct > 0.8;  // Adjust this threshold as needed
}
