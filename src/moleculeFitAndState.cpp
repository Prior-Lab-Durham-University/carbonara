#include "moleculeFitAndState.h"
#include "ktlMoleculeRandom.h"


moleculeFitAndState::moleculeFitAndState(std::vector<ktlMolecule> &molin, ModelParameters& params){
  //define the number of structures
  mol = molin;
  molDists.resize(mol.size());
  molSize.resize(mol.size());
  maxDistMol.resize(mol.size());
  contactPredPen.resize(mol.size());
  writhePenalty=0.0;
  connectionPenaltySet.resize(mol.size());
  
  // set the fixed fitting parameters

  closestApproachDist = params.closestApproachDist;
  rmin = params.rmin; rmax = params.rmax; 
  lmin = params.lmin;


  for(int i=0;i<mol.size();i++){
   writheFP wfp;
   std::vector<double> chainWrithes;
   for(int j=0;j<mol[i].noChains();j++){
     std::vector<std::vector<point> > crds =mol[i].getSubsecCoordinates(j);
     double wr =wfp.DIDownSampleAbsSingle(crds);
     chainWrithes.push_back(wr);
     //std::cout<<"initial abs writhe moelcule "<<i<<" chain "<<j<<" "<<wr<<"\n";
   }
   originalWrithes.push_back(chainWrithes);
  }
  currWrithes =originalWrithes;
  connectionPenalty=0.0;

  for(int i=0;i<mol.size();i++){
    //calculate distances
    //std::cout<<"for molecule "<<i<<"\n";
    std::vector<double> overlapDists= mol[i].checkOverlapWithRad(closestApproachDist);
    molDists[i] = mol[i].getDistSet();
    //get number of amino acids
    molSize[i]  = mol[i].getNoAminos();
    //sort the distances from largest to smallest for binning.
    std::sort(molDists[i].begin(),molDists[i].end());
    maxDistMol[i] = molDists[i][molDists[i].size()-1];
    // calculate any overlap distances
    overlapDistSet.push_back(overlapDists);
    double meanDist = mol[i].getMininumIntraMoelcularDistance();
     std::vector<double> minIntraChainDistances = mol[i].getMinimumintraMoleculaeDistancePerChain();
    for(int j=0;j<minIntraChainDistances.size();j++){
      if(minIntraChainDistances[j]<9999.9 && minIntraChainDistances[j]>5.0){
	// only penalise nmers not monomer
	  double dif = minIntraChainDistances[j] - 5.0;
	  connectionPenaltySet[i] = connectionPenaltySet[i] +0.0005*dif*dif*dif*dif*dif*dif*dif*dif*dif*dif;
      }
    }
    connectionPenaltySet[i]=connectionPenaltySet[i]/double(mol[i].noChains());
    connectionPenalty = connectionPenalty + connectionPenaltySet[i];
  }

  // to fill on first calc
  originalOverlapPenalty= 0.0;
  double overlapPenalty = applyOverlapPenalty();
  originalOverlapPenalty= overlapPenalty;
}

std::vector<ktlMolecule> moleculeFitAndState::getMolecule(){
  return mol;
}

void  moleculeFitAndState::updateMolecule(std::vector<ktlMolecule> &molNew){
  mol=molNew;
}

// version where a single changed section has been made and we update the distances

void moleculeFitAndState::calculateMoleculeDistances(ktlMolecule &molNew,int &i){
  std::vector<double> overlapDists= molNew.checkOverlapWithRad(closestApproachDist);
  molDists[i] = molNew.getDistSet();
  //get number of amino acids
  molSize[i]  = molNew.getNoAminos();
  //sort the distances from largest to smallest for binning.
  std::sort(molDists[i].begin(),molDists[i].end());
  maxDistMol[i] = molDists[i][molDists[i].size()-1];
  // calculate any overlap distances
  overlapDistSet[i]=overlapDists;
}



void  moleculeFitAndState::writeScatteringToFile(experimentalData &ed,std::vector<std::vector<double> > &mixtureVals,const char* filename){
  ed.writeScatteringToFile(mixtureVals,filename);
}

void  moleculeFitAndState::writeScatteringToFile_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureVals,const char* filename){
  ed.writeScatteringToFile_ChiSq(mixtureVals,filename);
}



double moleculeFitAndState::getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists){
  double distSumCurr=0.0;
  for(int l=0;l<overlapDists.size();l++){
    // Uncheck me to see what overalps we get...
    // std::cout<<l<<" "<<overlapDists[l]<<"\n";
    double dist = closestApproachDist-overlapDists[l];
    distSumCurr = distSumCurr + std::exp(dist)-1.0;
  }
  if(overlapDists.size()>0){
    distSumCurr =0.01*distSumCurr/overlapDists.size();
  }
  //std::cout<<"Distance penalty "<<distSumCurr<<"\n";
  return distSumCurr;
}

double moleculeFitAndState::applyOverlapPenalty(){
   double overlapPenalty = 0.0;
   for(int i=0;i<overlapDistSet.size();i++){
     // calculate any overlap distances
     overlapPenalty = overlapPenalty + getOverlapPenalty(closestApproachDist,overlapDistSet[i]);
  }
   //std::cout<<overlapPenalty<<" "<<originalOverlapPenalty<<"\n";
   if(overlapPenalty > originalOverlapPenalty){
     return overlapPenalty-originalOverlapPenalty;
   }else{
     return 0.0;
   }
}

double moleculeFitAndState::applyDistanceConstraints(){
 double contactPredPenTotal=0.0;
  for(int i=0;i<mol.size();i++){
    contactPredPen[i] = mol[i].getLennardJonesContact();
    contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
  }
  return contactPredPenTotal;
}

double moleculeFitAndState::applyDistanceConstraints(ktlMolecule &molNew,int &im){
 double contactPredPenTotal=0.0;
  for(int i=0;i<contactPredPen.size();i++){
    if(i==im){
      contactPredPen[i] = molNew.getLennardJonesContact();
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }else{
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }
  }
  return contactPredPenTotal;
}


void moleculeFitAndState::alterWritheSet(ktlMolecule &molNew,int &i){
   writheFP wfp;
    for(int j=0;j<molNew.noChains();j++){
     std::vector<std::vector<point> > crds =molNew.getSubsecCoordinates(j);
     double wr =wfp.DIDownSampleAbsSingle(crds);
     currWrithes[i][j] = wr;
    }
}

// calculate writhe lists
void moleculeFitAndState::applyWritheConstraint(){
  writhePenalty=0.0;
  for(int i=0;i<currWrithes.size();i++){
    for(int j=0;j<currWrithes[i].size();j++){
      double newWrithe =currWrithes[i][j];
      int index = j+1;
      double secLen = double(mol[i].getSubsecSize(index));
      double lowerBound = std::pow((secLen/7.5),1.6)-3.0;
      writhePenalty=  writhePenalty+1.0/(1.0+std::exp(20.0*(newWrithe-lowerBound)));
    }
  }
}

// the following funtion is for when we want to create nmers and keep them "connected"

void moleculeFitAndState::calculateConnectionPenalty(ktlMolecule &molNew,int &chInd){
  connectionPenaltySet[chInd]=0.0;
  std::vector<double> minIntraChainDistances= molNew.getMinimumintraMoleculaeDistancePerChain();
    //std::cout<<meanDist<<"\n";
   for(int j=0;j<minIntraChainDistances.size();j++){
     //std::cout<<" connie dif "<<minIntraChainDistances[j]<<"\n";
      if(minIntraChainDistances[j]<9999.9 && minIntraChainDistances[j]>5.0){
	// only penalise nmers not monomer
	  double dif =minIntraChainDistances[j] - 5.0;
	  connectionPenaltySet[chInd] = connectionPenaltySet[chInd] +0.0005*dif*dif*dif*dif*dif*dif*dif*dif*dif*dif;
      }
    }
   connectionPenaltySet[chInd] = connectionPenaltySet[chInd]/double(mol[chInd].noChains());
   connectionPenalty = 0.0;
   for(int i=0;i<connectionPenaltySet.size();i++){
     connectionPenalty = connectionPenalty + connectionPenaltySet[i];
   }
}



double moleculeFitAndState::calculateUserSpecifiedConnectionPenalty(ktlMolecule &molNew,std::vector<int> &chainSet1,std::vector<int> &chainSet2){
  std::vector<std::vector<double> > minIntraChainDistances= molNew.getMinimumintraMolecularDistances();
  double userConnectedPenalty;
  double minDist = 10000.0;
  //std::cout<<"here ? "<<minIntraChainDistances.size()<<"\n";
  for(int i=0;i<chainSet1.size();i++){
    for(int j=0;j<chainSet2.size();j++){
      if(minIntraChainDistances[chainSet1[i]][chainSet2[j]]<minDist){
	minDist = minIntraChainDistances[chainSet1[i]][chainSet2[j]];
      }
    }
  }
  if(minDist<9999.0 && minDist>5.0){
     double dif =minDist - 5.0;
     userConnectedPenalty = userConnectedPenalty +0.0005*dif*dif*dif*dif*dif*dif*dif*dif*dif*dif;
  }else{
    userConnectedPenalty = 0.0;
  }
  return userConnectedPenalty;
}


double moleculeFitAndState::calculateUserSpecifiedConnectionPenalty(int chInd,std::vector<int> &chainSet1,std::vector<int> &chainSet2){
  std::vector<std::vector<double> > minIntraChainDistances= mol[chInd].getMinimumintraMolecularDistances();
  double userConnectedPenalty;
  double minDist = 10000.0;
  for(int i=0;i<chainSet1.size();i++){
    for(int j=0;j<chainSet2.size();j++){
      if(minIntraChainDistances[chainSet1[i]][chainSet2[j]]<minDist){
	minDist = minIntraChainDistances[chainSet1[i]][chainSet2[j]];
      }
    }
  }
  if(minDist<9999.0 && minDist>5.0){
     double dif =minDist - 5.0;
     userConnectedPenalty = userConnectedPenalty +0.0005*dif*dif*dif*dif*dif*dif*dif*dif*dif*dif;
  }else{
    userConnectedPenalty = 0.0;
  }
  return userConnectedPenalty;
}



double moleculeFitAndState::getWrithePenalty(){
  return writhePenalty;
}

double moleculeFitAndState::getOverlapPenalty(){
  return applyOverlapPenalty();
}

double  moleculeFitAndState::getDistanceConstraints(){
  return applyDistanceConstraints();
}



std::pair<double,double> moleculeFitAndState::getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  // get the scattering
  double scatterAndHydrationConstraint = ed.calculateChiSquared(mol,kmin,kmax,mixtureList);
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.

  **************************************************************/
  //std::cout<<"scattering "<<scatterAndHydrationConstraint<<"\n";
  double overlapPenalty = applyOverlapPenalty();
  //std::cout<<"Overlap Constraints "<<overlapPenalty<<"\n";
  double distanceConstraints = applyDistanceConstraints();
  //std::cout<<"Distance Constraints "<<distanceConstraints<<"\n";
  applyWritheConstraint();
  //std::cout<<"Writhe penalty "<<writhePenalty<<"\n";
  //calculateConnectionPenalty(mol[0],0);
  double currFit = scatterAndHydrationConstraint +1.0*(distanceConstraints + writhePenalty+overlapPenalty);
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  // pass along best hydration C2 parameter
  return fitStats;
}

std::pair<double,double> moleculeFitAndState::getOverallFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  // get the scattering   
  double scatterAndHydrationConstraint = ed.calculateChiSquared_Weighted(mol,kmin,kmax,mixtureList);
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.

  **************************************************************/
  
  double overlapPenalty = applyOverlapPenalty();
  //std::cout<<"Overlap penalty "<<overlapPenalty<<"\n";
  double distanceConstraints = applyDistanceConstraints();
  //std::cout<<"Distance Constraints "<<distanceConstraints<<"\n";
  applyWritheConstraint();
  //std::cout<<"Writhe penalty "<<writhePenalty<<"\n";
  //calculateConnectionPenalty(mol[0],0);
  double currFit = scatterAndHydrationConstraint +100.0*(distanceConstraints + writhePenalty+overlapPenalty);
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  // pass along best hydration C2 parameter
  return fitStats;
}


std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  // get the scattering   
  double scatterAndHydrationConstraint = ed.calculateChiSquared(mol,kmin,kmax,mixtureList);;
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.

  **************************************************************/
  double overlapPenalty = applyOverlapPenalty();
  double distanceConstraints = applyDistanceConstraints();
  //std::cout<<"Distance Constraints "<<distanceConstraints<<"\n";
  applyWritheConstraint();
  //std::cout<<"Writhe penalty "<<writhePenalty<<"\n";
  //std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  for(int i=0;i<mol.size();i++){
    calculateConnectionPenalty(mol[i],i);
  }
  std::cout<<"original connection Pen "<<connectionPenalty<<"\n";
  // if the user has specidfed some sub chains to be connected.
  double currFit = scatterAndHydrationConstraint  +1.0*(distanceConstraints + writhePenalty+connectionPenalty);
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  return fitStats;
}

std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  // get the scattering   
  double scatterAndHydrationConstraint = ed.calculateChiSquared_Weighted(mol,kmin,kmax,mixtureList);;
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.

  **************************************************************/
  double overlapPenalty = applyOverlapPenalty();
  double distanceConstraints = applyDistanceConstraints();
  //std::cout<<"Distance Constraints "<<distanceConstraints<<"\n";
  applyWritheConstraint();
  //std::cout<<"Writhe penalty "<<writhePenalty<<"\n";
  //std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  for(int i=0;i<mol.size();i++){
    calculateConnectionPenalty(mol[i],i);
  }
  //std::cout<<"original connection Pen "<<connectionPenalty<<"\n";
  // if the user has specidfed some sub chains to be connected.
  double currFit = scatterAndHydrationConstraint  +100.0*(distanceConstraints +overlapPenalty+ writhePenalty+connectionPenalty);
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  return fitStats;
}



std::pair<double,double> moleculeFitAndState::getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i){
  // update the molecule distances for molecule i;
  calculateMoleculeDistances(molNew,i);
  double scatterAndHydrationConstraint = ed.calculateChiSquaredUpdate(molNew,i,kmin,kmax,mixtureList);
  // apply penalties
  //std::cout<<"updated "<<scatterAndHydrationConstraint<<"\n";
   double overlapPenalty = applyOverlapPenalty();
   // std::cout<<"Overlap Penalty update "<<overlapPenalty<<"\n";
   double distanceConstraints = applyDistanceConstraints(molNew,i);
   //std::cout<<"Distance constraints update "<<distanceConstraints<<"\n";
  alterWritheSet(molNew,i);
  applyWritheConstraint();
  // std::cout<<" writhe penalty  "<<writhePenalty<<"\n";
  //calculateConnectionPenalty(molNew,i);
  //std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  //std::cout<<" connection penalty  "<<connectionPenalty<<"\n";
  double currFit = scatterAndHydrationConstraint +1.0*(overlapPenalty +distanceConstraints + writhePenalty);
  //std::cout<<currFit<<"\n";
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;


  return fitStats;
}


std::pair<double,double> moleculeFitAndState::getOverallFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i){
  // update the molecule distances for molecule i;
  calculateMoleculeDistances(molNew,i);
  double scatterAndHydrationConstraint = ed.calculateChiSquaredUpdate_Weighted(molNew,i,kmin,kmax,mixtureList);
  // apply penalties
  //std::cout<<"updated "<<scatterAndHydrationConstraint<<"\n";
   double overlapPenalty = applyOverlapPenalty();
    //std::cout<<"Overlap Penalty "<<overlapPenalty<<"\n";
   double distanceConstraints = applyDistanceConstraints(molNew,i);
   // std::cout<<"Distance constraints "<<distanceConstraints<<"\n";
  alterWritheSet(molNew,i);
  applyWritheConstraint();
  //std::cout<<" writhe penalty  "<<writhePenalty<<"\n";
  //calculateConnectionPenalty(molNew,i);
  //std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  //std::cout<<" connection penalty  "<<connectionPenalty<<"\n";
  double currFit = scatterAndHydrationConstraint +100.0*(overlapPenalty +distanceConstraints + writhePenalty);
  //std::cout<<currFit<<"\n";
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;


  return fitStats;
}


std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i){
  // update the molecule distances for molecule i;
  calculateMoleculeDistances(molNew,i);
  int bestHelRatList=0;
  // now update the hydration shell
  double scatterAndHydrationConstraint =  ed.calculateChiSquaredUpdate(molNew,i,kmin,kmax,mixtureList);
  // apply penalties
  double overlapPenalty = applyOverlapPenalty();
  //std::cout<<"Overlap Penalty "<<overlapPenalty<<"\n";
  double distanceConstraints = applyDistanceConstraints(molNew,i);
  //std::cout<<"Distance constraints "<<distanceConstraints<<"\n";
  alterWritheSet(molNew,i);
  applyWritheConstraint();
  //std::cout<<" writhe penalty  "<<writhePenalty<<"\n";
  calculateConnectionPenalty(molNew,i);
  //std::cout<<" scattering  "<<scatterAndHydrationConstraint<<"\n";
  //std::cout<<" connection penalty  "<<connectionPenalty<<"\n";
  double currFit = scatterAndHydrationConstraint +1.0*(overlapPenalty +distanceConstraints + writhePenalty+connectionPenalty);
  //std::cout<<currFit<<"\n";
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  return fitStats;
}


std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i){
  // update the molecule distances for molecule i;
  calculateMoleculeDistances(molNew,i);
  int bestHelRatList=0;
  // now update the hydration shell
  double scatterAndHydrationConstraint =  ed.calculateChiSquaredUpdate_Weighted(molNew,i,kmin,kmax,mixtureList);
  // apply penalties
  double overlapPenalty = applyOverlapPenalty();
  //std::cout<<"Overlap Penalty change "<<overlapPenalty<<"\n";
  double distanceConstraints = applyDistanceConstraints(molNew,i);
  //std::cout<<"Distance constraints "<<distanceConstraints<<"\n";
  alterWritheSet(molNew,i);
  applyWritheConstraint();
  //std::cout<<" writhe penalty change "<<writhePenalty<<"\n";
  calculateConnectionPenalty(molNew,i);
  //std::cout<<" scattering change "<<scatterAndHydrationConstraint<<"\n";
  // std::cout<<" connection penalty change "<<connectionPenalty<<"\n";
  double currFit = scatterAndHydrationConstraint +100.0*(overlapPenalty +distanceConstraints + writhePenalty+connectionPenalty);
  //std::cout<<currFit<<"\n";
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  return fitStats;
}


/*When improved run this to update predicted scattering will make sure all I_model calculations are up to date with the current molecule*/

void moleculeFitAndState::updateScatteringFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  for(int i=0;i<mol.size();i++){
    calculateMoleculeDistances(mol[i],i);
  }
  double dummyScatter = ed.calculateChiSquared(mol,kmin,kmax,mixtureList);
  std::cout<<"double check baseLine "<<dummyScatter<<"\n";
}

void moleculeFitAndState::updateScatteringFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax){
  for(int i=0;i<mol.size();i++){
    calculateMoleculeDistances(mol[i],i);
  }
  double dummyScatter = ed.calculateChiSquared_Weighted(mol,kmin,kmax,mixtureList);
  //std::cout<<"double check baseLine "<<dummyScatter<<"\n";
}


double moleculeFitAndState::getBetaSheetReward() {

    double sheetRewards = 0.0;

    for(int i=0;i<mol.size();i++){

        double numSheets = mol[i].numBetaSheets;
        sheetRewards += mol[i].getBetaSheetProximityReward()/numSheets;

    }


    // for(const auto& molecule : mol) {
    //         sheetRewards += molecule.getBetaSheetProximityReward();
    //     }


    return sheetRewards;

}
