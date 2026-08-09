#include "moleculeFitAndState.h"
#include "ktlMoleculeRandom.h"


moleculeFitAndState::moleculeFitAndState(std::vector<ktlMolecule> &molin, ModelParameters& params){
  //define the number of structures
  mol = molin;
  molDists.resize(mol.size());
  molSize.resize(mol.size());
  maxDistMol.resize(mol.size());
  contactPredPen.resize(mol.size());
  perPairPenaltiesPerMol.resize(mol.size());
  hardOKPerMol.assign(mol.size(), true);
  hardMaxViolationPerMol.assign(mol.size(), 0.0);
  writhePenalty=0.0;
  connectionPenaltySet.resize(mol.size());
  maxWritheDiff = params.maxWritheDiff;
  writheDiffStride = params.writheDiffStride;
  writheDiffPenalty = 0.0;
  penaltyWeight = params.penaltyWeight;

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

// sums, across mixture states, the per-pair minimum penalty among pairs flagged
// ensemble-OR -- "satisfied by any one conformation" instead of "every conformation must
// satisfy it". Pairs not flagged ensemble-OR (and hard pairs) are 0.0 in every state's
// vector by construction (ktlMolecule::getSoftContactPenaltiesPerPair), so they
// contribute nothing here -- already counted via the regular per-state sum instead.
// Assumes every mixture state's constraint file is the identical replica
// (replicate_numbered_files' normal behaviour): same pair count, same order.
double moleculeFitAndState::ensembleOrPenalty(){
  double total = 0.0;
  if(perPairPenaltiesPerMol.empty() || perPairPenaltiesPerMol[0].empty()){ return total; }
  int numPairs = perPairPenaltiesPerMol[0].size();
  for(int k=0;k<numPairs;k++){
    double minAcrossStates = perPairPenaltiesPerMol[0][k];
    for(int i=1;i<perPairPenaltiesPerMol.size();i++){
      if(k<perPairPenaltiesPerMol[i].size() && perPairPenaltiesPerMol[i][k] < minAcrossStates){
        minAcrossStates = perPairPenaltiesPerMol[i][k];
      }
    }
    total += minAcrossStates;
  }
  return total;
}

double moleculeFitAndState::applyDistanceConstraints(){
 double contactPredPenTotal=0.0;
  for(int i=0;i<mol.size();i++){
    contactPredPen[i] = mol[i].getLennardJonesContact();
    perPairPenaltiesPerMol[i] = mol[i].getSoftContactPenaltiesPerPair();
    contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
  }
  contactPredPenTotal += ensembleOrPenalty();
  return contactPredPenTotal;
}

double moleculeFitAndState::applyDistanceConstraints(ktlMolecule &molNew,int &im){
 double contactPredPenTotal=0.0;
  for(int i=0;i<contactPredPen.size();i++){
    if(i==im){
      contactPredPen[i] = molNew.getLennardJonesContact();
      perPairPenaltiesPerMol[i] = molNew.getSoftContactPenaltiesPerPair();
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }else{
      contactPredPenTotal=contactPredPenTotal+contactPredPen[i];
    }
  }
  contactPredPenTotal += ensembleOrPenalty();
  return contactPredPenTotal;
}

// hard constraints are a feasibility filter, not a penalty -- these just (re)compute and
// cache per-molecule pass/fail + worst violation, mirroring applyDistanceConstraints'
// per-molecule caching pattern above.
void moleculeFitAndState::applyHardConstraints(){
  hardOKPerMol.assign(mol.size(), true);
  hardMaxViolationPerMol.assign(mol.size(), 0.0);
  for(int i=0;i<mol.size();i++){
    double v = 0.0;
    hardOKPerMol[i] = mol[i].hardConstraintsSatisfied(v);
    hardMaxViolationPerMol[i] = v;
  }
}

void moleculeFitAndState::applyHardConstraints(ktlMolecule &molNew,int &im){
  if(hardOKPerMol.size() != mol.size()){
    hardOKPerMol.assign(mol.size(), true);
    hardMaxViolationPerMol.assign(mol.size(), 0.0);
  }
  double v = 0.0;
  hardOKPerMol[im] = molNew.hardConstraintsSatisfied(v);
  hardMaxViolationPerMol[im] = v;
}

bool moleculeFitAndState::getHardConstraintsSatisfied(){
  for(int i=0;i<hardOKPerMol.size();i++){
    if(!hardOKPerMol[i]){ return false; }
  }
  return true;
}

double moleculeFitAndState::getHardConstraintsMaxViolation(){
  double m = 0.0;
  for(int i=0;i<hardMaxViolationPerMol.size();i++){
    if(hardMaxViolationPerMol[i] > m){ m = hardMaxViolationPerMol[i]; }
  }
  return m;
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

// Ensemble writhe-difference restraint: keeps mixture-ensemble members from
// diverging from each other into "one wild conformer, the rest untouched"
// rather than a coherent set of alternate states. Same soft-sigmoid pattern
// as applyWritheConstraint() above (added into currFit, not a hard reject),
// so a too-tight threshold makes divergence costly rather than stalling the
// search outright. maxWritheDiff <= 0 (the default) makes this a no-op --
// existing single-state runs and any caller not passing argv[20] are
// completely unaffected. The sigmoid's transition width is set relative to
// the threshold itself (rather than a separate hardcoded constant) so it
// stays sensible across a wide range of --max_writhe_diff choices; this is
// an easy first thing to hand-tune later if it doesn't behave as it should.
void moleculeFitAndState::applyWritheDiffConstraint(ktlMolecule &molNew,int &i){
  writheDiffPenalty = 0.0;
  if(maxWritheDiff <= 0.0 || mol.size() < 2) return;
  writheFP wfp;
  double steepness = 8.0/std::max(maxWritheDiff,1e-6);
  for(int k=0;k<mol.size();k++){
    if(k==i) continue;
    double diff = 0.0;
    int nCh = std::min(molNew.noChains(),mol[k].noChains());
    for(int c=0;c<nCh;c++){
      std::vector<point> coordsA = molNew.getCoordinatesSection(c);
      std::vector<point> coordsB = mol[k].getCoordinatesSection(c);
      std::vector<point> listA = wfp.strideList(coordsA,writheDiffStride);
      std::vector<point> listB = wfp.strideList(coordsB,writheDiffStride);
      diff = diff + wfp.writheMatrixAbsDiff(listA,listB);
    }
    writheDiffPenalty = writheDiffPenalty + 1.0/(1.0+std::exp(-steepness*(diff-maxWritheDiff)));
  }
}

// Baseline (no trial move yet) -- sums over every pair of ensemble members
// rather than just "moved state vs the rest", so the very first currFit is
// built from the same set of terms the per-move path uses afterwards.
void moleculeFitAndState::applyWritheDiffConstraint(){
  writheDiffPenalty = 0.0;
  if(maxWritheDiff <= 0.0 || mol.size() < 2) return;
  writheFP wfp;
  double steepness = 8.0/std::max(maxWritheDiff,1e-6);
  for(int a=0;a<mol.size();a++){
    for(int b=a+1;b<mol.size();b++){
      double diff = 0.0;
      int nCh = std::min(mol[a].noChains(),mol[b].noChains());
      for(int c=0;c<nCh;c++){
        std::vector<point> coordsA = mol[a].getCoordinatesSection(c);
        std::vector<point> coordsB = mol[b].getCoordinatesSection(c);
        std::vector<point> listA = wfp.strideList(coordsA,writheDiffStride);
        std::vector<point> listB = wfp.strideList(coordsB,writheDiffStride);
        diff = diff + wfp.writheMatrixAbsDiff(listA,listB);
      }
      writheDiffPenalty = writheDiffPenalty + 1.0/(1.0+std::exp(-steepness*(diff-maxWritheDiff)));
    }
  }
}

double moleculeFitAndState::getWritheDiffPenalty(){
  return writheDiffPenalty;
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



// Baseline (no trial move): full recompute across every mol[i]. Replaces the 4 no-arg
// getOverallFit*/getOverallFitForceConnection*[_ChiSq] functions -- those differed only in
// (a) which scattering function to call and (b) whether to add connectionPenalty, both now
// selected by `mode` instead of by which differently-named function you happened to call.
//
// Consolidating surfaced two real inconsistencies in the old 8-function version, fixed here
// as a natural consequence of there now being only one formula:
//   - the old non-ChiSq ForceConnection variant omitted overlapPenalty from currFit (every
//     other one of the 8 included it) -- currFit below always includes it.
//   - helpers.cpp's increaseKmax() used to hardcode the plain (non-ChiSq, non-ForceConnection)
//     variant regardless of the run's actual mode; it now takes a FitMode and passes the
//     run's real mode through.
std::pair<double,double> moleculeFitAndState::computeOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax,const FitMode &mode){
  double scatterAndHydrationConstraint = mode.weightedChiSq
    ? ed.calculateChiSquared_Weighted(mol,kmin,kmax,mixtureList)
    : ed.calculateChiSquared(mol,kmin,kmax,mixtureList);
  /***************************************************************

   apply penalties which are "un protein like". Currently we are using

     i) a very strict overlap penalty which exponetiallp penalises non local sections coming close than 4 A.
     ii) A distance constraint measure, which is only active if the user inputs a set of distance consrtrainst like contact predictions.
     iii) A writhe penalty to ensure the moelule doesn't become too disentangled.

  **************************************************************/
  double overlapPenalty = applyOverlapPenalty();
  double distanceConstraints = applyDistanceConstraints();
  applyHardConstraints();
  applyWritheConstraint();
  applyWritheDiffConstraint();

  double connectionTerm = 0.0;
  if(mode.forceConnection){
    for(int i=0;i<mol.size();i++){
      calculateConnectionPenalty(mol[i],i);
    }
    connectionTerm = connectionPenalty;
  }

  double currFit = scatterAndHydrationConstraint * (1.0 + penaltyWeight*(distanceConstraints + writhePenalty + overlapPenalty + writheDiffPenalty + connectionTerm));
  std::pair<double,double> fitStats;
  fitStats.first = currFit;
  fitStats.second = scatterAndHydrationConstraint;

  // pass along best hydration C2 parameter
  return fitStats;
}

// Trial move: incremental recompute for the single changed molecule i. Replaces the 4
// (molNew, i) getOverallFit*/getOverallFitForceConnection*[_ChiSq] functions -- same two
// axes (weightedChiSq, forceConnection), same consolidation.
std::pair<double,double> moleculeFitAndState::computeOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i,const FitMode &mode){
  calculateMoleculeDistances(molNew,i);
  double scatterAndHydrationConstraint = mode.weightedChiSq
    ? ed.calculateChiSquaredUpdate_Weighted(molNew,i,kmin,kmax,mixtureList)
    : ed.calculateChiSquaredUpdate(molNew,i,kmin,kmax,mixtureList);

  double overlapPenalty = applyOverlapPenalty();
  double distanceConstraints = applyDistanceConstraints(molNew,i);
  applyHardConstraints(molNew,i);
  alterWritheSet(molNew,i);
  applyWritheConstraint();
  applyWritheDiffConstraint(molNew,i);

  double connectionTerm = 0.0;
  if(mode.forceConnection){
    calculateConnectionPenalty(molNew,i);
    connectionTerm = connectionPenalty;
  }

  double currFit = scatterAndHydrationConstraint * (1.0 + penaltyWeight*(overlapPenalty + distanceConstraints + writhePenalty + writheDiffPenalty + connectionTerm));
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
