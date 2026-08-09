#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "writheFP.h"
#include "experimentalData.h"
#include "parameters.h"

// Which fitting procedure to use for a currFit computation, replacing what used to be
// 8 separately-named, hand-duplicated functions (every combination of these two flags,
// crossed with baseline-vs-trial-move). Select the procedure at the call site by passing
// a FitMode instead of picking a differently-named function -- adding a new fitting
// behaviour in future is "add a field here", not "copy-paste a 9th/10th function and hope
// every call site remembers to use it".
struct FitMode {
  // true = use the weighted chi-squared scattering calculation (calculateChiSquared_Weighted
  // / calculateChiSquaredUpdate_Weighted) -- was the "_ChiSq" suffixed functions.
  bool weightedChiSq = true;
  // true = also add the inter-chain connection penalty to currFit (calculateConnectionPenalty)
  // -- was the "ForceConnection" functions.
  bool forceConnection = false;
};

class moleculeFitAndState{
public:

  // best fitting hydration parameter
  double C2;

  moleculeFitAndState(std::vector<ktlMolecule> &mol, ModelParameters& params);
  std::vector<ktlMolecule> getMolecule();
  void updateMolecule(std::vector<ktlMolecule> &molNew);
  void writeScatteringToFile(experimentalData &ed,std::vector<std::vector<double> > &mixtureVals,const char* filename);
  void writeScatteringToFile_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureVals,const char* filename);
  double getOverlapPenalty(double &closestApproachDist,std::vector<double> &overlapDists);
  double applyOverlapPenalty();
  double applyDistanceConstraints();
  double applyDistanceConstraints(ktlMolecule &molNew,int &i);
  // sums, across mixture states, the per-pair minimum penalty among ensemble-OR-flagged
  // pairs -- "satisfied by any one conformation" rather than "every conformation must
  // satisfy it". See ktlMolecule::getSoftContactPenaltiesPerPair.
  double ensembleOrPenalty();
  // hard constraints (disulfides, strict posts, etc.) -- true feasibility filter, not a
  // penalty. applyHardConstraints(...) (re)computes and caches pass/fail + worst violation
  // per molecule, mirroring applyDistanceConstraints' caching pattern; the getters below
  // read that cache. Caller must check getHardConstraintsSatisfied() alongside chi2
  // improvement before accepting a move -- no chi2 gain can substitute for it.
  void applyHardConstraints();
  void applyHardConstraints(ktlMolecule &molNew,int &i);
  bool getHardConstraintsSatisfied();
  double getHardConstraintsMaxViolation();
  void calculateMoleculeDistances(ktlMolecule &molNew,int &i);
  double calculateUserSpecifiedConnectionPenalty(ktlMolecule &molNew,std::vector<int> &chainSet1,std::vector<int> &chainSet2);
  double calculateUserSpecifiedConnectionPenalty(int chInd,std::vector<int> &chainSet1,std::vector<int> &chainSet2);
  void applyWritheConstraint();
  void calculateConnectionPenalty(ktlMolecule &molNew, int &chInd);
  double getWrithePenalty();
  // Ensemble writhe-difference restraint (mixture_states > 1 fitting only;
  // a no-op whenever maxWritheDiff <= 0, which is the default). The no-arg
  // overload is for the one-off baseline fit (sums over every pair of
  // current ensemble members); the (molNew, i) overload is what the
  // per-trial-move path uses (molNew, the candidate for state i, against
  // every other *unchanged* mol[k]).
  void applyWritheDiffConstraint();
  void applyWritheDiffConstraint(ktlMolecule &molNew, int &i);
  double getWritheDiffPenalty();
  double getOverlapPenalty();
  double getDistanceConstraints();
  void alterWritheSet(ktlMolecule &molNew,int &i);
  // baseline: full recompute across every mol[i]. Replaces getOverallFit / getOverallFit_ChiSq
  // / getOverallFitForceConnection / getOverallFitForceConnection_ChiSq (no-arg versions).
  std::pair<double,double> computeOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax,const FitMode &mode);
  // trial move: incremental recompute for the single changed molecule i. Replaces the
  // (molNew, i) versions of the same 4 functions.
  std::pair<double,double> computeOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i,const FitMode &mode);
  void updateScatteringFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  double currFit;
  void updateScatteringFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  double getBetaSheetReward();



private:
  std::vector<std::vector<double> > molDists;
  std::vector<std::vector<double> > solDists;
  std::vector<std::vector<double> > solMolDists;
  std::vector<std::vector<double> > overlapDistSet;
  std::vector<int> molSize;
  std::vector<int> noSol;
  double maxDist;
  double hydroPhobicPacking;
  std::vector<std::vector<double> > originalWrithes;
  std::vector<std::vector<double> > currWrithes;
  double maxWritheDiff;
  int writheDiffStride;
  double writheDiffPenalty;
  double penaltyWeight;
  std::vector<double> maxDistMol;
  std::vector<double> maxDistSol;
  std::vector<double> contactPredPen;
  // per mixture state, per constraint pair index (same order as loaded contactPairList):
  // capped soft penalty, populated only for ensemble-OR-flagged pairs (0.0 otherwise).
  // See ensembleOrPenalty.
  std::vector<std::vector<double>> perPairPenaltiesPerMol;
  std::vector<bool> hardOKPerMol;
  std::vector<double> hardMaxViolationPerMol;
  double writhePenalty;
  double originalOverlapPenalty;
  double Rin,Rout,RShell,ntrivs,closestApproachDist;
  double solventsPerLink,rmin,rmax,lmin;
  std::vector<double> percentageCombinations;
  std::vector<ktlMolecule> mol;
  double connectionPenalty;
  std::vector<double> connectionPenaltySet;
};

#endif
