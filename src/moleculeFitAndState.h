#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "writheFP.h"
#include "experimentalData.h"
#include "parameters.h"

#include <vector>
#include <utility>

class moleculeFitAndState {
public:
  // best fitting hydration parameter (unchanged usage in your codebase)
  double C2 = 0.0;

  // lifecycle
  moleculeFitAndState(const std::vector<ktlMolecule>& mol, ModelParameters& params);
  ~moleculeFitAndState() = default;

  moleculeFitAndState(const moleculeFitAndState&) = default;
  moleculeFitAndState& operator=(const moleculeFitAndState&) = default;

  moleculeFitAndState(moleculeFitAndState&&) noexcept = default;
  moleculeFitAndState& operator=(moleculeFitAndState&&) noexcept = default;

  // molecule access
  //const std::vector<ktlMolecule>& getMolecule() const { return mol; }
  void updateMolecule(const std::vector<ktlMolecule>& molNew);  // copy
  void updateMolecule(std::vector<ktlMolecule>&& molNew);        // move

  // IO
  void writeScatteringToFile(experimentalData& ed, std::vector<std::vector<double>>& mixtureVals, const char* filename);
  void writeScatteringToFile_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureVals, const char* filename);

  // penalties / constraints
  double getOverlapPenalty(double& closestApproachDist, std::vector<double>& overlapDists);
  double applyOverlapPenalty();
  double applyDistanceConstraints();
  double applyDistanceConstraints(ktlMolecule& molNew, int i);
  void   applyWritheConstraint();

  // per-molecule updates
  void   calculateMoleculeDistances(ktlMolecule& molNew, int i);
  void   calculateConnectionPenalty(ktlMolecule& molNew, int chInd);
  void   alterWritheSet(ktlMolecule& molNew, int i);

  // user-specified connectivity
  double calculateUserSpecifiedConnectionPenalty(ktlMolecule& molNew, std::vector<int>& chainSet1, std::vector<int>& chainSet2);
  double calculateUserSpecifiedConnectionPenalty(int chInd, std::vector<int>& chainSet1, std::vector<int>& chainSet2);

  // reporting
  double getWrithePenalty() { return writhePenalty; }
  double getOverlapPenalty() { return applyOverlapPenalty(); }
  double getDistanceConstraints() { return applyDistanceConstraints(); }
  double getBetaSheetReward();

  // overall fit APIs (unchanged names; some ints by value)
  std::pair<double,double> getOverallFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);
  std::pair<double,double> getOverallFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);
  std::pair<double,double> getOverallFitForceConnection(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);
  std::pair<double,double> getOverallFitForceConnection_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);

  std::pair<double,double> getOverallFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i);
  std::pair<double,double> getOverallFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i);
  std::pair<double,double> getOverallFitForceConnection(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i);
  std::pair<double,double> getOverallFitForceConnection_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i);

  void updateScatteringFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);
  void updateScatteringFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax);

  double currFit = 0.0;

  // memory hygiene (optional)
  void shrink_temporaries();
  void replaceMoleculeAt(int idx, const ktlMolecule& m) {
    if (idx >= 0 && idx < (int)mol.size()) mol[idx] = m;
  };
  const std::vector<ktlMolecule>& getMolecule() const { return mol; }; // keep const

private:
  // cached per-molecule info
  std::vector<std::vector<double>> molDists;
  std::vector<std::vector<double>> solDists;
  std::vector<std::vector<double>> solMolDists;
  std::vector<std::vector<double>> overlapDistSet;
  std::vector<int>                 molSize;
  std::vector<int>                 noSol;

  double maxDist = 0.0;
  double hydroPhobicPacking = 0.0;

  std::vector<std::vector<double>> originalWrithes;
  std::vector<std::vector<double>> currWrithes;

  std::vector<double> maxDistMol;
  std::vector<double> maxDistSol;
  std::vector<double> contactPredPen;

  double writhePenalty = 0.0;
  double originalOverlapPenalty = 0.0;

  double Rin = 0.0, Rout = 0.0, RShell = 0.0, ntrivs = 0.0, closestApproachDist = 0.0;
  double solventsPerLink = 0.0, rmin = 0.0, rmax = 0.0, lmin = 0.0;

  std::vector<double> percentageCombinations;

  std::vector<ktlMolecule> mol;

  double connectionPenalty = 0.0;
  std::vector<double> connectionPenaltySet;

  // scratch buffer to avoid per-call allocations when we must create sub-sections
  std::vector<std::vector<point>> scratchSubcoords;
};

#endif
