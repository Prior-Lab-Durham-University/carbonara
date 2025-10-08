#include "moleculeFitAndState.h"
#include "ktlMoleculeRandom.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// ctor
// -----------------------------------------------------------------------------
moleculeFitAndState::moleculeFitAndState(const std::vector<ktlMolecule>& molin, ModelParameters& params)
: mol(molin)
{
  // fixed fitting params
  closestApproachDist = params.closestApproachDist;
  rmin = params.rmin; rmax = params.rmax; lmin = params.lmin;

  // size-dependent containers
  const int M = static_cast<int>(mol.size());
  molDists.resize(M);
  molSize.resize(M);
  maxDistMol.resize(M);
  contactPredPen.resize(M);
  connectionPenaltySet.assign(M, 0.0);
  writhePenalty = 0.0;
  connectionPenalty = 0.0;

  // precompute writhe per chain (reuse scratchSubcoords per chain to reduce churn)
  for (int i = 0; i < M; ++i) {
    writheFP wfp;
    std::vector<double> chainWrithes;
    chainWrithes.reserve(mol[i].noChains());

    const auto& coords    = mol[i].getCoordinates();
    const auto& chainList = mol[i].getChainList();

    for (int j = 0; j < mol[i].noChains(); ++j) {
      const int first = chainList[j].first;
      const int last  = chainList[j].second;

      scratchSubcoords.clear();
      scratchSubcoords.insert(scratchSubcoords.end(), coords.begin() + first, coords.begin() + last + 1);

      double wr = wfp.DIDownSampleAbsSingle(scratchSubcoords);
      chainWrithes.push_back(wr);
    }
    originalWrithes.push_back(chainWrithes);
  }
  currWrithes = originalWrithes;

  // distance/overlap baselines per molecule
  for (int i = 0; i < M; ++i) {
    std::vector<double> overlapDists = mol[i].checkOverlapWithRad(closestApproachDist);
    molDists[i]  = mol[i].getDistSet(); // copy from const-ref is fine (you cache your own bins)
    molSize[i]   = mol[i].getNoAminos();

    std::sort(molDists[i].begin(), molDists[i].end());
    if (!molDists[i].empty()) {
      maxDistMol[i] = molDists[i].back();
    } else {
      maxDistMol[i] = 0.0;
    }

    overlapDistSet.push_back(std::move(overlapDists));

    // connectivity penalty baseline
    std::vector<double> minIntraChainDistances = mol[i].getMinimumintraMoleculaeDistancePerChain();
    for (int j = 0; j < static_cast<int>(minIntraChainDistances.size()); ++j) {
      const double md = minIntraChainDistances[j];
      if (md < 9999.9 && md > 5.0) {
        const double dif = md - 5.0;
        connectionPenaltySet[i] += 0.0005 * std::pow(dif, 10.0);
      }
    }
    connectionPenaltySet[i] /= std::max(1, mol[i].noChains());
    connectionPenalty += connectionPenaltySet[i];
  }

  originalOverlapPenalty = 0.0;
  const double overlapPenalty = applyOverlapPenalty();
  originalOverlapPenalty = overlapPenalty;
}

// -----------------------------------------------------------------------------
// molecule access
// -----------------------------------------------------------------------------
void moleculeFitAndState::updateMolecule(const std::vector<ktlMolecule>& molNew) {
  mol = molNew;
}
void moleculeFitAndState::updateMolecule(std::vector<ktlMolecule>&& molNew) {
  mol = std::move(molNew);
}

// -----------------------------------------------------------------------------
// IO
// -----------------------------------------------------------------------------
void  moleculeFitAndState::writeScatteringToFile(experimentalData& ed, std::vector<std::vector<double>>& mixtureVals, const char* filename) {
  ed.writeScatteringToFile(mixtureVals, filename);
}
void  moleculeFitAndState::writeScatteringToFile_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureVals, const char* filename) {
  ed.writeScatteringToFile_ChiSq(mixtureVals, filename);
}

// -----------------------------------------------------------------------------
// penalties / constraints
// -----------------------------------------------------------------------------
double moleculeFitAndState::getOverlapPenalty(double& closestApproachDist_, std::vector<double>& overlapDists) {
  double distSumCurr = 0.0;
  for (int l = 0; l < static_cast<int>(overlapDists.size()); ++l) {
    const double dist = closestApproachDist_ - overlapDists[l];
    distSumCurr += std::exp(dist) - 1.0;
  }
  if (!overlapDists.empty()) distSumCurr = 0.01 * distSumCurr / overlapDists.size();
  return distSumCurr;
}

double moleculeFitAndState::applyOverlapPenalty() {
  double overlapPenalty = 0.0;
  for (int i = 0; i < static_cast<int>(overlapDistSet.size()); ++i) {
    overlapPenalty += getOverlapPenalty(closestApproachDist, overlapDistSet[i]);
  }
  return (overlapPenalty > originalOverlapPenalty) ? (overlapPenalty - originalOverlapPenalty) : 0.0;
}

double moleculeFitAndState::applyDistanceConstraints() {
  double contactPredPenTotal = 0.0;
  for (int i = 0; i < static_cast<int>(mol.size()); ++i) {
    contactPredPen[i] = mol[i].getLennardJonesContact();
    contactPredPenTotal += contactPredPen[i];
  }
  return contactPredPenTotal;
}

double moleculeFitAndState::applyDistanceConstraints(ktlMolecule& molNew, int im) {
  double contactPredPenTotal = 0.0;
  for (int i = 0; i < static_cast<int>(contactPredPen.size()); ++i) {
    if (i == im) {
      contactPredPen[i] = molNew.getLennardJonesContact();
    }
    contactPredPenTotal += contactPredPen[i];
  }
  return contactPredPenTotal;
}

void moleculeFitAndState::applyWritheConstraint() {
  writhePenalty = 0.0;
  for (int i = 0; i < static_cast<int>(currWrithes.size()); ++i) {
    for (int j = 0; j < static_cast<int>(currWrithes[i].size()); ++j) {
      const double newWrithe = currWrithes[i][j];
      const double secLen = double(mol[i].getSubsecSize(j + 1)); // getSubsecSize: 1-based
      const double lowerBound = std::pow((secLen / 7.5), 1.6) - 3.0;
      writhePenalty += 1.0 / (1.0 + std::exp(20.0 * (newWrithe - lowerBound)));
    }
  }
}

// -----------------------------------------------------------------------------
// per-molecule updates
// -----------------------------------------------------------------------------
void moleculeFitAndState::calculateMoleculeDistances(ktlMolecule& molNew, int i) {
  std::vector<double> overlapDists = molNew.checkOverlapWithRad(closestApproachDist);
  molDists[i] = molNew.getDistSet(); // copy cached set
  molSize[i]  = molNew.getNoAminos();

  std::sort(molDists[i].begin(), molDists[i].end());
  maxDistMol[i] = molDists[i].empty() ? 0.0 : molDists[i].back();

  overlapDistSet[i] = std::move(overlapDists);
}

void moleculeFitAndState::alterWritheSet(ktlMolecule& molNew, int i) {
  writheFP wfp;
  const auto& coords    = molNew.getCoordinates();
  const auto& chainList = molNew.getChainList();

  // (re)size currWrithes[i] if needed
  if (i >= static_cast<int>(currWrithes.size())) currWrithes.resize(i + 1);
  currWrithes[i].resize(molNew.noChains());

  for (int j = 0; j < molNew.noChains(); ++j) {
    const int first = chainList[j].first;
    const int last  = chainList[j].second;

    scratchSubcoords.clear();
    scratchSubcoords.insert(scratchSubcoords.end(), coords.begin() + first, coords.begin() + last + 1);

    const double wr = wfp.DIDownSampleAbsSingle(scratchSubcoords);
    currWrithes[i][j] = wr;
  }
}

void moleculeFitAndState::calculateConnectionPenalty(ktlMolecule& molNew, int chInd) {
  connectionPenaltySet[chInd] = 0.0;
  std::vector<double> minIntraChainDistances = molNew.getMinimumintraMoleculaeDistancePerChain();

  for (int j = 0; j < static_cast<int>(minIntraChainDistances.size()); ++j) {
    const double md = minIntraChainDistances[j];
    if (md < 9999.9 && md > 5.0) {
      const double dif = md - 5.0;
      connectionPenaltySet[chInd] += 0.0005 * std::pow(dif, 10.0);
    }
  }
  connectionPenaltySet[chInd] /= std::max(1, mol[chInd].noChains());

  connectionPenalty = 0.0;
  for (double v : connectionPenaltySet) connectionPenalty += v;
}

// -----------------------------------------------------------------------------
// user specified connectivity
// -----------------------------------------------------------------------------
double moleculeFitAndState::calculateUserSpecifiedConnectionPenalty(ktlMolecule& molNew, std::vector<int>& chainSet1, std::vector<int>& chainSet2) {
  std::vector<std::vector<double>> minIntra = molNew.getMinimumintraMolecularDistances();
  double userConnectedPenalty = 0.0;
  double minDist = 10000.0;

  for (int a : chainSet1) {
    for (int b : chainSet2) {
      if (minIntra[a][b] < minDist) minDist = minIntra[a][b];
    }
  }
  if (minDist < 9999.0 && minDist > 5.0) {
    const double dif = minDist - 5.0;
    userConnectedPenalty += 0.0005 * std::pow(dif, 10.0);
  }
  return userConnectedPenalty;
}

double moleculeFitAndState::calculateUserSpecifiedConnectionPenalty(int chInd, std::vector<int>& chainSet1, std::vector<int>& chainSet2) {
  std::vector<std::vector<double>> minIntra = mol[chInd].getMinimumintraMolecularDistances();
  double userConnectedPenalty = 0.0;
  double minDist = 10000.0;

  for (int a : chainSet1) {
    for (int b : chainSet2) {
      if (minIntra[a][b] < minDist) minDist = minIntra[a][b];
    }
  }
  if (minDist < 9999.0 && minDist > 5.0) {
    const double dif = minDist - 5.0;
    userConnectedPenalty += 0.0005 * std::pow(dif, 10.0);
  }
  return userConnectedPenalty;
}

// -----------------------------------------------------------------------------
// overall fits
// -----------------------------------------------------------------------------
std::pair<double,double> moleculeFitAndState::getOverallFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  double scatter = ed.calculateChiSquared(mol, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints();
  applyWritheConstraint();

  const double curr = scatter + 1.0 * (distanceConstraints + writhePenalty + overlapPenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  double scatter = ed.calculateChiSquared_Weighted(mol, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints();
  applyWritheConstraint();

  const double curr = scatter + 100.0 * (distanceConstraints + writhePenalty + overlapPenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  double scatter = ed.calculateChiSquared(mol, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints();
  applyWritheConstraint();

  for (int i = 0; i < static_cast<int>(mol.size()); ++i) calculateConnectionPenalty(mol[i], i);

  const double curr = scatter + 1.0 * (distanceConstraints + writhePenalty + connectionPenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  double scatter = ed.calculateChiSquared_Weighted(mol, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints();
  applyWritheConstraint();

  for (int i = 0; i < static_cast<int>(mol.size()); ++i) calculateConnectionPenalty(mol[i], i);

  const double curr = scatter + 100.0 * (distanceConstraints + overlapPenalty + writhePenalty + connectionPenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i) {
  calculateMoleculeDistances(molNew, i);
  double scatter = ed.calculateChiSquaredUpdate(molNew, i, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints(molNew, i);
  alterWritheSet(molNew, i);
  applyWritheConstraint();

  const double curr = scatter + 1.0 * (overlapPenalty + distanceConstraints + writhePenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i) {
  calculateMoleculeDistances(molNew, i);
  double scatter = ed.calculateChiSquaredUpdate_Weighted(molNew, i, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints(molNew, i);
  alterWritheSet(molNew, i);
  applyWritheConstraint();

  const double curr = scatter + 100.0 * (overlapPenalty + distanceConstraints + writhePenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i) {
  calculateMoleculeDistances(molNew, i);
  double scatter = ed.calculateChiSquaredUpdate(molNew, i, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints(molNew, i);
  alterWritheSet(molNew, i);
  applyWritheConstraint();
  calculateConnectionPenalty(molNew, i);

  const double curr = scatter + 1.0 * (overlapPenalty + distanceConstraints + writhePenalty + connectionPenalty);
  return {curr, scatter};
}

std::pair<double,double> moleculeFitAndState::getOverallFitForceConnection_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, ktlMolecule& molNew, double& kmin, double& kmax, int i) {
  calculateMoleculeDistances(molNew, i);
  double scatter = ed.calculateChiSquaredUpdate_Weighted(molNew, i, kmin, kmax, mixtureList);

  const double overlapPenalty      = applyOverlapPenalty();
  const double distanceConstraints = applyDistanceConstraints(molNew, i);
  alterWritheSet(molNew, i);
  applyWritheConstraint();
  calculateConnectionPenalty(molNew, i);

  const double curr = scatter + 100.0 * (overlapPenalty + distanceConstraints + writhePenalty + connectionPenalty);
  return {curr, scatter};
}

// -----------------------------------------------------------------------------
// scattering updates
// -----------------------------------------------------------------------------
void moleculeFitAndState::updateScatteringFit(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  for (int i = 0; i < static_cast<int>(mol.size()); ++i) calculateMoleculeDistances(mol[i], i);
  const double dummyScatter = ed.calculateChiSquared(mol, kmin, kmax, mixtureList);
  std::cout << "double check baseLine " << dummyScatter << "\n";
}

void moleculeFitAndState::updateScatteringFit_ChiSq(experimentalData& ed, std::vector<std::vector<double>>& mixtureList, double& kmin, double& kmax) {
  for (int i = 0; i < static_cast<int>(mol.size()); ++i) calculateMoleculeDistances(mol[i], i);
  (void)ed.calculateChiSquared_Weighted(mol, kmin, kmax, mixtureList);
}

// -----------------------------------------------------------------------------
// beta sheet reward (unchanged math; minor cleanup)
// -----------------------------------------------------------------------------
double moleculeFitAndState::getBetaSheetReward() {
  double sheetRewards = 0.0;
  for (int i = 0; i < static_cast<int>(mol.size()); ++i) {
    const double numSheets = mol[i].numBetaSheets;
    sheetRewards += (numSheets > 0.0) ? (mol[i].getBetaSheetProximityReward() / numSheets) : 0.0;
  }
  return sheetRewards;
}

// -----------------------------------------------------------------------------
// memory hygiene
// -----------------------------------------------------------------------------
void moleculeFitAndState::shrink_temporaries() {
  for (auto& v : molDists)          v.shrink_to_fit();
  for (auto& v : solDists)          v.shrink_to_fit();
  for (auto& v : solMolDists)       v.shrink_to_fit();
  for (auto& v : overlapDistSet)    v.shrink_to_fit();
  maxDistMol.shrink_to_fit();
  maxDistSol.shrink_to_fit();
  contactPredPen.shrink_to_fit();
  molSize.shrink_to_fit();
  noSol.shrink_to_fit();
  connectionPenaltySet.shrink_to_fit();
  scratchSubcoords.shrink_to_fit();
}
