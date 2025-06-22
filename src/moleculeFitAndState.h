#ifndef ALG_ROUTE
#define ALG_ROUTE

#include "ktlMoleculeRandom.h"
#include "writheFP.h"
#include "experimentalData.h"
#include "parameters.h"

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
  void calculateMoleculeDistances(ktlMolecule &molNew,int &i);
  double calculateUserSpecifiedConnectionPenalty(ktlMolecule &molNew,std::vector<int> &chainSet1,std::vector<int> &chainSet2);
  double calculateUserSpecifiedConnectionPenalty(int chInd,std::vector<int> &chainSet1,std::vector<int> &chainSet2);
  void applyWritheConstraint();
  void calculateConnectionPenalty(ktlMolecule &molNew, int &chInd);
  double getWrithePenalty();
  double getOverlapPenalty();
  double getDistanceConstraints();
  void alterWritheSet(ktlMolecule &molNew,int &i);
  std::pair<double,double> getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  std::pair<double,double> getOverallFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  std::pair<double,double>  getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  std::pair<double,double>  getOverallFitForceConnection_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,double &kmin,double &kmax);
  std::pair<double,double> getOverallFit(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  std::pair<double,double> getOverallFit_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  std::pair<double,double>  getOverallFitForceConnection(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
  std::pair<double,double>  getOverallFitForceConnection_ChiSq(experimentalData &ed,std::vector<std::vector<double> > &mixtureList,ktlMolecule &molNew,double &kmin,double &kmax,int &i);
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
  std::vector<double> maxDistMol;
  std::vector<double> maxDistSol;
  std::vector<double> contactPredPen;
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
