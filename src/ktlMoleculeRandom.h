#ifndef KTL_MOL
#define KTL_MOL

#include "point.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <dirent.h>
#include "polyHelix.h"
#include <random>
#include "randomMolGen.h"
#include <tuple>

class ktlMolecule{
public:
  ktlMolecule();
  //ktlMolecule(ktlMolecule &ktlm);
  void setParams(double &rminIn,double &rmaxIn,double &lminIn);
  std::vector<int> getUnitNos();
  std::vector<std::vector<point> > getTangents();
  std::vector<std::vector<point> > getNormals();
  std::vector<std::vector<point> > getBinormals();
  std::vector<std::vector<point> > getCoordinates();
  std::vector<std::vector<std::string> > getAminoList();
  std::vector<point> getCoordinatesSection(int i);
  int getSubsecSize(int sec);
  std::vector<std::vector<point> > getSubsecCoordinates(int &sec);
  std::vector<std::pair<std::string,int> > getNameSizeListOfSection(int &sec);
  //std::vector<point> getEulerAngles();
  std::vector<double> getDistChanges();
  double getCurvatureJoined(int index);
  double getTorsionJoined(int index);
  double getLengthJoined(int index);
  int getUnitNo(int index);
  double maxNeighbourDistSec(int &sec);
  point getTangent(int mindex,int subindex);
  point getNormal(int mindex,int subindex);
  point getBinormal(int mindex,int subindex);
  point getCoordinate(int mindex,int subindex);
  double getAlbeadJoined(int index);
  std::string getType(int &index);
  std::string getType(int &chainNo,int &index);
  //double getDistChange(int index);
  double getMaxDistChange();
  std::pair<double,double> getMaxPossibleLength();
  int molSize();
  int noChains();
  int noSecSize();
  int getNoAminos();
  void createSyntheticMolecule(int nLinks);
  void readInSequence(const char* filename,double &rmin,double &rmax,double &lmin);
  void readInCoordinates(const char* filename);
  void getHydrophobicResidues();
  std::vector<double> getHydrophobicDistance(std::vector<std::vector<point> > &solventList,double &maxSolDist);
  void getCoiledCoilResidues();
  void getPositiveResidues();
  void getNegativeResidues();
  void getPolarResidues();
  double coiledCoilPotential();
  double coiledCoilPotentialBetween(int &secNo);
  double coiledCoilPotentialBetween();
  point getCentreOfMass(std::vector<std::vector<point> > &cdSet);
  std::vector<int> checkOverlap(std::vector<std::vector<point> > &cdsIN);
  std::vector<double> checkOverlapWithRad(double &wRad,int &sec);
  std::vector<double> checkOverlapWithRad(double &wRad);
  double getMininumIntraMoelcularDistance();
  std::vector<double> getMinimumintraMoleculaeDistancePerChain();
  std::vector<std::vector<double> > getMinimumintraMolecularDistances();
  std::vector<double> getDistSet();
  double compareDistances(std::vector<std::vector<point> > &coords2);
  bool checkCalphas(std::vector<std::vector<point> > &coordsIn);
  bool checkCalphas();
  bool checkCalphas(int &index);
  bool checkCalphas(int &secNo,ktlMolecule &ktl);
  void changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn,std::vector<std::pair<std::string,int> > &nameSizeSubList);
  int getRandomMolecule();
  int getRandomMoleculeReset();
  void resetRandomMolecule();
  void changeMoleculeSingleMulti(int &index,int sec);
  void changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec);
  void replicateMolecule(int &noReplications);
  void rotation3D(point &p,point &centre,point &k,double &cosangle,double &sinangle);
  void rotateSection(std::vector<std::vector<point> >  &section,point &centre,point &k,double &angle,point &transVec);
  void writeMoleculeToFile(const char* filename);
  std::vector<std::pair<double,double> > getKapTauVals();
  double getPairDistance(std::pair<int,int> &index1,std::pair<int,int> &index2);
  double lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff);
  std::vector<double> solMolDists(std::vector<std::vector<point> > &pts1);
  void loadContactPredictions(const char* contactloc);
  double getLennardJonesContact();
  void loadFixedSections(const char* fixedsecloc);

  double getBetaSheetProximityReward();
  bool areSheetsPairallelOrAntiparallel(std::vector<point>& sheet1, std::vector<point>& sheet2);
  double numBetaSheets;

private:
  std::vector<double> minimumintraMoleculaeDistancePerChain;
  std::vector<std::vector<double> > minimumintraMolecularDistances;
  double minimumintraMolecularDistanceMean;
  std::vector<int> noPts;
  std::vector<std::vector<point> > coords;
  std::vector<std::vector<std::string> > aminoList;
  std::vector<std::vector<point> > tanlist;
  std::vector<std::vector<point> > normlist;
  std::vector<std::vector<point> > binormlist;
  polyHelix ph;
  std::vector<std::pair<int,int> > chainList;
  std::vector<double> distChanges;
  std::vector<std::pair<std::string,int> > nameSizeList;
  std::vector<std::pair<int,int> > hydroPhobicList;
  std::vector<std::pair<int,int> > polarList;
  std::vector<std::pair<int,int> > posChargeList;
  std::vector<std::pair<int,int> > negChargeList;
  std::vector<std::pair<int,int> > coiledCoilList;
  randomMol rmg;
  std::vector<std::vector<double> > distSetsSecs;
  std::vector<double> distSets;
  double kapvallink,kapvalbeta,kapvalalpha,tauvallink,tauvalbeta,tauvalalpha,alvallink,alvalbeta,alvalalpha,maxDistChange;
  std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,std::pair<double,double> > > contactPairList;
  std::vector<int> unchangedSections;
};

#endif
