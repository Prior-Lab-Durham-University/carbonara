#ifndef KTL_MOL
#define KTL_MOL

#include "point.h"
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <dirent.h>
#include "polyHelix.h"
#include <random>
#include "randomMolGen.h"
#include <tuple>
#include <vector>
#include <string>
#include <utility>

class ktlMolecule {
public:
  // lifecycle
  ktlMolecule();
  ~ktlMolecule() = default;

  ktlMolecule(const ktlMolecule&) = default;
  ktlMolecule& operator=(const ktlMolecule&) = default;

  ktlMolecule(ktlMolecule&&) noexcept = default;
  ktlMolecule& operator=(ktlMolecule&&) noexcept = default;

  // params / setup
  void setParams(double& rminIn, double& rmaxIn, double& lminIn);

  // ---- ZERO-COPY GETTERS (const refs) ----
  const std::vector<int>&                               getUnitNos()        const { return noPts; }
  const std::vector<std::vector<point>>&                getTangents()       const { return tanlist; }
  const std::vector<std::vector<point>>&                getNormals()        const { return normlist; }
  const std::vector<std::vector<point>>&                getBinormals()      const { return binormlist; }
  const std::vector<std::vector<point>>&                getCoordinates()    const { return coords; }
  const std::vector<std::vector<std::string>>&          getAminoList()      const { return aminoList; }
  const std::vector<double>&                            getDistChanges()    const { return distChanges; }
  const std::vector<std::pair<std::string,int>>&        getNameSizeList()   const { return nameSizeList; }
  const std::vector<std::pair<int,int>>&                getChainList()      const { return chainList; }
  const std::vector<double>&                            getDistSet()        const { return distSets; }

  // single-section (zero-copy) view
  const std::vector<point>& getCoordinatesSection(int i) const { return coords[i]; }

  // NOTE: These return copies (API-compatible), prefer range-based alternatives below.
  std::vector<std::vector<point>>        getSubsecCoordinates(int sec) const;
  std::vector<std::pair<std::string,int>> getNameSizeListOfSection(int sec) const;

  // section/range helpers (no copies)
  int  getSubsecSize(int sec) const; // sec is 1-based chain index
  int  noSecSize()            const { return static_cast<int>(coords.size()); }
  int  noChains()             const { return static_cast<int>(chainList.size()); }
  int  molSize()              const;  // weighted count (loops fully, helix/strand as 1)
  int  getNoAminos()          const;

  // geometry getters (small-by-value)
  int         getUnitNo(int index) const;
  point       getTangent(int mindex, int subindex)   const { return tanlist[mindex][subindex]; }
  point       getNormal(int mindex, int subindex)    const { return normlist[mindex][subindex]; }
  point       getBinormal(int mindex, int subindex)  const { return binormlist[mindex][subindex]; }
  point       getCoordinate(int mindex, int subindex) const;
  std::string getType(int index)                     const;
  std::string getType(int chainNo, int index)        const;
  double      getCurvatureJoined(int index)          const;
  double      getTorsionJoined(int index)            const;
  double      getAlbeadJoined(int index)             const;
  double      maxNeighbourDistSec(int sec)           const;
  double      getMaxDistChange()                     const { return maxDistChange; }
  std::pair<double,double> getMaxPossibleLength()    const;

  // IO / build
  void readInSequence(const char* filename, double& rmin, double& rmax, double& lmin);
  void readInCoordinates(const char* filename);
  void writeMoleculeToFile(const char* filename);

  // phys/chem annotations
  void getHydrophobicResidues();
  void getCoiledCoilResidues();
  void getPositiveResidues();
  void getNegativeResidues();
  void getPolarResidues();

  std::vector<double> getHydrophobicDistance(std::vector<std::vector<point>>& solventList, double& maxSolDist);

  // scoring
  double coiledCoilPotential();
  double coiledCoilPotentialBetween(int secNo);
  double coiledCoilPotentialBetween();
  point  getCentreOfMass(std::vector<std::vector<point>>& cdSet) const;
  point  getCentreOfMassRange(int first, int last) const; // [first,last] inclusive in coords

  std::vector<int>    checkOverlap(std::vector<std::vector<point>>& cdsIN);
  std::vector<double> checkOverlapWithRad(double& wRad, int sec);
  std::vector<double> checkOverlapWithRad(double& wRad);

  double getMininumIntraMoelcularDistance() const { return minimumintraMolecularDistanceMean; }
  std::vector<double>              getMinimumintraMoleculaeDistancePerChain() const { return minimumintraMoleculaeDistancePerChain; }
  std::vector<std::vector<double>> getMinimumintraMolecularDistances() const { return minimumintraMolecularDistances; }

  double compareDistances(std::vector<std::vector<point>>& coords2) const;

  bool checkCalphas(std::vector<std::vector<point>>& coordsIn) const;
  bool checkCalphas() const;
  bool checkCalphas(int secNo) const;
  bool checkCalphas(int secNo, const ktlMolecule& ktl) const;

  // editing / transforms
  void changeMoleculeSingle(int& index, std::vector<std::vector<point>>& cdsIn, std::vector<std::pair<std::string,int>>& nameSizeSubList);
  int  getRandomMolecule();
  int  getRandomMoleculeReset(); // (unchanged impl elsewhere)
  void resetRandomMolecule();    // (unchanged impl elsewhere)

  void changeMoleculeSingleMulti(int index, int sec);
  void changeMoleculeMultiRotate(double& angle, point& k, int secIn, point& transVec);

  void replicateMolecule(int& noReplications);
  void rotation3D(point& p, point& centre, point& k, double& cosangle, double& sinangle);
  void rotateSection(std::vector<std::vector<point>>& section, point& centre, point& k, double& angle, point& transVec);

  // scattering helpers
  std::vector<std::pair<double,double>> getKapTauVals();

  // distances & potentials
  double getPairDistance(std::pair<int,int>& index1, std::pair<int,int>& index2); // (declared earlier)
  double lennardJones(double& idealDist, double& currDist, int noConPred, double& weightCoeff);
  std::vector<double> solMolDists(std::vector<std::vector<point>>& pts1);

  // constraints
  void   loadContactPredictions(const char* contactloc);
  double getLennardJonesContact();
  void   loadFixedSections(const char* fixedsecloc);

  // beta-sheet reward
  double getBetaSheetProximityReward();
  bool   areSheetsPairallelOrAntiparallel(std::vector<point>& sheet1, std::vector<point>& sheet2);
  double numBetaSheets = 0.0;

  // memory hygiene (optional)
  void shrink_temporaries();

private:
  // data
  std::vector<double> minimumintraMoleculaeDistancePerChain;
  std::vector<std::vector<double>> minimumintraMolecularDistances;
  double minimumintraMolecularDistanceMean = 0.0;

  std::vector<int> noPts;
  std::vector<std::vector<point>> coords;
  std::vector<std::vector<std::string>> aminoList;
  std::vector<std::vector<point>> tanlist;
  std::vector<std::vector<point>> normlist;
  std::vector<std::vector<point>> binormlist;

  polyHelix ph;
  std::vector<std::pair<int,int>> chainList; // per chain: [firstIdx, lastIdx] into coords/nameSizeList
  std::vector<double> distChanges;
  std::vector<std::pair<std::string,int>> nameSizeList;

  std::vector<std::pair<int,int>> hydroPhobicList;
  std::vector<std::pair<int,int>> polarList;
  std::vector<std::pair<int,int>> posChargeList;
  std::vector<std::pair<int,int>> negChargeList;
  std::vector<std::pair<int,int>> coiledCoilList;

  randomMol rmg;
  std::vector<std::vector<double>> distSetsSecs;
  std::vector<double> distSets;

  double kapvallink = 0.0, kapvalbeta = 0.0, kapvalalpha = 0.0;
  double tauvallink = 0.0, tauvalbeta = 0.0, tauvalalpha = 0.0;
  double alvallink  = 0.0, alvalbeta  = 0.0, alvalalpha  = 0.0;
  double maxDistChange = 0.0;

  std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,std::pair<double,double>>> contactPairList;
  std::vector<int> unchangedSections;

  // scratch to avoid per-call allocations when we must materialize subranges
  std::vector<std::vector<point>> scratchSubcoords;
};

#endif

