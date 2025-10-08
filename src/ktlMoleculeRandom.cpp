/* Carbonara Version: 0.2.0 (ktlMolecule – memory-friendly) */

#include "ktlMoleculeRandom.h" // (your original include)
#include "point.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>

// --- ctor, params ------------------------------------------------------------

ktlMolecule::ktlMolecule() {
  kapvallink  = 0.36932;
  kapvalbeta  = 0.343782;
  kapvalalpha = 0.380928;
  tauvallink  = 0.13439;
  tauvalbeta  = 0.123222;
  tauvalalpha = 0.145831;
  alvallink   = 4.22128;
  alvalbeta   = 4.13944;
  alvalalpha  = 4.2329;
}

void ktlMolecule::setParams(double& rminIn, double& rmaxIn, double& lminIn) {
  rmg.setParams(rminIn, rmaxIn, lminIn);
}

// --- simple getters / counters ----------------------------------------------

int ktlMolecule::getSubsecSize(int sec) const {
  // sec is 1-based; chainList is 0-based
  return chainList[sec - 1].second - chainList[sec - 1].first + 1;
}

int ktlMolecule::getUnitNo(int index) const {
  return static_cast<int>(coords[index].size());
}

point ktlMolecule::getCoordinate(int mindex, int subindex) const {
  if (subindex == -1) {
    return coords[mindex - 1][coords[mindex - 1].size() - 1];
  } else {
    return coords[mindex][subindex];
  }
}

std::string ktlMolecule::getType(int index) const {
  return nameSizeList[index].first;
}

std::string ktlMolecule::getType(int chainNo, int index) const {
  int fullIndex = chainList[chainNo - 1].first;
  return nameSizeList[fullIndex + index].first;
}

double ktlMolecule::getCurvatureJoined(int index) const {
  const int sz = static_cast<int>(coords[index].size());
  if (sz <= 2) return kapvallink;
  if (3 < sz && sz < 7) return kapvalbeta;
  return kapvalalpha;
}

double ktlMolecule::getTorsionJoined(int index) const {
  const int sz = static_cast<int>(coords[index].size());
  if (sz <= 2) return tauvallink;
  if (3 < sz && sz < 7) return tauvalbeta;
  return tauvalalpha;
}

double ktlMolecule::getAlbeadJoined(int index) const {
  const int sz = static_cast<int>(coords[index].size());
  if (sz <= 2) return alvallink;
  if (3 < sz && sz < 7) return alvalbeta;
  return alvalalpha;
}

int ktlMolecule::molSize() const {
  int sz = 0;
  for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
    if (nameSizeList[i].first == "Helix" || nameSizeList[i].first == "Strand") {
      sz += 1;
    } else {
      sz += static_cast<int>(coords[i].size());
    }
  }
  return sz;
}

int ktlMolecule::getNoAminos() const {
  int sz = 0;
  for (const auto& sec : coords) sz += static_cast<int>(sec.size());
  return sz;
}

double ktlMolecule::maxNeighbourDistSec(int sec) const {
  // sec is 1-based
  double dmax = 0.0;
  const int first = chainList[sec - 1].first;
  const int last  = chainList[sec - 1].second;

  for (int i = first; i <= last; ++i) {
    for (int j = 0; j < static_cast<int>(coords[i].size()); ++j) {
      double d = 0.0;
      if (j == static_cast<int>(coords[i].size()) - 1 && i < last) {
        d = coords[i][j].eDist(coords[i + 1][0]);
      } else if (j == static_cast<int>(coords[i].size()) - 1 && i == last) {
        d = 0.0;
      } else {
        d = coords[i][j].eDist(coords[i][j + 1]);
      }
      if (d > dmax) dmax = d;
    }
  }
  return dmax;
}

std::pair<double,double> ktlMolecule::getMaxPossibleLength() const {
  double maxL = 0.0;
  double totLength = 0.0;
  double dst = 0.0;
  for (int i = 0; i < static_cast<int>(coords.size()) - 1; ++i) {
    for (int j = i + 1; j < static_cast<int>(coords.size()); ++j) {
      dst = coords[i][0].eDist(coords[j][0]);
      if (j == i + 1) totLength += dst;
      if (dst > maxL) maxL = dst;
    }
  }
  maxL = maxL * 2.0;
  return {totLength, maxL};
}

// --- copy-returning sub-section helpers (compat) -----------------------------

std::vector<std::vector<point>> ktlMolecule::getSubsecCoordinates(int sec) const {
  const int first  = chainList[sec].first;
  const int second = chainList[sec].second;
  return std::vector<std::vector<point>>(coords.begin() + first, coords.begin() + second + 1);
}

std::vector<std::pair<std::string,int>> ktlMolecule::getNameSizeListOfSection(int sec) const {
  const int first  = chainList[sec].first;
  const int second = chainList[sec].second;
  return std::vector<std::pair<std::string,int>>(nameSizeList.begin() + first, nameSizeList.begin() + second + 1);
}

// --- IO (unchanged except minor cleanups) ------------------------------------

void ktlMolecule::readInSequence(const char* filename, double& rmin, double& rmax, double& lmin) {
  std::ifstream myfile(filename);
  std::string line;

  int noChains = 0;
  if (myfile.is_open()) {
    std::getline(myfile, line);
    std::stringstream ss(line);
    ss >> noChains;
    std::cout << "number of chains " << noChains << "\n";

    distSetsSecs.resize(noChains);
    for (int ch = 1; ch <= noChains; ++ch) {
      std::pair<int,int> p;
      p.first = static_cast<int>(nameSizeList.size());

      std::string empty, sequence, predictions;
      std::getline(myfile, empty);
      std::getline(myfile, sequence);
      std::getline(myfile, empty);
      std::getline(myfile, predictions);

      int n = 0;
      std::string type, prevType;
      std::vector<std::string> aminoType;

      for (size_t i = 0; i < predictions.size(); ++i) {
        if (i == 0) {
          if (predictions.compare(i,1,"-") == 0)       type = "loop";
          else if (predictions.compare(i,1,"E") == 0 ||
                   predictions.compare(i,1,"S") == 0)   type = "Strand";
          else                                          type = "Helix";
          n = 1;
          prevType = type;
          aminoType.push_back(sequence.substr(i,1));
        } else {
          if (predictions.compare(i,1,"-") == 0)       type = "loop";
          else if (predictions.compare(i,1,"E") == 0 ||
                   predictions.compare(i,1,"S") == 0)   type = "Strand";
          else                                          type = "Helix";

          if (prevType == type) {
            ++n;
            aminoType.push_back(sequence.substr(i,1));
          } else {
            nameSizeList.emplace_back(prevType, n);
            distChanges.push_back(0.0);
            prevType = type;
            aminoList.push_back(aminoType);
            aminoType.clear();
            aminoType.push_back(sequence.substr(i,1));
            n = 1;
          }
        }
      }

      nameSizeList.emplace_back(prevType, n);
      distChanges.push_back(0.0);
      aminoList.push_back(aminoType);
      aminoType.clear();

      p.second = static_cast<int>(nameSizeList.size()) - 1;
      chainList.push_back(p);
    }
    rmg.setParams(rmin, rmax, lmin);
  } else {
    std::cout << "no sequence provided for chain number header\n";
  }
}

void ktlMolecule::readInCoordinates(const char* filename) {
  std::ifstream coordFile(filename);
  std::string line;

  if (coordFile.is_open()) {
    std::vector<point> secondarySec;
    for (int iv = 0; iv < static_cast<int>(nameSizeList.size()); ++iv) {
      bool waitToFill = true;
      const int noInSection = nameSizeList[iv].second;
      int currNoMols = 0;

      while ((coordFile.eof() == false) && waitToFill == true) {
        std::getline(coordFile, line);
        if (line.find("chain") == std::string::npos) {
          if (line.size() > 1) {
            point p(line);
            currNoMols++;
            secondarySec.push_back(p);
            if (currNoMols == noInSection) {
              coords.push_back(secondarySec);
              secondarySec.clear();
              waitToFill = false;
            }
          }
        }
      }
    }
  } else {
    std::cout << "there is no valid coordinate file supplied\n";
  }
}

void ktlMolecule::writeMoleculeToFile(const char* filename) {
  std::ofstream ofile(filename);
  if (ofile.is_open()) {
    for (int k = 0; k < static_cast<int>(chainList.size()); ++k) {
      for (int i = chainList[k].first; i <= chainList[k].second; ++i) {
        for (int j = 0; j < static_cast<int>(coords[i].size()); ++j) {
          ofile << coords[i][j].getX() << " "
                << coords[i][j].getY() << " "
                << coords[i][j].getZ() << "\n";
        }
        ofile << "\n";
      }
      ofile << "End chain " << k + 1 << "\n";
    }
  } else {
    ofile << "cannot write molecule to file";
  }
}

// --- phys/chem lists (unchanged) --------------------------------------------

void ktlMolecule::getHydrophobicResidues() {
  for (int i = 0; i < static_cast<int>(aminoList.size()); ++i) {
    for (int j = 0; j < static_cast<int>(aminoList[i].size()); ++j) {
      const auto& a = aminoList[i][j];
      if (a=="A"||a=="I"||a=="L"||a=="M"||a=="F"||a=="V"||a=="P"||a=="G") {
        hydroPhobicList.emplace_back(i,j);
      }
    }
  }
}

void ktlMolecule::getCoiledCoilResidues() {
  for (int i = 0; i < static_cast<int>(aminoList.size()); ++i) {
    for (int j = 0; j < static_cast<int>(aminoList[i].size()); ++j) {
      const auto& a = aminoList[i][j];
      if (a=="L"||a=="I"||a=="G") {
        if (nameSizeList[i].first == "Helix") {
          coiledCoilList.emplace_back(i,j);
        }
      }
    }
  }
}

void ktlMolecule::getPolarResidues() {
  for (int i = 0; i < static_cast<int>(aminoList.size()); ++i) {
    for (int j = 0; j < static_cast<int>(aminoList[i].size()); ++j) {
      const auto& a = aminoList[i][j];
      if (a=="Q"||a=="N"||a=="H"||a=="S"||a=="T"||a=="Y"||a=="C") {
        polarList.emplace_back(i,j);
      }
    }
  }
}

void ktlMolecule::getPositiveResidues() {
  for (int i = 0; i < static_cast<int>(aminoList.size()); ++i) {
    for (int j = 0; j < static_cast<int>(aminoList[i].size()); ++j) {
      const auto& a = aminoList[i][j];
      if (a=="K"||a=="Y"||a=="R") {
        posChargeList.emplace_back(i,j);
      }
    }
  }
}

void ktlMolecule::getNegativeResidues() {
  for (int i = 0; i < static_cast<int>(aminoList.size()); ++i) {
    for (int j = 0; j < static_cast<int>(aminoList[i].size()); ++j) {
      const auto& a = aminoList[i][j];
      if (a=="D"||a=="E") {
        negChargeList.emplace_back(i,j);
      }
    }
  }
}

// --- distances / overlaps ----------------------------------------------------

std::vector<int> ktlMolecule::checkOverlap(std::vector<std::vector<point>>& cdsIN) {
  std::vector<int> overlappedSecs;
  for (int i = 0; i < static_cast<int>(cdsIN.size()) - 1; ++i) {
    bool triggered = false;
    for (int j = i + 1; j < static_cast<int>(cdsIN.size()); ++j) {
      for (int k = 0; k < static_cast<int>(cdsIN[i].size()); ++k) {
        for (int l = 0; l < static_cast<int>(cdsIN[j].size()); ++l) {
          double dist = cdsIN[i][k].eDist(cdsIN[j][l]);
          if (k == static_cast<int>(cdsIN[i].size()) - 1 && l == 0 && j == i + 1) dist = 10.0;
          if (dist < 5.0) triggered = true;
        }
      }
    }
    if (triggered) overlappedSecs.push_back(i);
  }
  return overlappedSecs;
}

std::vector<double> ktlMolecule::checkOverlapWithRad(double& wRad, int sec) {
  std::vector<double> overlappedSecs;
  distSetsSecs[sec - 1].clear();

  const int first = chainList[sec - 1].first;
  const int last  = chainList[sec - 1].second;

  for (int i = first; i < last; ++i) {
    bool triggered = false;
    for (int j = i + 1; j <= last; ++j) {
      for (int k = 0; k < static_cast<int>(coords[i].size()); ++k) {
        for (int l = 0; l < static_cast<int>(coords[j].size()); ++l) {
          double dist = coords[i][k].eDist(coords[j][l]);
          distSetsSecs[sec - 1].push_back(dist);

          if (k == static_cast<int>(coords[i].size()) - 1 && l == 0 && j == i + 1) dist = 10.0;
          if (dist < wRad) {
            triggered = true;
            overlappedSecs.push_back(dist);
          }
        }
      }
    }
    (void)triggered; // kept for parity; you used to ignore it
  }
  return overlappedSecs;
}

std::vector<double> ktlMolecule::checkOverlapWithRad(double& wRad) {
  std::vector<double> overlappedSecs;
  distSets.clear();
  minimumintraMolecularDistances.clear();
  minimumintraMolecularDistances.resize(chainList.size(), std::vector<double>(chainList.size(), 10000.0));
  minimumintraMolecularDistanceMean = 0;
  minimumintraMoleculaeDistancePerChain.clear();

  int chainIndex1 = 0, chainIndex2 = 0;

  for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
    if (i > chainList[chainIndex1].second) chainIndex1++;

    chainIndex2 = chainIndex1;
    for (int j = i; j < static_cast<int>(coords.size()); ++j) {
      if (j > chainList[chainIndex2].second) chainIndex2++;

      for (int k = 0; k < static_cast<int>(coords[i].size()); ++k) {
        for (int l = 0; l < static_cast<int>(coords[j].size()); ++l) {
          double dist = coords[i][k].eDist(coords[j][l]);
          if (i == j) {
            if (l > k) distSets.push_back(dist);
          } else {
            distSets.push_back(dist);
          }

          if (chainIndex1 != chainIndex2) {
            auto& a = minimumintraMolecularDistances[chainIndex1][chainIndex2];
            auto& b = minimumintraMolecularDistances[chainIndex2][chainIndex1];
            if (dist < a) { a = dist; b = dist; }
          }

          if ((k == static_cast<int>(coords[i].size()) - 1 && l == 0 && j == i + 1) || j == i) dist = 10.0;
          if (dist < wRad) overlappedSecs.push_back(dist);
        }
      }
    }
  }

  for (int i = 0; i < static_cast<int>(minimumintraMolecularDistances.size()); ++i) {
    auto ptr = std::min_element(minimumintraMolecularDistances[i].begin(), minimumintraMolecularDistances[i].end());
    const double min = *ptr;
    if (min > minimumintraMolecularDistanceMean && min < 9999.9) {
      minimumintraMolecularDistanceMean = min;
    }
    minimumintraMoleculaeDistancePerChain.push_back(min);
  }
  return overlappedSecs;
}

// --- misc geometry -----------------------------------------------------------

double ktlMolecule::compareDistances(std::vector<std::vector<point>>& coords2) const {
  double distDifTot = 0.0;
  double n = 0.0;
  for (int i = 0; i < static_cast<int>(coords.size()) - 1; ++i) {
    for (int j = i + 1; j < static_cast<int>(coords.size()); ++j) {
      for (int k = 0; k < static_cast<int>(coords[i].size()); ++k) {
        for (int l = 0; l < static_cast<int>(coords[j].size()); ++l) {
          const double dist  = coords[i][k].eDist(coords[j][l]);
          const double dist2 = coords2[i][k].eDist(coords2[j][l]);
          distDifTot += std::abs(dist2 - dist);
          n += 1.0;
        }
      }
    }
  }
  return (n > 0.0) ? (distDifTot / n) : 0.0;
}

bool ktlMolecule::checkCalphas(std::vector<std::vector<point>>& coordsIn) const {
  for (int i = 0; i < static_cast<int>(coordsIn.size()); ++i) {
    for (int j = 0; j < static_cast<int>(coordsIn[i].size()); ++j) {
      double dist = 0.0;
      if ((j == static_cast<int>(coordsIn[i].size()) - 1) && (i < static_cast<int>(coordsIn.size()) - 1)) {
        dist = coordsIn[i][j].eDist(coordsIn[i + 1][0]);
      } else if ((j == static_cast<int>(coordsIn[i].size()) - 1) && (i == static_cast<int>(coordsIn.size()) - 1)) {
        dist = 3.7;
      } else {
        dist = coordsIn[i][j].eDist(coordsIn[i][j + 1]);
      }
      if (dist > 4.0 || dist < 3.3) return true;
    }
  }
  return false;
}

bool ktlMolecule::checkCalphas() const {
  for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
    for (int j = 0; j < static_cast<int>(coords[i].size()); ++j) {
      double dist = 0.0;
      if ((j == static_cast<int>(coords[i].size()) - 1) && (i < static_cast<int>(coords.size()) - 1)) {
        dist = coords[i][j].eDist(coords[i + 1][0]);
      } else if ((j == static_cast<int>(coords[i].size()) - 1) && (i == static_cast<int>(coords.size()) - 1)) {
        dist = 3.7;
      } else {
        dist = coords[i][j].eDist(coords[i][j + 1]);
      }
      if (dist > 4.0 || dist < 3.3) return true;
    }
  }
  return false;
}

bool ktlMolecule::checkCalphas(int secNo) const {
  const int sec = secNo - 1;
  const int first = chainList[sec].first;
  const int last  = chainList[sec].second;

  for (int i = first; i <= last; ++i) {
    for (int j = 0; j < static_cast<int>(coords[i].size()); ++j) {
      double dist = 0.0;
      if ((j == static_cast<int>(coords[i].size()) - 1) && (i < last)) {
        dist = coords[i][j].eDist(coords[i + 1][0]);
      } else if ((j == static_cast<int>(coords[i].size()) - 1) && (i == last)) {
        dist = 3.7;
      } else {
        dist = coords[i][j].eDist(coords[i][j + 1]);
      }
      if (dist > 4.6 || dist < 2.6) return true;
    }
  }
  return false;
}

bool ktlMolecule::checkCalphas(int secNo, const ktlMolecule& ktl) const {
  const int sec = secNo - 1;

  const auto& coordsOld = ktl.getCoordinates();

  const int first = chainList[sec].first;
  const int last  = chainList[sec].second;

  for (int i = 0; i <= (last - first); ++i) {
    const auto& curr = coords[first + i];
    const auto& old  = coordsOld[first + i];
    for (int j = 0; j < static_cast<int>(curr.size()); ++j) {
      double dist    = 0.0;
      double distOld = 0.0;

      if ((j == static_cast<int>(curr.size()) - 1) && (i < (last - first))) {
        dist    = curr[j].eDist(coords[first + i + 1][0]);
        distOld = old[j].eDist(coordsOld[first + i + 1][0]);
      } else if ((j == static_cast<int>(curr.size()) - 1) && (i == (last - first))) {
        dist = distOld = 3.7;
      } else {
        dist    = curr[j].eDist(curr[j + 1]);
        distOld = old[j].eDist(old[j + 1]);
      }

      if (std::abs(dist - distOld) > 0.9) return true;
    }
  }
  return false;
}

// --- transforms --------------------------------------------------------------

point ktlMolecule::getCentreOfMass(std::vector<std::vector<point>>& cdSet) const {
  point mean(0.0, 0.0, 0.0);
  int noMols = 0;
  for (auto& sec : cdSet) {
    for (auto& p : sec) { mean = mean + p; noMols++; }
  }
  return mean * (1.0 / double(noMols));
}

point ktlMolecule::getCentreOfMassRange(int first, int last) const {
  point mean(0.0, 0.0, 0.0);
  int noMols = 0;
  for (int i = first; i <= last; ++i) {
    for (const auto& p : coords[i]) { mean = mean + p; noMols++; }
  }
  return mean * (1.0 / double(noMols));
}

void ktlMolecule::rotation3D(point& p, point& centre, point& k, double& cosangle, double& sinangle) {
  point out = p - centre;
  point term1 = out * cosangle;
  point term2 = k.cross(out); term2 = term2 * sinangle;
  double dotProd = k.dotprod(out);
  point term3 = k * dotProd; term3 = term3 * (1.0 - cosangle);
  out = term1 + term2 + term3;
  out = out + centre;
  p = out;
}

void ktlMolecule::rotateSection(std::vector<std::vector<point>>& section, point& centre, point& k, double& angle, point& transVec) {
  double c = std::cos(angle), s = std::sin(angle);
  for (auto& sec : section) {
    for (auto& p : sec) {
      rotation3D(p, centre, k, c, s);
      p = p + transVec;
    }
  }
}

void ktlMolecule::changeMoleculeMultiRotate(double& angle, point& k, int secIn, point& transVec) {
  const int sec = secIn - 1;
  const int first = chainList[sec].first;
  const int last  = chainList[sec].second;

  point com = getCentreOfMassRange(first, last);
  double c = std::cos(angle), s = std::sin(angle);
  for (int i = first; i <= last; ++i) {
    for (auto& p : coords[i]) {
      rotation3D(p, com, k, c, s);
      p = p + transVec;
    }
  }
}

void ktlMolecule::changeMoleculeSingle(int& index, std::vector<std::vector<point>>& cdsIn, std::vector<std::pair<std::string,int>>& nameSizeSubList) {
  point sp(0.0, 0.0, 0.0);
  bool suc = true;
  maxDistChange = rmg.reshapeMol(nameSizeSubList, cdsIn, index, sp, suc);
}

void ktlMolecule::changeMoleculeSingleMulti(int index, int secIn) {
  const int sec = secIn - 1;
  point sp(0.0, 0.0, 0.0);

  const int first = chainList[sec].first;
  const int last  = chainList[sec].second;

  // build sub-lists into scratch (reused capacity)
  scratchSubcoords.clear();
  scratchSubcoords.insert(scratchSubcoords.end(), coords.begin() + first, coords.begin() + last + 1);

  std::vector<std::pair<std::string,int>> subns(nameSizeList.begin() + first, nameSizeList.begin() + last + 1);

  bool suc = true;

  // pre/post size check
  auto countPoints = [](const std::vector<std::vector<point>>& v) {
    size_t s = 0; for (const auto& x : v) s += x.size(); return s;
  };
  const size_t preSize = countPoints(scratchSubcoords);

  if (index < static_cast<int>(scratchSubcoords.size()) - 2) {
    if (subns[index].first == "loop" && subns[index+1].first == "Helix" && subns[index+2].first == "Helix") {
      int noHelices = 0; int k = index + 1;
      while (k < static_cast<int>(scratchSubcoords.size()) && subns[k].first == "Helix") { noHelices++; k++; }
      maxDistChange = rmg.reshapeMolLoopThenHelixSet(subns, scratchSubcoords, index, noHelices, sp, suc);
    } else if (subns[index].first == "Helix" && subns[index+1].first == "Helix" && subns[index+2].first == "Helix") {
      int noHelices = 0; int k = index + 1;
      while (k < static_cast<int>(scratchSubcoords.size()) && subns[k].first == "Helix") { noHelices++; k++; }
      maxDistChange = rmg.reshapeMolHelixSet(subns, scratchSubcoords, index, noHelices, sp, suc);
    } else {
      maxDistChange = rmg.reshapeMol(subns, scratchSubcoords, index, sp, suc);
    }
  } else {
    maxDistChange = rmg.reshapeMol(subns, scratchSubcoords, index, sp, suc);
  }

  const size_t postSize = countPoints(scratchSubcoords);
  if (preSize == postSize) {
    std::copy_n(scratchSubcoords.begin(), (last - first + 1), &coords[first]);
  }
}

// --- replicate ---------------------------------------------------------------

void ktlMolecule::replicateMolecule(int& noReplications) {
  std::random_device rdev{};
  std::default_random_engine gen{rdev()};
  std::uniform_real_distribution<> rotAng(0.0, 0.5);
  std::uniform_real_distribution<> theAng(0.0, 3.14159265359);
  std::uniform_real_distribution<> phiAng(0.0, 6.28318530718);

  point com = getCentreOfMass(const_cast<std::vector<std::vector<point>>&>(coords)); // reuse existing function
  double maxRad = 0.0;
  for (const auto& sec : coords)
    for (const auto& p : sec)
      maxRad = std::max(maxRad, com.eDist(p));

  std::vector<std::vector<point>> replicatedCoords = coords;
  std::vector<std::pair<std::string,int>> replicatedNameSizeList = nameSizeList;
  std::vector<std::pair<int,int>> chainListCopy = chainList;
  std::vector<std::vector<std::string>> aminoListCopy = aminoList;

  replicatedCoords.reserve(coords.size() * (noReplications + 1));
  replicatedNameSizeList.reserve(nameSizeList.size() * (noReplications + 1));
  chainListCopy.reserve(chainList.size() * (noReplications + 1));
  aminoListCopy.reserve(aminoList.size() * (noReplications + 1));

  for (int i = 1; i <= noReplications; ++i) {
    std::vector<std::vector<point>> newPts = coords;
    double angle = rotAng(gen);
    double theta = theAng(gen);
    double phi   = phiAng(gen);
    point k(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    theta = theAng(gen); phi = phiAng(gen);
    point translate(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    translate.scalarMult(2.2 * maxRad);
    rotateSection(newPts, com, k, angle, translate);

    replicatedCoords.insert(replicatedCoords.end(), newPts.begin(), newPts.end());

    std::vector<std::pair<std::string,int>> copyNameList(nameSizeList.begin(), nameSizeList.end());
    replicatedNameSizeList.insert(replicatedNameSizeList.end(), copyNameList.begin(), copyNameList.end());

    std::pair<int,int> prInts;
    prInts.first  = chainListCopy.back().second + 1;
    prInts.second = chainListCopy.back().second + (chainList.back().second - chainList.back().first + 1);
    chainListCopy.push_back(prInts);

    aminoListCopy.insert(aminoListCopy.end(), aminoList.begin(), aminoList.end());
  }

  aminoList  = std::move(aminoListCopy);
  nameSizeList = std::move(replicatedNameSizeList);
  chainList  = std::move(chainListCopy);
  coords     = std::move(replicatedCoords);
}

// --- kap/tau -----------------------------------------------------------------

std::vector<std::pair<double,double>> ktlMolecule::getKapTauVals() {
  std::vector<std::pair<double,double>> ktllst;
  for (int k = 0; k < static_cast<int>(chainList.size()); ++k) {
    for (int i = chainList[k].first; i <= chainList[k].second; ++i) {
      for (int j = 0; j < static_cast<int>(coords[i].size()); ++j) {
        if (i != chainList[k].second) {
          if (j < static_cast<int>(coords[i].size()) - 3) {
            ktllst.push_back(rmg.kapTau(coords[i][j], coords[i][j+1], coords[i][j+2], coords[i][j+3]));
          } else if (j == static_cast<int>(coords[i].size()) - 3) {
            ktllst.push_back(rmg.kapTau(coords[i][j], coords[i][j+1], coords[i][j+2], coords[i+1][0]));
          } else if (j == static_cast<int>(coords[i].size()) - 2) {
            ktllst.push_back(rmg.kapTau(coords[i][j], coords[i][j+1], coords[i+1][0], coords[i+1][1]));
          } else if (j == static_cast<int>(coords[i].size()) - 1) {
            ktllst.push_back(rmg.kapTau(coords[i][j], coords[i][0], coords[i+1][1], coords[i+1][2]));
          }
        } else {
          if (j < static_cast<int>(coords[i].size()) - 3) {
            ktllst.push_back(rmg.kapTau(coords[i][j], coords[i][j+1], coords[i][j+2], coords[i][j+3]));
          }
        }
      }
    }
  }
  return ktllst;
}

// --- potentials / distances --------------------------------------------------

double ktlMolecule::lennardJones(double& idealDist, double& currDist, int noConPred, double& weightCoeff) {
  const double dif = (idealDist - currDist) / idealDist;
  return weightCoeff * dif * dif;
}

std::vector<double> ktlMolecule::solMolDists(std::vector<std::vector<point>>& pts1) {
  std::vector<double> distSet;
  for (auto& sec : pts1) {
    for (auto& p1 : sec) {
      for (auto& csec : coords) {
        for (auto& p2 : csec) {
          distSet.push_back(p1.eDist(p2));
        }
      }
    }
  }
  return distSet;
}

// --- constraints -------------------------------------------------------------

void ktlMolecule::loadContactPredictions(const char* contactloc) {
  std::ifstream cpfile(contactloc);
  if (cpfile.is_open()) {
    std::string output;
    while (!cpfile.eof()) {
      int ind1, ind2, ind3, ind4; double distance, percentage;
      std::getline(cpfile, output);
      std::stringstream ss(output);
      if (output.length() > 0) {
        ss >> ind1; ss.ignore(); ss >> ind2;
        std::pair<int,int> pr1{ind1, ind2};
        ss.ignore(); ss >> ind3; ss.ignore(); ss >> ind4; ss.ignore();
        std::pair<int,int> pr2{ind3, ind4};
        ss >> distance; ss.ignore(); ss >> percentage;
        std::pair<double,double> pr3{distance, percentage};
        contactPairList.emplace_back(pr1, pr2, pr3);
      }
    }
  }
}

double ktlMolecule::getLennardJonesContact() {
  double ljval = 0.0;
  for (const auto& tp : contactPairList) {
    const auto& pr1 = std::get<0>(tp);
    const auto& pr2 = std::get<1>(tp);
    const auto& targ = std::get<2>(tp);

    point cd1 = coords[pr1.first][pr1.second];
    point cd2 = coords[pr2.first][pr2.second];
    const double prWiseDist = cd1.eDist(cd2);

    const double distFrac = (prWiseDist - targ.first) / (targ.first);
    const double distFracWeighted = distFrac / (targ.second);
    ljval += 0.0001 * distFracWeighted * distFracWeighted * distFracWeighted * distFracWeighted;
  }
  return ljval;
}

void ktlMolecule::loadFixedSections(const char* fixedsecloc) {
  std::ifstream fsfile(fixedsecloc);
  if (fsfile.is_open()) {
    std::string output;
    while (!fsfile.eof()) {
      std::getline(fsfile, output);
      std::stringstream ss(output);
      int section = 0; ss >> section;
      if (section != 0) unchangedSections.push_back(section);
    }
  }
}

// --- beta-sheet reward -------------------------------------------------------

double ktlMolecule::getBetaSheetProximityReward() {
  std::vector<std::vector<point>> betaSheets;
  for (int sec = 0; sec < static_cast<int>(coords.size()); ++sec) {
    if (nameSizeList[sec].first == "Strand") betaSheets.push_back(coords[sec]);
  }

  numBetaSheets = static_cast<double>(betaSheets.size());
  if (betaSheets.size() < 2) return 0.0;

  double totalReward = 0.0;
  for (int i = 0; i < static_cast<int>(betaSheets.size()) - 1; ++i) {
    for (int j = i + 1; j < static_cast<int>(betaSheets.size()); ++j) {
      double minDistance = 10000.0;
      for (auto& p1 : betaSheets[i]) {
        for (auto& p2 : betaSheets[j]) {
          double d = p1.eDist(p2);
          if (d < minDistance) minDistance = d;
        }
      }

      const bool isAligned = areSheetsPairallelOrAntiparallel(betaSheets[i], betaSheets[j]);
      if (isAligned) {
        double reward = 0.0;
        const double Dopt = 5.0;
        const double Dmax = 10.0;
        const double Dmin = 3.8;

        if (minDistance <= Dmin) {
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
  point dir1 = sheet1.back() - sheet1.front();
  point dir2 = sheet2.back() - sheet2.front();
  dir1.normalise(); dir2.normalise();
  double dotProduct = std::abs(dir1.dotprod(dir2));
  return dotProduct > 0.8;
}

// --- memory hygiene ----------------------------------------------------------

void ktlMolecule::shrink_temporaries() {
  distSets.shrink_to_fit();
  for (auto& v : distSetsSecs) v.shrink_to_fit();
  scratchSubcoords.shrink_to_fit();
}
