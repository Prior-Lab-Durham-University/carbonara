
#include "helpers.h"
#include <sstream>
#include <fstream>
#include <algorithm>

// Loads in structural data to referenced mol class
void readInStructures(const char* argv[], std::vector<ktlMolecule>& mol, ModelParameters& params) {

    int noStructures = std::atoi(argv[6]);

    if (noStructures <= 0) {
        std::cerr << "Invalid number of structures specified." << std::endl;
        return;
    }

    if (strcmp(argv[3], "True") == 0) {
        // Restart from a fit obtained on a previous run
        if (argv[15] == nullptr) {
            std::cerr << "No file string provided for restart." << std::endl;
            return;
        }

        std::stringstream filestring(argv[15]);
        std::string segment;
        std::vector<std::string> seglist;

        while (std::getline(filestring, segment, '+')) {
            seglist.push_back(segment);
        }

        if (seglist.empty()) {
            std::cerr << "No segments found for restart." << std::endl;
            return;
        }

        for (size_t i = 0; i < seglist.size(); i++) {
            ktlMolecule molTmp;
            std::string sequenceLoc = std::string(argv[2]) + "fingerPrint" + std::to_string(i + 1) + ".dat";
            molTmp.readInSequence(sequenceLoc.c_str(), params.rmin, params.rmax, params.lmin);
            molTmp.readInCoordinates(seglist[i].c_str());
            molTmp.getHydrophobicResidues();
            mol.push_back(std::move(molTmp));
        }

    } else {
        // Fresh start
        mol.reserve(noStructures);
        for (int i = 0; i < noStructures; i++) {
            ktlMolecule molTmp;
            std::string sequenceLoc   = std::string(argv[2]) + "fingerPrint" + std::to_string(i + 1) + ".dat";
            std::string coordinateLoc = std::string(argv[2]) + "coordinates"  + std::to_string(i + 1) + ".dat";
            molTmp.readInSequence(sequenceLoc.c_str(), params.rmin, params.rmax, params.lmin);
            molTmp.readInCoordinates(coordinateLoc.c_str());
            molTmp.getHydrophobicResidues();
            mol.push_back(std::move(molTmp));
        }
    }
}

// Loads in the allowed varying sections
void determineVaryingSections(const char* argv[], std::vector<std::vector<int>>& vary_sec_list_list) {

    int noStructures = std::atoi(argv[6]);
    vary_sec_list_list.clear();
    vary_sec_list_list.reserve(noStructures);

    for (int i = 0; i < noStructures; i++) {

        std::ifstream vary_sec_file;
        std::vector<int> vary_sec_list;
        std::string vary_sec_loc = std::string(argv[2]) + "varyingSectionSecondary" + std::to_string(i + 1) + ".dat";
        vary_sec_file.open(vary_sec_loc.c_str());
        std::string line;
        int index;

        if (vary_sec_file.is_open()) {

            while (std::getline(vary_sec_file, line)) {
                if (line.empty()) continue;
                std::stringstream ss(line);
                while (ss >> index) vary_sec_list.push_back(index);
                // ignore trailing non-numbers silently
            }

            vary_sec_list_list.push_back(std::move(vary_sec_list));

        } else {
            std::cerr << "Failed to open varying section file: " << vary_sec_loc << std::endl;
        }
        vary_sec_file.close();
    }
}

// Loads in (if available) contact constraints
void readFixedDistancesConstraints(const char* argv[], std::vector<ktlMolecule>& mol) {

    int noStructures = std::atoi(argv[6]);

    if (strcmp(argv[4], "True") == 0) {
      for (int i = 0; i < noStructures; i++) {
        std::string contactPredictions = std::string(argv[2]) + "fixedDistanceConstraints" + std::to_string(i + 1) + ".dat";
        std::cout << contactPredictions << "\n";
        mol[i].loadContactPredictions(contactPredictions.c_str());
      }
    }
}

// Loads in mixture file into parameters
void readPermissibleMixtures(const char* argv[], ModelParameters& params) {

    std::string filePath = argv[14];

    std::ifstream permissibleMixtureFile(filePath);
    if (!permissibleMixtureFile) {
        std::cerr << "Failed to open mixture file: " << filePath << std::endl;
        return;
    }

    std::vector<std::vector<double>> mixtureList;
    std::string line;

    while (std::getline(permissibleMixtureFile, line)) {
        if (line.empty()) continue;

        std::vector<double> mixtureSet;
        std::stringstream lineStream(line);
        double value;

        while (lineStream >> value) mixtureSet.push_back(value);

        if (!lineStream.eof()) {
            std::cerr << "Warning: Encountered non-numeric data in: " << line << std::endl;
            continue;
        }

        if (!mixtureSet.empty()) {
            mixtureList.push_back(std::move(mixtureSet));
        } else {
            std::cerr << "Warning: Empty numeric line found in file: " << filePath << std::endl;
        }
    }

    if (mixtureList.empty()) {
      std::cerr << "No valid mixtures were read from the file: " << filePath << std::endl;
    }

    params.mixtureList = std::move(mixtureList);
}

// number of chains for each structure
std::vector<int> findNumberSections(const std::vector<ktlMolecule>& mol) {
    std::vector<int> noSections;
    noSections.reserve(mol.size());
    for (const auto& m : mol) {
        noSections.push_back(m.noChains());
    }
    return noSections;
}

// Add the original molState N times to the historical set
std::vector<moleculeFitAndState> makeHistoricalStateSet(const moleculeFitAndState& molState, ModelParameters& params){
    std::vector<moleculeFitAndState> molStateSet;
    molStateSet.reserve(params.noHistoricalFits);
    for (int i = 0; i < params.noHistoricalFits; i++) {
        molStateSet.push_back(molState); // one-time copies at start
    }
    return molStateSet;
}

void increaseKmax(std::pair<double,double>& scatterFit, std::vector<moleculeFitAndState>& molFitAndStateSet,
                  experimentalData& ed,  ModelParameters& params, Logger& logger) {

    params.kmaxCurr = std::min(params.kmax, params.kmaxCurr + 0.01);
    logger.consoleChange("krangeIncrease", params);
    params.improvementIndexTest = 0;

    // recompute from the first historical state (already in the vector)
    scatterFit = molFitAndStateSet[0].getOverallFit(ed, params.mixtureList, params.kmin, params.kmaxCurr);
}

bool modifyMolecule(ktlMolecule& newMol, const ktlMolecule& existingMol, int indexCh, int section) {
    (void)existingMol; // used only for post-check below
    newMol.changeMoleculeSingleMulti(indexCh, section);
    return newMol.checkCalphas(section, const_cast<ktlMolecule&>(existingMol));
}


void updateAndLog(int& improvementIndex, const ktlMolecule& newMol,
                  moleculeFitAndState& molState, moleculeFitAndState& newMolState,
                  std::pair<double,double>& overallFit, const std::pair<double,double>& newOverallFit,
                  Logger& logger, int l, int k, experimentalData& ed, ModelParameters& params) {

    // 1) Adopt the accepted state (penalties, caches, etc.)
    molState = newMolState;

    molState = newMolState;
    molState.replaceMoleculeAt(l, newMol);


    // 3) Commit the accepted fit
    overallFit = newOverallFit;

    // 4) Output
    std::string moleculeNameMain = write_molecules(params.basePath, improvementIndex, molState, "default");
    std::string scatterNameMain  = write_scatter(params.basePath, improvementIndex, molState, ed,
                                                 params.kmin, params.kmaxCurr, params.mixtureList);

    logger.logEntry(improvementIndex, k, overallFit.first, molState.getWrithePenalty(),
                    molState.getOverlapPenalty(), molState.getDistanceConstraints(),
                    params.kmaxCurr, scatterNameMain, moleculeNameMain, molState.C2);
}

void updateAndLog_ChiSq(int& improvementIndex, const ktlMolecule& newMol,
                  moleculeFitAndState& molState, moleculeFitAndState& newMolState,
                  std::pair<double,double>& overallFit, const std::pair<double,double>& newOverallFit,
                  Logger& logger, int l, int k, experimentalData& ed, ModelParameters& params) {

    molState = newMolState;

    molState = newMolState;
    molState.replaceMoleculeAt(l, newMol);
    overallFit = newOverallFit;

    overallFit = newOverallFit;

    std::string moleculeNameMain = write_molecules(params.basePath, improvementIndex, molState, "default");
    std::string scatterNameMain  = write_scatter_ChiSq(params.basePath, improvementIndex, molState, ed,
                                                       params.kmin, params.kmaxCurr, params.mixtureList);

    logger.logEntry(improvementIndex, k, overallFit.first, molState.getWrithePenalty(),
                    molState.getOverlapPenalty(), molState.getDistanceConstraints(),
                    params.kmaxCurr, scatterNameMain, moleculeNameMain, molState.C2);
}

std::string constructMoleculeName(const std::string& basePath, const std::string& /*prefix*/, const std::string& extension,
                                  int submol, int improvementIndex, const std::string& body) {

    std::stringstream ss;
    if (body == "initial")      { ss << basePath << "_sub_" << submol << "_initial_xyz" << extension; }
    else if (body == "end")     { ss << basePath << "_sub_" << submol << "_end_xyz"     << extension; }
    else                        { ss << basePath << "_sub_" << submol << "_step_" << improvementIndex << "_xyz" << extension; }
    return ss.str();
}

std::string constructScatterName(const std::string& basePath, const std::string& /*prefix*/, const std::string& extension,
                                 int improvementIndex, const std::string& body) {

    std::stringstream ss;
    if (body == "initial")      { ss << basePath << "_initial_scatter" << extension; }
    else if (body == "end")     { ss << basePath << "_end_scatter"     << extension; }
    else                        { ss << basePath << "_step_" << improvementIndex << "_scatter" << extension; }
    return ss.str();
}

// writes all sub-molecules of current state to files (no extra large copies)
std::string write_molecules(const std::string& basePath, int improvementIndex, moleculeFitAndState& molState, const std::string& body) {
    auto& mol = const_cast<std::vector<ktlMolecule>&>(molState.getMolecule());
    std::string moleculeName;
    for (int i = 0; i < static_cast<int>(mol.size()); i++) {
        moleculeName = constructMoleculeName(basePath, "xyz", ".dat", i, improvementIndex, body);
        mol[i].writeMoleculeToFile(moleculeName.c_str());
    }
    return moleculeName;
}

std::string write_scatter(const std::string& basePath, int improvementIndex, moleculeFitAndState& molFit,
                          experimentalData& ed, double kmin, double kmaxCurr, std::vector<std::vector<double>>& mixtureList, const std::string& body) {

    std::string scatterName = constructScatterName(basePath, "scatter", ".dat", improvementIndex, body);
    molFit.writeScatteringToFile(ed, mixtureList, scatterName.c_str());
    return scatterName;
}

std::string write_scatter_ChiSq(const std::string& basePath, int improvementIndex, moleculeFitAndState& molFit,
                          experimentalData& ed, double kmin, double kmaxCurr, std::vector<std::vector<double>>& mixtureList, const std::string& body) {

    std::string scatterName = constructScatterName(basePath, "scatter", ".dat", improvementIndex, body);
    molFit.writeScatteringToFile_ChiSq(ed, mixtureList, scatterName.c_str());
    return scatterName;
}

bool checkTransition(double &chiSqVal, double &chiSqCurr, double &uniformProb, int /*index*/, int &/*maxSteps*/) {
  (void)uniformProb;
  return (chiSqVal < chiSqCurr);
}

void sortVec(std::vector<moleculeFitAndState> &mfs){
  std::sort(mfs.begin(), mfs.end(),[](const moleculeFitAndState &x, const moleculeFitAndState &y) {
    return x.currFit < y.currFit;
  });
}

void tokenize(std::string &str, const char delim, std::vector<std::string> &out) {
    std::stringstream ss(str);
    std::string s;
    while (std::getline(ss, s, delim)) out.push_back(s);
}

double getHydrophobicPackingPenalty(double &packValue){
  return 0.00001*std::exp(3.0*(packValue-1.6));
}

// RNG -------------------------------------------------------------------------

RandomGenerator::RandomGenerator()
    : generator(rdev()),
      distTran(-10.0, 10.0),
      rotAng(0.0, 2.0),
      theAng(0.0, 3.14159265359),
      phiAng(0.0, 6.28318530718),
      distributionR(0.0, 1.0) {}

double RandomGenerator::getDistTran() { return distTran(generator); }
double RandomGenerator::getRotAng()   { return rotAng(generator); }
double RandomGenerator::getTheAng()   { return theAng(generator); }
double RandomGenerator::getPhiAng()   { return phiAng(generator); }
double RandomGenerator::getDistributionR() { return distributionR(generator); }

int RandomGenerator::getChangeIndexProbability(int k, ModelParameters& params) {
    double p = 0.7 - 0.6 * (static_cast<double>(k) / std::max(1, params.noScatterFitSteps));
    if (p < 0.0) p = 0.0; if (p > 1.0) p = 1.0;
    std::binomial_distribution<> changeIndexProbability(std::max(0, params.noHistoricalFits - 1), p);
    return changeIndexProbability(generator);
}
