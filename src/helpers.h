#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <random>

#include "ktlMoleculeRandom.h"
#include "moleculeFitAndState.h"
#include "experimentalData.h"
#include "parameters.h"
#include "Logger.h"

// Loads in structural data to referenced mol class
void readInStructures(const char* argv[], std::vector<ktlMolecule>& mol, ModelParameters& params);

// Loads in the allowed varying sections to referenced mol class
void determineVaryingSections(const char* argv[], std::vector<std::vector<int>>& vary_sec_list_list);

// Loads in (if available) contact constraints to referenced mol class
void readFixedDistancesConstraints(const char* argv[], std::vector<ktlMolecule>& mol);

// Loads in mixture file into parameters
void readPermissibleMixtures(const char* argv[], ModelParameters& params);

// find number of subsections in each molecule
std::vector<int> findNumberSections(const std::vector<ktlMolecule>& mol);

// construct the historical states set
std::vector<moleculeFitAndState> makeHistoricalStateSet(const moleculeFitAndState& molState, ModelParameters& params);

// Increase kmax range logic - changes made in place with reference
void increaseKmax(std::pair<double,double>& scatterFit, std::vector<moleculeFitAndState>& molFitAndStateSet,
                  experimentalData& ed,  ModelParameters& params, Logger& logger);

// make referenced change to the newMol, return bool for c-alpha check
bool modifyMolecule(ktlMolecule& newMol, const ktlMolecule& existingMol, int indexCh, int section);

// update molState (internal molecules), log, and persist files
void updateAndLog(int& improvementIndex, const ktlMolecule& newMol,
                  moleculeFitAndState& molState, moleculeFitAndState& newMolState,
                  std::pair<double,double>& overallFit, const std::pair<double,double>& newOverallFit,
                  Logger& logger, int l, int k, experimentalData& ed, ModelParameters& params);

void updateAndLog_ChiSq(int& improvementIndex, const ktlMolecule& newMol,
                  moleculeFitAndState& molState, moleculeFitAndState& newMolState,
                  std::pair<double,double>& overallFit, const std::pair<double,double>& newOverallFit,
                  Logger& logger, int l, int k, experimentalData& ed, ModelParameters& params);

// Construct molecule file name
std::string constructMoleculeName(const std::string& basePath, const std::string& prefix, const std::string& extension,
                                  int submol, int improvementIndex, const std::string& body);

// Construct scattering file name
std::string constructScatterName(const std::string& basePath, const std::string& prefix, const std::string& extension,
                                 int improvementIndex, const std::string& body);

// writes all sub-molecules of molState to files (no extra copy)
std::string write_molecules(const std::string& basePath, int improvementIndex, moleculeFitAndState& molState, const std::string& body);

// writes simulated scattering to file and returns the scatterName
std::string write_scatter(const std::string& basePath, int improvementIndex, moleculeFitAndState& molFit,
                          experimentalData& ed, double kmin, double kmaxCurr, std::vector<std::vector<double>>& mixtureList, const std::string& body = "main");

std::string write_scatter_ChiSq(const std::string& basePath, int improvementIndex, moleculeFitAndState& molFit,
                          experimentalData& ed, double kmin, double kmaxCurr, std::vector<std::vector<double>>& mixtureList, const std::string& body = "main");

// accept / reject
bool checkTransition(double &chiSqVal, double &chiSqCurr, double &uniformProb, int index, int &maxSteps);

void sortVec(std::vector<moleculeFitAndState> &mfs);

void tokenize(std::string &str, const char delim, std::vector<std::string> &out);

double getHydrophobicPackingPenalty(double &packValue);

// RNG
class RandomGenerator {
private:
    std::random_device rdev;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distTran;
    std::uniform_real_distribution<double> rotAng;
    std::uniform_real_distribution<double> theAng;
    std::uniform_real_distribution<double> phiAng;
    std::uniform_real_distribution<double> distributionR;

public:
    RandomGenerator();
    double getDistTran();
    double getRotAng();
    double getTheAng();
    double getPhiAng();
    double getDistributionR();
    int getChangeIndexProbability(int k, ModelParameters& params); // by value
};

#endif
