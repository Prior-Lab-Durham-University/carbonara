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

// Loads in (if availible) contanct constraints to referenced mol class
void readFixedDistancesConstraints(const char* argv[], std::vector<ktlMolecule>& mol);

// Loads in mixture file into parameters
void readPermissibleMixtures(const char* argv[], ModelParameters& params);

// find number of subsections in each molecule - e.g. for a monomer/dimer mixture noSections[0]=1,noSections[1]=2.
std::vector<int> findNumberSections(std::vector<ktlMolecule>& mol);

// construct the historical states set
std::vector<moleculeFitAndState> makeHistoricalStateSet(moleculeFitAndState molState, ModelParameters& params);


// Increase kmax range logic - changes made in place with reference
void increaseKmax(std::pair<double,double>& scatterFit, std::vector<moleculeFitAndState>& molFitAndStateSet,
                  experimentalData& ed,  ModelParameters& params, Logger& logger);


// make referenced change to the newMol, return bool for c-alpha check
bool modifyMolecule(ktlMolecule& newMol, ktlMolecule& existingMol, int indexCh, int section);

// update mol, molState, overallFit and write result to log
void updateAndLog(int& improvementIndex, std::vector<ktlMolecule>& mol, ktlMolecule& newMol,
                  moleculeFitAndState& molState, moleculeFitAndState& newMolState,
                  std::pair<double,double>& overallFit, std::pair<double,double>& newOverallFit,
                  Logger& logger, int l, int k, experimentalData& ed, ModelParameters& params);


// Construct molecule file name - accounting for sub structure
std::string constructMoleculeName(const std::string& basePath, const std::string& prefix, const std::string& extension,
                                  const int& submol, const int& improvementIndex, const std::string& body);

// Construct scattering file name
std::string constructScatterName(const std::string& basePath, const std::string& prefix, const std::string& extension,
                                 const int& improvementIndex, const std::string& body);

// writes all sub-molecules of vector[ktlMolecule object] to files
std::string write_molecules(const std::string& basePath, const int& improvementIndex, std::vector<ktlMolecule>& mol, const std::string& body);

// writes simulated scattering to file and returns the scatterName (scatterName needed for logging!)
std::string write_scatter(const std::string& basePath, const int& improvementIndex, moleculeFitAndState& molFit,
                          experimentalData& ed, double kmin, double kmaxCurr,std::vector<std::vector<double> > & mixtureList, const std::string& body = "main");

// Do we accept or reject a new fitting - method implemented in here
bool checkTransition(double &chiSqVal, double &chiSqCurr,double &uniformProb,int index,int &maxSteps);

void sortVec(std::vector<moleculeFitAndState> &mfs);

void tokenize(std::string &str, const char delim, std::vector<std::string> &out);

double getHydrophobicPackingPenalty(double &packValue);


class RandomGenerator {
private:

    std::random_device rdev;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distTran;
    std::uniform_real_distribution<double> rotAng;
    std::uniform_real_distribution<double> theAng;
    std::uniform_real_distribution<double> phiAng;
    std::uniform_real_distribution<double> distributionR;  // Uniform real distribution between 0.0 and 1.0

public:

    RandomGenerator();
    double getDistTran();
    double getRotAng();
    double getTheAng();
    double getPhiAng();
    double getDistributionR();
    int getChangeIndexProbability(int& k, ModelParameters& params);
};

#endif
