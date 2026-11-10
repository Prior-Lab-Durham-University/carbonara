#ifndef PARAMS
#define PARAMS

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

/**
 * \struct ModelParameters
 * \brief A structure representing the parameters used by Carbonara.
 *
 * Some parameters are set to default values, while others are user-definable.
 */
struct ModelParameters {

    double lmin = 4.0;
    double rmin = 3.7;
    double rmax = 3.9;
    double closestApproachDist = 3.9;
    double Rin = 6.0;
    double Rout = 7.0;
    double RShell = 5.5;
    int ntrivs = 6;
    int solventsPerLink = 1;
    double q_lim = 0.15;

    // user definable
    double kmin = 0.0;
    double kmax = 0.0;
    double kmaxCurr = 0.0;
    int noScatterFitSteps = 5000;
    bool affineTrans=false;

    std::vector<double> helRatList;
    std::vector< std::vector<double> > mixtureList;
    std::string basePath;

    int improvementIndexTest = 0;
    int noHistoricalFits = 1;
    
};

/**
 * \brief Implementation of the loadParameters function.
 * 
 * This function takes an array of command-line arguments and uses them to fill a ModelParameters object.
 * 
 * \param argv Carbonara command-line arguments.
 * \return ModelParameters object.
 */
ModelParameters loadParameters(const char* argv[]);

#endif