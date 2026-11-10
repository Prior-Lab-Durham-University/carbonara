/* Carbonara Version: 0.2.0 */

#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "moleculeFitAndState.h"
#include <cstring>
#include <chrono>
#include <tuple>

#include "Logger.h"
#include "helpers.h"

using namespace std::chrono;

/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

  argv[ 1] scattering data file
  argv[ 2] sequence file location
  argv[ 3] restart tag (use to start from existing prediction)
  argv[ 4] paired distances file (can be empty)
  argv[ 5] fixed sections file (again can be empty)
  argv[ 6] number of structures
  argv[ 7] request to apply hydrophobic covering WITHIN monomers -- Currently not used
  argv[ 8] request to apply hydrophobic covering BETWEEN monomers -- currently not used
  argv[ 9] kmin
  argv[10] kmax
  argv[11] Max number of fitting steps
  argv[12] prediction file - mol[i] in the fitting folder
  argv[13] scattering output file
  argv[14] mixture list file, a list of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
  argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
  argv[16] log file location
  argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
  argv[18] is true if we want to apply affine rotations, false if not.


 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- */

double fitToScatteringMultiple(
    std::vector<std::vector<double>>& molDists,
    std::vector<std::vector<double>>& solDists,
    std::vector<std::vector<double>>& solMolDists,
    std::vector<int>& molSizes,
    std::vector<int>& solSizes,
    std::vector<double>& percentageCombinations) {

    /*****************************************
         Helper function for binning distances
    ******************************************/
    auto binDistances = [&](std::vector<std::vector<double>>& dists, std::vector<int>& sizes, std::vector<double>& distNo) {
        for (size_t j = 0; j < dists.size(); ++j) {
            size_t k = 0;
            for (size_t i = 0; i < distBins.size(); ++i) {  // distBins is not defined
                const auto& binPr = distBins[i];  // distBins is not defined
                size_t startK = k;
                while (dists[j][k] > binPr.first && dists[j][k] <= binPr.second) {
                    ++k;
                }
                size_t endK = k;
                if (i != 0) {
                    distNo[i] += percentageCombinations[j] * (2.0 * double(endK - startK));
                } else {
                    distNo[i] += percentageCombinations[j] * (2.0 * double(endK - startK) + double(sizes[j]));
                }
            }
        }
    };

    std::vector<double> molDistNo(distBins.size(), 0.0);  // distBins is not defined
    std::vector<double> solDistNo(distBins.size(), 0.0);  // distBins is not defined
    std::vector<double> solMolDistNo(distBins.size(), 0.0);  // distBins is not defined

    binDistances(molDists, molSizes, molDistNo);
    binDistances(solDists, solSizes, solDistNo);

    // Handle solMolDists separately because it does not use solSizes
    for (size_t j = 0; j < solMolDists.size(); ++j) {
        size_t k = 0;
        for (size_t i = 0; i < distBins.size(); ++i) {  // distBins is not defined
            const auto& binPr = distBins[i];  // distBins is not defined
            size_t startK = k;
            while (solMolDists[j][k] > binPr.first && solMolDists[j][k] <= binPr.second) {
                ++k;
            }
            size_t endK = k;
            solMolDistNo[i] += percentageCombinations[j] * 2.0 * double(endK - startK);
        }
    }

    /*****************************************
         Construct the scattering formula
    ******************************************/
    double bestPred = 1000.0;

    for (int hi = 0; hi < 5; ++hi) {
        double hydRat = 0.75 + 0.5 * hi / 4.0;
        std::vector<double> scatVals;
        
        for (size_t i = 0; i < scatPhases[0].size(); ++i) {  // scatPhases is not defined
            double scatk = 0.0;
            double aminoScat = gaussianFixed(qvals[i]);  // gaussianFixed and qvals are not defined
            double solScat = gaussianHydFixed(qvals[i]);  // gaussianHydFixed and qvals are not defined

            for (size_t j = 0; j < scatPhases.size(); ++j) {  // scatPhases is not defined
                scatk += scatPhases[j][i] * (aminoScat * aminoScat * molDistNo[j] +
                    hydRat * hydRat * solScat * solScat * solDistNo[j] +
                    hydRat * solScat * aminoScat * solMolDistNo[j]);
            }
            scatVals.push_back(scatk);
        }

        std::vector<double> logDiffs;
        double logDiffMean = 0.0;
        int countMean = 0;

        for (size_t i = 0; i < scatVals.size(); ++i) {
            double logScatDiff = std::log(scatVals[i]) - std::log(experimentalIntensity[i]);  // experimentalIntensity is not defined
            logDiffs.push_back(logScatDiff);

            if (qvals[i] < 0.1) {  // qvals is not defined
                logDiffMean += logScatDiff;
                ++countMean;
            }
        }

        logDiffMean /= double(countMean);
        double predTemp = 0.0;

        for (const auto& logDiff : logDiffs) {
            double scatDiff = logDiff - logDiffMean;
            predTemp += scatDiff * scatDiff;
        }

        if (predTemp < bestPred) {
            bestPred = predTemp;
        }
    }

    return bestPred / (scatPhases[0].size() - 1);  // scatPhases is not defined
}



int main(int argc, const char* argv[]) {

  /* initialise the log file */
  Logger logger(argv[16]);

  /* Set up model parameters */
  ModelParameters params = loadParameters(argv);

  /* Initialise the molecule(s) vector */
  std::vector<ktlMolecule> moleculeStructures;
  readInStructures(argv, moleculeStructures, params);

  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);

  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  /* initialise the state of mol vector */
  moleculeFitAndState molState(moleculeStructures, params);

  std::pair<double, double> overallFit;
  overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);

  logger.logMetadata(argv[16], params);

  int improvementIndex = 0;

  std::string scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, "initial");
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, moleculeStructures, "initial");

  // log starting point
  logger.logEntry(0, 0, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial);

  logger.consoleInitial(overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), molState.getDistanceConstraints());

} // end of main
