

/* Carbonara Version: 0.2.0 */

#include "ktlMoleculeRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "moleculeFitAndState.h"
#include <cstring>
#include <chrono>
#include <tuple>
#include <utility>   // std::swap

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
  argv[ 8] kmin
  argv[ 9] kmax
  argv[10] kmax_start
  argv[11] Max number of fitting steps
  argv[12] prediction file - mol[i] in the fitting folder
  argv[13] scattering output file
  argv[14] mixture list file, a list of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
  argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
  argv[16] log file location
  argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
  argv[18] is true if we want to apply affine rotations, false if not.
  argv[19] is true if the user wants to use error weightings
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- */


int main(int argc, const char* argv[]) {

  /* initialise the log file */
  Logger logger(argv[16]);

  /* Set up model parameters */
  ModelParameters params = loadParameters(argv);

  /* Initialise the molecule(s) vector */
  std::vector<ktlMolecule> initialMolecules;
  readInStructures(argv, initialMolecules, params);

  /* Determine which sections are being altered */
  std::vector<std::vector<int>> vary_sec_list_list;
  determineVaryingSections(argv, vary_sec_list_list);

  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, initialMolecules);

  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);

  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  /* Random generator */
  RandomGenerator rng;

  /* initialise the state of mol vector (takes ownership of initialMolecules by copy) */
  moleculeFitAndState molState(initialMolecules, params);

  int improvementIndex = 0;
  if ((strcmp(argv[3], "True") == 0)) {
    improvementIndex = std::atoi(argv[17]);
  }

  std::pair<double, double> overallFit;
  if (params.affineTrans == true) {
    if ((strcmp(argv[19], "True") == 0)) {
      overallFit = molState.getOverallFitForceConnection_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    } else {
      overallFit = molState.getOverallFitForceConnection(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }
  } else {
    if ((strcmp(argv[19], "True") == 0)) {
      overallFit = molState.getOverallFit_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    } else {
      overallFit = molState.getOverallFit(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }
  }

  logger.logMetadata(argv[16], params);

  // Initial outputs (write from molState directly — no extra big vector copies)
  std::string scatterNameInitial;
  if (strcmp(argv[19], "True") == 0) {
    scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, params.mixtureList, "initial");
  } else {
    scatterNameInitial = write_scatter_ChiSq(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, params.mixtureList, "initial");
  }
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, molState, "initial");

  logger.logEntry(0, 0, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial, molState.C2);
  logger.consoleInitial(overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), molState.getDistanceConstraints());

  // numberOfChainsInEachStructure
  const std::vector<int> numberOfChainsInEachStructure = findNumberSections(molState.getMolecule());

  /* historical states */
  std::vector<moleculeFitAndState> molStateSet = makeHistoricalStateSet(molState, params);

  int fitStep = 0;

  while (fitStep < params.noScatterFitSteps) {

    if (overallFit.second < 0.0002) {
      increaseKmax(overallFit, molStateSet, ed, params, logger);
    }

    params.improvementIndexTest += 1;

    // pick a historical state without copying the whole object
    int historicFitIndex = rng.getChangeIndexProbability(fitStep, params);
    std::swap(molState, molStateSet[historicFitIndex]);

    // Alias the current working molecule vector inside molState (no copy).
    auto& molecules = const_cast<std::vector<ktlMolecule>&>(molState.getMolecule());

    for (int structureIndex = 0; structureIndex < static_cast<int>(molecules.size()); ++structureIndex) {
      int netIndex = 0;

      for (int chainNumber = 1; chainNumber <= numberOfChainsInEachStructure[structureIndex]; ++chainNumber) {

        // Selected transformation option?
        if (params.affineTrans == true) {
          ktlMolecule molCopyR = molecules[structureIndex]; // copy only this sub-structure

          double angle = rng.getRotAng();
          double theta = rng.getTheAng();
          double phi   = rng.getPhiAng();
          point kv(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));

          point tranVec(rng.getDistTran(), rng.getDistTran(), rng.getDistTran());

          molCopyR.changeMoleculeMultiRotate(angle, kv, chainNumber, tranVec);
          bool cacaDist = molCopyR.checkCalphas(chainNumber, molecules[structureIndex]);
          if (cacaDist == false) {

            // evaluate by copying the state (cheapest safe way with current API)
            moleculeFitAndState newMolState = molState;
            std::pair<double, double> newOverallFit;
            if ((strcmp(argv[19], "True") == 0)) {
              newOverallFit = newMolState.getOverallFitForceConnection_ChiSq(ed, params.mixtureList, molCopyR, params.kmin, params.kmaxCurr, structureIndex);
            } else {
              newOverallFit = newMolState.getOverallFitForceConnection(ed, params.mixtureList, molCopyR, params.kmin, params.kmaxCurr, structureIndex);
            }

            double uProb = rng.getDistributionR();
            if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {
              improvementIndex++;
              if ((strcmp(argv[19], "True") == 0)) {
                updateAndLog_ChiSq(improvementIndex, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
              } else {
                updateAndLog(improvementIndex, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
              }
              logger.consoleChange("fitImprove", params);
              if ((strcmp(argv[19], "True") == 0)) {
                molState.updateScatteringFit_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
              } else {
                molState.updateScatteringFit(ed, params.mixtureList, params.kmin, params.kmaxCurr);
              }
            }
          }
        } // affine

        if (chainNumber > 1) {
          netIndex += molecules[structureIndex].getSubsecSize(chainNumber - 1);
        }

        const bool doAll = false;

        for (int secondarySectionIndex = 0; secondarySectionIndex < molecules[structureIndex].getSubsecSize(chainNumber) - 1; ++secondarySectionIndex) {

          int totalIndex = netIndex + secondarySectionIndex;

          if ((doAll == true) ||
              (std::find(vary_sec_list_list[structureIndex].begin(),
                         vary_sec_list_list[structureIndex].end(),
                         totalIndex) != vary_sec_list_list[structureIndex].end())) {

            int indexCh = totalIndex - netIndex;
            ktlMolecule newMol = molecules[structureIndex]; // copy only this sub-structure
            bool cacaDist = modifyMolecule(newMol, molecules[structureIndex], indexCh, chainNumber);
            cacaDist = newMol.checkCalphas(chainNumber, molecules[structureIndex]);

            if (cacaDist == false) {
              moleculeFitAndState newmolState = molState;
              std::pair<double, double> newOverallFit;
              if ((strcmp(argv[19], "True") == 0)) {
                newOverallFit = newmolState.getOverallFit_ChiSq(ed, params.mixtureList, newMol, params.kmin, params.kmaxCurr, structureIndex);
              } else {
                newOverallFit = newmolState.getOverallFit(ed, params.mixtureList, newMol, params.kmin, params.kmaxCurr, structureIndex);
              }

              double uProb = rng.getDistributionR();
              if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {
                improvementIndex++;
                if ((strcmp(argv[19], "True") == 0)) {
                  updateAndLog_ChiSq(improvementIndex, newMol, molState, newmolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
                } else {
                  updateAndLog(improvementIndex, newMol, molState, newmolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
                }
                logger.consoleChange("fitImprove", params);
                if ((strcmp(argv[19], "True") == 0)) {
                  molState.updateScatteringFit_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
                } else {
                  molState.updateScatteringFit(ed, params.mixtureList, params.kmin, params.kmaxCurr);
                }
              }
            }
          }
        } // secondarySectionIndex
      }   // chainNumber
    }     // structureIndex

    // Put the (possibly updated) state back without copying
    std::swap(molState, molStateSet[historicFitIndex]);
    sortVec(molStateSet);

    logger.consoleFitAttempt(fitStep, improvementIndex, params, overallFit.first, overallFit.second);
    fitStep++;
  }

  improvementIndex++;

  // final outputs (from current best working state in molStateSet[0] after sort)
  // The best is at index 0 after sort; if desired, you can std::swap into molState first.
  std::swap(molState, molStateSet[0]);

  std::string moleculeNameEnd = write_molecules(argv[12], improvementIndex, molState, "end");
  std::string scatterNameEnd;
  if ((strcmp(argv[19], "True") == 0)) {
    scatterNameEnd = write_scatter_ChiSq(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, params.mixtureList, "end");
  } else {
    scatterNameEnd = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, params.mixtureList, "end");
  }
  std::cout << "\n best overall mol name: " << moleculeNameEnd << "\n";
  std::cout << " overallFitBest fit: " << overallFit.first << "\n";

  logger.logEntry(improvementIndex, fitStep, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd, molState.C2);

} // end of main
