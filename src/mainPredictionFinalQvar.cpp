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

// combination of all structures = moleculeStructures
// each structure is 

// note: this version showing funky behaviour with getFit()'s - not always consistent when recalled!

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

void molProbe(std::vector<ktlMolecule>& mol, ModelParameters params, experimentalData ed, std::string output) {

  moleculeFitAndState copy_molState(mol, params);
  // find the overall fit of the updated copy
  std::pair<double,double> copy_molOverallFit_wol = copy_molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  std::cout << output << copy_molOverallFit_wol.first << "\n";
}

int main(int argc, const char* argv[]) {

  /* initialise the log file */
  Logger logger(argv[16]);

  /* Set up model parameters */
  ModelParameters params = loadParameters(argv);

  /* Determine initial model: Two options no initial prediction, we must generate a structure
   or some initial structure provided. Actually we need a half-half option */

  /* Initialise the molecule(s) vector */
  std::vector<ktlMolecule> moleculeStructures;
  readInStructures(argv, moleculeStructures, params);

  /* Determine which sections are being altered */
  std::vector<std::vector<int>> vary_sec_list_list;
  determineVaryingSections(argv, vary_sec_list_list);

  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, moleculeStructures);

  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);

  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  /* Random generator */
  RandomGenerator rng;

  /* initialise the state of mol vector */
  moleculeFitAndState molState(moleculeStructures, params);

  int improvementIndex = 0;
  // If we resume from previous run - argv[3] restart True/False
  if ((strcmp(argv[3], "True") == 0)) {
    improvementIndex = std::atoi(argv[17]);
  }

  std::pair<double, double> overallFit;
  if (params.affineTrans == true) {
    overallFit = molState.getOverallFitForceConnection(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  } else {
    overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);
  }

  logger.logMetadata(argv[16], params);

  std::string scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, "initial");
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, moleculeStructures, "initial");

  // log starting point
  logger.logEntry(0, 0, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial);

  logger.consoleInitial(overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), molState.getDistanceConstraints());
  
  /* Main algorithm */

  // numberOfChainsInEachStructure vector tells us how many chains are in each structure
  // e.g. for a monomer/dimer mixture numberOfChainsInEachStructure[0]=1, numberOfChainsInEachStructure[1]=2.
  std::vector<int> numberOfChainsInEachStructure = findNumberSections(moleculeStructures);
  
  /* initialise the set of historical states - currently basic, but used to save previous fit stages */
  std::vector<moleculeFitAndState> molStateSet = makeHistoricalStateSet(molState, params);

  // loop number
  int fitStep = 0;
  // This is a monster while loop - strap in chaps
  while (fitStep < params.noScatterFitSteps) {
    // Increasing the kmax if we have a good enough fit, consider a little more of the experimental data!
    if (overallFit.second < 0.0005 || (params.improvementIndexTest > std::round(params.noScatterFitSteps / 5) && overallFit.second < 0.0007)) {
      increaseKmax(overallFit, molStateSet, ed, params, logger);
    }

    params.improvementIndexTest = params.improvementIndexTest + 1;

    // pick a 'random' molState from the historical molStateSet
    // to become update function
    int historicFitIndex = rng.getChangeIndexProbability(fitStep, params);
    molState = molStateSet[historicFitIndex];
    moleculeStructures = molState.getMolecule();
    overallFit = molState.getFit();

    for (int structureIndex = 0; structureIndex < moleculeStructures.size(); structureIndex++) {
      int netIndex = 0;
      
      // loop over the sections of the given molecule (i.e. if its a monomer this loop is tivial, but not for a multimer
      // another monster looooooop
      for (int chainNumber = 1; chainNumber <= numberOfChainsInEachStructure[structureIndex]; chainNumber++) {
        
        // Selected transformation option?
        if (params.affineTrans == true) {
          
          ktlMolecule molCopyR = moleculeStructures[structureIndex];

          double angle = rng.getRotAng();
          double theta = rng.getTheAng();
          double phi = rng.getPhiAng();
          point kv(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));

          double xtran = rng.getDistTran();
          double ytran = rng.getDistTran();
          double ztran = rng.getDistTran();
          point tranVec(xtran, ytran, ztran);

          molCopyR.changeMoleculeMultiRotate(angle, kv, chainNumber, tranVec);
          bool cacaDist = molCopyR.checkCalphas(chainNumber);

          if (cacaDist == false) {

            // calculate the new fit for this
            moleculeFitAndState newMolState = molState;
            std::pair<double, double> newOverallFit = newMolState.getOverallFit(ed, params.mixtureList, params.helRatList, molCopyR, params.kmin, params.kmaxCurr, structureIndex);
            double uProb = rng.getDistributionR();

            if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {

              improvementIndex++;
              updateAndLog(improvementIndex, moleculeStructures, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);

              logger.consoleChange("fitImprove", params);
            }
          }
        } // rotate/translate section ends

        // net index tells us how far we are through the whole molecule
        if (chainNumber > 1) {
          
          netIndex = netIndex + moleculeStructures[structureIndex].getSubsecSize(chainNumber - 1);
        }

        bool doAll = false;

        // Now loop over the secondary structures of the given unit or section
        for (int secondarySectionIndex = 0; secondarySectionIndex < moleculeStructures[structureIndex].getSubsecSize(chainNumber) - 1; secondarySectionIndex++) {
          
          int totalIndex = netIndex + secondarySectionIndex;

          // in this if statement we check which secondary sections are being changed
          if ((doAll == true) || (std::find(vary_sec_list_list[structureIndex].begin(), vary_sec_list_list[structureIndex].end(), totalIndex) != vary_sec_list_list[structureIndex].end())) {
            int indexCh = totalIndex - netIndex;
            ktlMolecule newMol = moleculeStructures[structureIndex];
            bool cacaDist = modifyMolecule(newMol, moleculeStructures[structureIndex], indexCh, chainNumber);

            if (cacaDist == false) {

              moleculeFitAndState newmolState = molState;

              // calculate the fitting of changed molecule
              std::pair<double, double> newOverallFit = newmolState.getOverallFit(ed, params.mixtureList, params.helRatList, newMol, params.kmin, params.kmaxCurr, structureIndex);

              double uProb = rng.getDistributionR();

              if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {

                // Success! Add to the update index
                improvementIndex++;
                updateAndLog(improvementIndex, moleculeStructures, newMol, molState, newmolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
                logger.consoleChange("fitImprove", params);

              }
            }

          } // totalIndex an allowed varying section?

        } // structureIndex
      } // chainNumber
    } // structureIndex

    // Assign the new 'improved' molecule state to the historical tracker
    molStateSet[historicFitIndex] = molState;
    molStateSet[historicFitIndex].updateMolecule(moleculeStructures);
    sortVec(molStateSet);

    // Print out to terminal window
    logger.consoleFitAttempt(fitStep, improvementIndex, params, overallFit.first, overallFit.second);

    fitStep++;
  }

  improvementIndex++;

  std::string moleculeNameEnd = write_molecules(argv[12], improvementIndex, moleculeStructures, "end");
  std::string scatterNameEnd = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr, "end");

  std::cout << "\n best overall mol name: " << moleculeNameEnd << "\n";
  std::cout << " overallFitBest fit: " << overallFit.first << "\n";

  logger.logEntry(improvementIndex, fitStep, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd);

} // end of main
