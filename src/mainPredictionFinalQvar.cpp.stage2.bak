/* Carbonara Version: 0.2.0 */

#include "ktlMoleculeRandom.h"
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
    if ((strcmp(argv[19], "True") == 0)) {
      overallFit = molState.getOverallFitForceConnection_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }else{
      overallFit = molState.getOverallFitForceConnection(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }
  } else {
     if ((strcmp(argv[19], "True") == 0)) {
      overallFit = molState.getOverallFit_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }else{
      overallFit = molState.getOverallFit(ed, params.mixtureList, params.kmin, params.kmaxCurr);
    }
  }
  logger.logMetadata(argv[16], params);
  std::string scatterNameInitial;
  if(strcmp(argv[19], "True") == 0) {
    scatterNameInitial = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr,params.mixtureList, "initial");
    }else{
    scatterNameInitial = write_scatter_ChiSq(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr,params.mixtureList, "initial");
    }
  std::string xyzNameInitial = write_molecules(argv[12], improvementIndex, moleculeStructures, "initial");
   // log starting point
  logger.logEntry(0, 0, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial, molState.C2);
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

//    if (overallFit.second < 0.0002 || (params.improvementIndexTest > std::round(params.noScatterFitSteps / 5) && overallFit.second < 0.0007)) {

      if (overallFit.second < 0.0002) {

      increaseKmax(overallFit, molStateSet, ed, params, logger);
    }

    params.improvementIndexTest = params.improvementIndexTest + 1;

    // pick a 'random' molState from the historical molStateSet
    // to become update function
    int historicFitIndex = rng.getChangeIndexProbability(fitStep, params);
    molState = molStateSet[historicFitIndex];
    moleculeStructures = molState.getMolecule();
    //overallFit = molState.getFit();

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
          bool cacaDist= molCopyR.checkCalphas(chainNumber,moleculeStructures[structureIndex]);
          if (cacaDist == false) {

            // calculate the new fit for this
            moleculeFitAndState newMolState = molState;
            std::pair<double, double> newOverallFit;
            if ((strcmp(argv[19], "True") == 0)) {
            newOverallFit = newMolState.getOverallFitForceConnection_ChiSq(ed, params.mixtureList, molCopyR, params.kmin, params.kmaxCurr, structureIndex);
            }else{
              newOverallFit = newMolState.getOverallFitForceConnection(ed, params.mixtureList, molCopyR, params.kmin, params.kmaxCurr, structureIndex);
            }
            double uProb = rng.getDistributionR();
            if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {
              improvementIndex++;
               if ((strcmp(argv[19], "True") == 0)) {
                updateAndLog_ChiSq(improvementIndex, moleculeStructures, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
              }else{
                updateAndLog(improvementIndex, moleculeStructures, molCopyR, molState, newMolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
              }
              logger.consoleChange("fitImprove", params);
              if ((strcmp(argv[19], "True") == 0)) {
	         molState.updateScatteringFit_ChiSq(ed,params.mixtureList,params.kmin, params.kmaxCurr);
               }else{
                 molState.updateScatteringFit(ed,params.mixtureList,params.kmin, params.kmaxCurr);
               }    
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
	    cacaDist= newMol.checkCalphas(chainNumber,moleculeStructures[structureIndex]);
            if (cacaDist == false) {

              moleculeFitAndState newmolState = molState;

              // calculate the fitting of changed molecule
              std::pair<double, double> newOverallFit;
              if ((strcmp(argv[19], "True") == 0)) {
               newOverallFit= newmolState.getOverallFit_ChiSq(ed, params.mixtureList, newMol, params.kmin, params.kmaxCurr, structureIndex);
               }else{
		newOverallFit= newmolState.getOverallFit(ed, params.mixtureList, newMol, params.kmin, params.kmaxCurr, structureIndex);
               }
	      //std::cout<<"improve ever ? "<<indexCh<<" "<<newOverallFit.second<<" "<<newOverallFit.first<<" "<<overallFit.first<<"\n";
              double uProb = rng.getDistributionR();
              if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps)) {

                // Success! Add to the update index
                improvementIndex++;
                if ((strcmp(argv[19], "True") == 0)) {
                  updateAndLog_ChiSq(improvementIndex, moleculeStructures, newMol, molState, newmolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
                  }else{
                    updateAndLog(improvementIndex, moleculeStructures, newMol, molState, newmolState, overallFit, newOverallFit, logger, structureIndex, fitStep, ed, params);
                  }
                logger.consoleChange("fitImprove", params);
                if ((strcmp(argv[19], "True") == 0)) {
		   molState.updateScatteringFit_ChiSq(ed, params.mixtureList, params.kmin, params.kmaxCurr);     
                }else{
                   molState.updateScatteringFit(ed, params.mixtureList, params.kmin, params.kmaxCurr); 
                }
		// std::cout << "Hydration density parameter C2: " << newmolState.C2 << " \n";

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
  std::string scatterNameEnd;
  if ((strcmp(argv[19], "True") == 0)) {
    scatterNameEnd = write_scatter_ChiSq(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr,params.mixtureList, "end");
  }else{
    scatterNameEnd = write_scatter(argv[12], improvementIndex, molState, ed, params.kmin, params.kmaxCurr,params.mixtureList, "end");
  }
  std::cout << "\n best overall mol name: " << moleculeNameEnd << "\n";
  std::cout << " overallFitBest fit: " << overallFit.first << "\n";

  logger.logEntry(improvementIndex, fitStep, overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd, molState.C2);

} // end of main
