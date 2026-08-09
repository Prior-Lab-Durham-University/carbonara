/* Carbonara Version: 0.2.0 */

#include "ktlMoleculeRandom.h"
#include "experimentalData.h"
#include <string.h>
#include "moleculeFitAndState.h"
#include <cstring>
#include <cstdlib>
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

  Optional, appended -- absent (argc <= 20/21) leaves this restraint disabled,
  so every pre-existing invocation behaves exactly as before:
  argv[20] max writhe-difference allowed between ensemble members (mixture_states > 1
           fits only); <= 0 or omitted disables it. See moleculeFitAndState::
           applyWritheDiffConstraint.
  argv[21] CA-backbone stride used for that calculation (e.g. 2 or 4); omitted
           or <= 0 means no striding (every residue).
  argv[22] penalty weight for the fit objective: currFit = chi2 * (1 +
           penaltyWeight * (overlap + distance + writhe + writheDiff
           penalties)). Proportional rather than additive, so the
           penalties' share of the objective stays roughly constant across
           the whole search instead of vanishing whenever chi2 is large.
           Omitted/empty keeps ModelParameters' default (5.0). How strict
           this should be is dataset-dependent (a high chi2 can mean "very
           wrong structure" or just "noisy SAXS data"), so tune per system.
  argv[23] cap on a single *soft* distance-constraint pair's penalty contribution
           (ktlMolecule::getLennardJonesContact). Identical to the uncapped quartic
           for small violations, saturates for large ones. Does not apply to
           *hard*-flagged pairs (disulfides, strict posts) -- those are a true
           feasibility filter (moleculeFitAndState::applyHardConstraints), rejected
           outright if violated regardless of chi2, and unaffected by this value.
           Omitted/empty keeps ModelParameters' default (50.0).
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- */

int main(int argc, const char* argv[]) {

  /* initialise the log file */
  Logger logger(argv[16]);

  /* Set up model parameters */
  ModelParameters params = loadParameters(argv, argc);

  /* Determine initial model: Two options no initial prediction, we must generate a structure
   or some initial structure provided. Actually we need a half-half option */

  /* Initialise the molecule(s) vector */
  std::vector<ktlMolecule> moleculeStructures;
  readInStructures(argv, moleculeStructures, params);

  /* Determine which sections are being altered */
  std::vector<std::vector<int>> vary_sec_list_list;
  determineVaryingSections(argv, vary_sec_list_list);

  /* Read in any fixed distances constraints (contact predictions/sulfide bonds) */
  readFixedDistancesConstraints(argv, moleculeStructures, params);

  /* Read in the permissible mixture list */
  readPermissibleMixtures(argv, params);

  /* Read in the scattering and set up the scattering model */
  experimentalData ed(argv[1]);

  // Which fitting procedure(s) this run uses -- replaces the old pattern of picking among
  // 8 differently-named getOverallFit* functions at each call site via ad-hoc if/else on
  // these same two flags. Three distinct FitModes are genuinely needed here (not one
  // reused everywhere): forceConnection tracks the *move type*, not just affineTrans --
  // the per-loop-reshape trial move never uses ForceConnection even when affineTrans is on
  // (matches the original code's actual behaviour, confirmed by reading it directly).
  const bool weightedChiSq = (strcmp(argv[19], "True") == 0);
  FitMode baseMode{weightedChiSq, params.affineTrans};   // baseline (initial + kmax-increase) computation
  FitMode rotationMode{weightedChiSq, true};             // rigid-rotation trial moves (only reached when affineTrans==true)
  FitMode reshapeMode{weightedChiSq, false};             // per-loop-reshape trial moves -- never ForceConnection

  /* Random generator */
  RandomGenerator rng;

  /* initialise the state of mol vector */
  moleculeFitAndState molState(moleculeStructures, params);

  int improvementIndex = 0;
  // If we resume from previous run - argv[3] restart True/False
  if ((strcmp(argv[3], "True") == 0)) {
    improvementIndex = std::atoi(argv[17]);
  }

  std::pair<double, double> overallFit = molState.computeOverallFit(ed, params.mixtureList, params.kmin, params.kmaxCurr, baseMode);
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
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameInitial, xyzNameInitial, molState.C2,
                  molState.getWritheDiffPenalty(), overallFit.second,
                  molState.getHardConstraintsSatisfied(), molState.getHardConstraintsMaxViolation());
  logger.consoleInitial(overallFit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(), molState.getDistanceConstraints(),
                        overallFit.second, molState.getHardConstraintsSatisfied(), molState.getHardConstraintsMaxViolation());
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

      increaseKmax(overallFit, molStateSet, ed, params, logger, baseMode);
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
            std::pair<double, double> newOverallFit = newMolState.computeOverallFit(ed, params.mixtureList, molCopyR, params.kmin, params.kmaxCurr, structureIndex, rotationMode);
            double uProb = rng.getDistributionR();
            // Hard constraints (disulfides, strict posts) are a feasibility filter, not a
            // penalty: no chi2 improvement can buy past one, so it's required in addition to
            // (not instead of) the normal combined-objective acceptance check.
            if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps) && newMolState.getHardConstraintsSatisfied()) {
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

              // calculate the fitting of changed molecule (never ForceConnection, regardless of affineTrans)
              std::pair<double, double> newOverallFit = newmolState.computeOverallFit(ed, params.mixtureList, newMol, params.kmin, params.kmaxCurr, structureIndex, reshapeMode);
	      //std::cout<<"improve ever ? "<<indexCh<<" "<<newOverallFit.second<<" "<<newOverallFit.first<<" "<<overallFit.first<<"\n";
	      // Diagnostic only, zero-cost/zero-behavior-change unless the env var is
	      // set: dumps every trial evaluation (accepted or not) so the relative
	      // scale of chi2 vs. each penalty term can be inspected across a run.
	      if (std::getenv("CARBONARA_DEBUG_FIT") != nullptr) {
	        std::cerr << "TRIAL fitStep=" << fitStep << " structureIndex=" << structureIndex
	                  << " chi2=" << newOverallFit.second
	                  << " overlap=" << newmolState.getOverlapPenalty()
	                  << " distance=" << newmolState.getDistanceConstraints()
	                  << " writhe=" << newmolState.getWrithePenalty()
	                  << " writheDiff=" << newmolState.getWritheDiffPenalty()
	                  << " currFit=" << newOverallFit.first
	                  << " prevFit=" << overallFit.first
	                  << " hardOK=" << newmolState.getHardConstraintsSatisfied()
	                  << " hardMaxViolation=" << newmolState.getHardConstraintsMaxViolation()
	                  // matches the real acceptance check below exactly -- both conditions,
	                  // not just the chi2/combined-objective comparison, since a hard
	                  // violation rejects a move regardless of how much currFit improved.
	                  << " willAccept=" << ((newOverallFit.first < overallFit.first && newmolState.getHardConstraintsSatisfied()) ? 1 : 0)
	                  << "\n";
	      }
              double uProb = rng.getDistributionR();
              // Hard constraints (disulfides, strict posts) are a feasibility filter, not a
              // penalty: no chi2 improvement can buy past one, so it's required in addition to
              // (not instead of) the normal combined-objective acceptance check.
              if (checkTransition(newOverallFit.first, overallFit.first, uProb, fitStep, params.noScatterFitSteps) && newmolState.getHardConstraintsSatisfied()) {

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
    logger.consoleFitAttempt(fitStep, improvementIndex, params, overallFit.first, overallFit.second,
                             molState.getHardConstraintsSatisfied(), molState.getHardConstraintsMaxViolation());

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
                  molState.getDistanceConstraints(), params.kmaxCurr, scatterNameEnd, moleculeNameEnd, molState.C2,
                  molState.getWritheDiffPenalty(), overallFit.second,
                  molState.getHardConstraintsSatisfied(), molState.getHardConstraintsMaxViolation());

} // end of main
