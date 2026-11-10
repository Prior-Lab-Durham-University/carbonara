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


std::string constructMoleculeName_search(const std::string& basePath, int varying_index, const int& structure_number) {

    std::stringstream ss;

    ss << basePath << "_xyz_changeindex_" << varying_index << "_structure_" << structure_number << ".dat"; 
    
    return ss.str();

}

std::string constructScatterName_search(const std::string& basePath, int varying_index, const int& structure_number) {

    std::stringstream ss;

    ss << basePath << "_scatter_changeindex_" << varying_index << "_structure_" << structure_number << ".dat"; 

    return ss.str();

}


std::string write_molecules_search(const std::string& basePath, const int& structure_num, ktlMolecule mol, int varying_index) {
    
    std::string moleculeName; 

    // for(int i=0; i<mol.size(); i++) {  

    moleculeName = constructMoleculeName_search(basePath, varying_index, structure_num);
    mol.writeMoleculeToFile(moleculeName.c_str());
    // } 

    return moleculeName;
}


std::string write_scatter_search(const std::string& basePath, const int& structure_num, int varying_index, 
                                 moleculeFitAndState& molFit, experimentalData& ed, double kmin, double kmaxCurr) {
    
    std::string scatterName;

    scatterName = constructScatterName_search(basePath, varying_index, structure_num);

    molFit.writeScatteringToFile(ed, kmin, kmaxCurr, scatterName.c_str());

    return scatterName;

}



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

  std::pair<double, double> overallFit;
  overallFit = molState.getOverallFit(ed, params.mixtureList, params.helRatList, params.kmin, params.kmaxCurr);

  /* Main algorithm */

  std::vector<int> numberOfChainsInEachStructure = findNumberSections(moleculeStructures);

  // This is a monster while loop - strap in chaps

    for (int structureIndex = 0; structureIndex < moleculeStructures.size(); structureIndex++) {
      int netIndex = 0;

      for (int chainNumber = 1; chainNumber <= numberOfChainsInEachStructure[structureIndex]; chainNumber++) {

        // net index tells us how far we are through the whole molecule
        if (chainNumber > 1) {
          netIndex = netIndex + moleculeStructures[structureIndex].getSubsecSize(chainNumber - 1);
        }

        // Now loop over the secondary structures of the given unit or section
        for (int secondarySectionIndex = 0; secondarySectionIndex < moleculeStructures[structureIndex].getSubsecSize(chainNumber) - 1; secondarySectionIndex++) {

          int totalIndex = netIndex + secondarySectionIndex;

          // in this if statement we check which secondary sections are being changed
          if (std::count(vary_sec_list_list[structureIndex].begin(), vary_sec_list_list[structureIndex].end(), totalIndex)) {
            
            int num_new_structure = 0;

            while (num_new_structure < 200) {

                int change_section = totalIndex - netIndex;
                // std::vector<ktlMolecule> moleculeStructuresCopy = moleculeStructures;
                ktlMolecule newMol = moleculeStructures[structureIndex];

                bool cacaDist = modifyMolecule(newMol, moleculeStructures[structureIndex], change_section, chainNumber);

                if (cacaDist == false) {

                    moleculeFitAndState newmolState = molState;

                    // moleculeStructuresCopy[structureIndex] = newMol;

                    // calculate the fitting of changed molecule
                    std::pair<double, double> newOverallFit = newmolState.getOverallFit(ed, params.mixtureList, params.helRatList, newMol, params.kmin, params.kmaxCurr, structureIndex);

                    std::string moleculeName = write_molecules_search(argv[12], num_new_structure, newMol, change_section);
                    std::string scatterName = write_scatter_search(argv[12], num_new_structure, change_section, molState, ed, params.kmin, params.kmaxCurr);

                    // add to a text file the moleculeName, scatterName and newOverallFit, not using the logger class
                    std::ofstream logFile;
                    logFile.open(argv[16], std::ios::app);
                    logFile << num_new_structure << " " << change_section << " " <<  moleculeName << " " << scatterName << " " << newOverallFit.first << " " << newOverallFit.second << "\n";
                    logFile.close();

                    num_new_structure++;

                }
            }

          } // totalIndex an allowed varying section?

        } // structureIndex
      } // chainNumber
    } // structureIndex
}

