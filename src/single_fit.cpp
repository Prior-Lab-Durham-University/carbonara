#include <iostream>
#include <vector>
#include <string>
#include "ktlMoleculeRandom.h"
#include "experimentalData.h"
#include "moleculeFitAndState.h"
#include "helpers.h"

// Function to create a horizontal line
std::string horizontal_line(int width = 50) {
    return std::string(width, '-');
}

// Function to center a string
std::string center_string(const std::string& str, int width) {
    int padding = width - str.length();
    int left_padding = padding / 2;
    return std::string(left_padding, ' ') + str;
}

int main(int argc, const char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <scattering_data_file> <fingerprint_file> <coordinate_file> <output_prefix> <qmin> <q_max>" << std::endl;
        return 1;
    }

    // Set up model parameters
    ModelParameters params;
    params.kmin = std::stod(argv[5]);
    params.kmax = std::stod(argv[6]);
    params.kmaxCurr = params.kmax;
    params.rmin = 3.7;
    params.rmax = 3.9;
    params.lmin = 4.0;




    // Read experimental data
    experimentalData ed(argv[1]);

    // Initialize molecule
    ktlMolecule molecule;
    molecule.readInSequence(argv[2], params.rmin, params.rmax, params.lmin);
    molecule.readInCoordinates(argv[3]);
    molecule.getHydrophobicResidues();

    // Create vector of molecules (in this case, just one)
    std::vector<ktlMolecule> molecules = {molecule};

    // Initialize molecule fit and state
    moleculeFitAndState molState(molecules, params);

    // Perform single fit
    std::vector<std::vector<double>> dummyMixtureList = {{1.0}}; // 100% of single structure
    std::vector<double> dummyHelRatList = {0.5}; // Dummy helix ratio
    std::pair<double, double> fit = molState.getOverallFit_ChiSq(ed, dummyMixtureList, params.kmin, params.kmaxCurr);

    std::cout<<"fit quality "<<fit.second<<"\n";
			    
			    
    /* Output results
    std::cout << horizontal_line() << std::endl;
    std::cout << center_string("Fit Results", 50) << std::endl;
    std::cout << horizontal_line() << std::endl;
    std::cout << std::left << std::setw(20) << "Overall fit:" << std::setw(30) << fit.first << std::endl;
    std::cout << std::left << std::setw(20) << "Scattering fit:" << std::setw(30) << fit.second << std::endl;
    std::cout << horizontal_line() << std::endl;
    */

    std::string outputPrefix = argv[4];
    std::string logFile = outputPrefix + "initialScatter.dat";
    ed.writeScatteringToFile_ChiSq(dummyMixtureList,logFile.c_str());
    /* Write fitted structure and scattering
    //std::string moleculeName = write_molecules(outputPrefix, 0, molecules, "fitted");
    std::string scatterName = write_scatter(outputPrefix, 0, molState, ed, params.kmin, params.kmaxCurr,params.mixtureList, "fitted");

    // Log it
    std::string logFile = outputPrefix + "_logFile.dat";
    Logger logger(logFile);
    logger.logEntry(0, 0, fit.first, molState.getWrithePenalty(), molState.getOverlapPenalty(),
                    molState.getDistanceConstraints(), params.kmaxCurr, scatterName, moleculeName,
                    molState.C2);

    // Output file information
    std::cout << center_string("Output Files", 50) << std::endl;
    std::cout << horizontal_line() << std::endl;
    std::cout << std::left << std::setw(25) << "Fitted structure:" << moleculeName << std::endl;
    std::cout << std::left << std::setw(25) << "Fitted scattering:" << scatterName << std::endl;
    std::cout << horizontal_line() << std::endl;

    */
    return 0;
}
