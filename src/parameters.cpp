#include "parameters.h"

ModelParameters loadParameters(const char* argv[], int argc) {

    ModelParameters params;

    // command-line arguments for scatter range
    params.kmin = std::atof(argv[8]);
    params.kmax = std::atof(argv[9]);
    params.kmaxCurr = std::atof(argv[10]);

    // If resuming from a previous run, adjust parameters accordingly.
    if (strcmp(argv[3], "True") == 0) {
        params.kmaxCurr = std::atof(argv[24]);  // Ensure this is the correct argument index for kmaxCurr.
        std::cout << "kmax curr is " << params.kmaxCurr << "\n";
        // params.improvementIndex = std::atoi(argv[17]);  // Ensure this is the correct argument index for improvementIndex.
    }

    // command-line arguments for max number of steps
    params.noScatterFitSteps = std::atoi(argv[11]);

    // command-line arguments for rigid body transformations
    if(strcmp(argv[18],"True") == 0){
        params.affineTrans=true;
    } else {
        params.affineTrans=false;
    }

    // logic for helRatList --change?
    params.helRatList.push_back(0.5);

    params.basePath = argv[12];

    // argv[20]/argv[21]: optional ensemble writhe-difference restraint.
    // Guarded on argc so any caller built against the old 19-argument
    // convention (main_multi.cpp, older scripts) is completely unaffected
    // and gets the disabled default set in ModelParameters.
    if (argc > 20 && std::strlen(argv[20]) > 0) {
        params.maxWritheDiff = std::atof(argv[20]);
    }
    if (argc > 21 && std::strlen(argv[21]) > 0) {
        int stride = std::atoi(argv[21]);
        params.writheDiffStride = stride > 0 ? stride : 1;
    }

    // argv[22]: optional penalty weight for the proportional fit objective
    // (currFit = chi2 * (1 + penaltyWeight * penalties)). Absent/empty keeps
    // the ModelParameters default above.
    if (argc > 22 && std::strlen(argv[22]) > 0) {
        params.penaltyWeight = std::atof(argv[22]);
    }

    // argv[23]: optional cap on a single soft distance-constraint pair's penalty
    // contribution. Absent/empty keeps the ModelParameters default above.
    if (argc > 23 && std::strlen(argv[23]) > 0) {
        params.distanceConstraintCap = std::atof(argv[23]);
    }

    return params;
}
