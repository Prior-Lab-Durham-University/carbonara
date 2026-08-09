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

    // Ensemble writhe-difference restraint (mixture_states > 1 only). <= 0
    // disables it -- the default, and what every invocation lacking argv[20]
    // (older scripts, main_multi.cpp) gets, so this is a no-op unless a
    // caller opts in explicitly. See moleculeFitAndState::applyWritheDiffConstraint.
    double maxWritheDiff = -1.0;
    int writheDiffStride = 1;

    // Weight on the soft constraint penalties (overlap, distance, writhe,
    // writheDiff) relative to chi2 in the fit objective. Combination is
    // proportional -- currFit = chi2 * (1 + penaltyWeight * penalties) --
    // so the penalties' share of the objective stays roughly constant
    // across the whole search instead of being swamped whenever chi2 is
    // large (early in a hard fit, or with a noisy/high-chi2 dataset) and
    // only starting to bind once chi2 has already dropped near convergence.
    // How "wrong" a given raw penalty sum actually is varies a lot by
    // system (number of chains/sections, how many constraints are active),
    // so this is deliberately left tunable rather than baked in -- there is
    // no single correct value across datasets.
    double penaltyWeight = 5.0;

    // Ceiling on a single *soft* distance-constraint pair's penalty contribution
    // (see ktlMolecule::getLennardJonesContact). Keeps a strict tolerance from being able
    // to produce an unbounded penalty that swamps chi2 unpredictably -- identical to the
    // uncapped quartic for small violations, saturates for large ones. Does not apply to
    // *hard* pairs (disulfides, strict posts), which are enforced as a true feasibility
    // filter instead (see moleculeFitAndState::applyHardConstraints) and are unaffected by
    // this value.
    double distanceConstraintCap = 50.0;

};

/**
 * \brief Implementation of the loadParameters function.
 * 
 * This function takes an array of command-line arguments and uses them to fill a ModelParameters object.
 * 
 * \param argv Carbonara command-line arguments.
 * \return ModelParameters object.
 */
ModelParameters loadParameters(const char* argv[], int argc);

#endif