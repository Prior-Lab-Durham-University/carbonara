#ifndef LOGGER
#define LOGGER

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "parameters.h"
#include "moleculeFitAndState.h"


class Logger {

private:

    std::ofstream logFile;
    std::string escapeJSON(const std::string& s);
    std::chrono::high_resolution_clock::time_point start;

    // Rolling window of recent *accepted* Chi2 values, most recent last, for the
    // console sparkline in consoleFitAttempt. Capped at sparklineWidth entries --
    // this is display-only state, never written to the JSON log.
    std::vector<double> chi2History;
    static const int sparklineWidth = 30;
    std::string renderSparkline(const std::vector<double>& values) const;
    // ANSI colour is only ever emitted when stdout is an actual terminal (isatty),
    // so redirecting a run's output to a file or another program never ends up
    // full of escape codes.
    bool colorEnabled() const;
    std::string colorize(const std::string& text, const char* ansiCode) const;

public:

    Logger(const std::string& filePath);
    ~Logger();

    long long getElapsedTime() const;

    // scatterFitFirst is the fully combined, penalty-weighted objective
    // (chi2 * (1 + penaltyWeight * (overlap + distance + writhe + writheDiff))) --
    // NOT raw chi2, despite the name (kept as-is for backward compatibility with
    // existing log readers). chi2 is the actual unpenalized scattering fit; pass
    // -1.0 (the default) if unavailable, so old call sites remain valid.
    // hardSatisfied/hardMaxViolation report the *hard* (feasibility-filter) distance
    // constraints -- disulfides, strict posts -- which never contribute to
    // scatterFitFirst/distanceConstraints at all (see applyHardConstraints); default
    // true/0.0 so old call sites remain valid.
    void logEntry(const int& improvementIndex, const int& fitStep, const double& scatterFitFirst,
                  const double& writhePenalty, const double& overlapPenalty, const double& distanceConstraints,
                  const double& kmaxCurr, const std::string& scatterPath, const std::string& moleculePath,
                  const double& C2, const double& writheDiffPenalty = 0.0, const double& chi2 = -1.0,
                  const bool& hardSatisfied = true, const double& hardMaxViolation = 0.0);

    void logMetadata(const std::string& run, ModelParameters params);

    void consoleInitial(const double& scatterFitFirst, const double& writhePenalty,
                        const double& overlapPenalty, const double& distanceConstraints,
                        const double& chi2 = -1.0, const bool& hardSatisfied = true,
                        const double& hardMaxViolation = 0.0);

    void consoleCurrentStep(int step, int index, double currFit);

    void consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond,
                           const bool& hardSatisfied = true, const double& hardMaxViolation = 0.0);

    void consoleChange(std::string updateType, ModelParameters& params);

};

#endif
