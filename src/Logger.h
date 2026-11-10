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

public:

    Logger(const std::string& filePath);
    ~Logger();

    long long getElapsedTime() const;

    void logEntry(const int& improvementIndex, const int& fitStep, const double& scatterFitFirst,
                  const double& writhePenalty, const double& overlapPenalty, const double& distanceConstraints,
                  const double& kmaxCurr, const std::string& scatterPath, const std::string& moleculePath,
                  const double& C2);

    void logMetadata(const std::string& run, ModelParameters params);

    void consoleInitial(const double& scatterFitFirst, const double& writhePenalty,
                        const double& overlapPenalty, const double& distanceConstraints);

    void consoleCurrentStep(int step, int index, double currFit);

    void consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond);

    void consoleChange(std::string updateType, ModelParameters& params);

};

#endif
