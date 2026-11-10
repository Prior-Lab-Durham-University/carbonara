#include "Logger.h"

Logger::Logger(const std::string& filePath) {

    logFile.open(filePath, std::ios::out | std::ios::app);

    if (!logFile.is_open()) {
        std::cerr << "Failed to open log file: " << filePath << std::endl;
    }

    start = std::chrono::high_resolution_clock::now();
}

Logger::~Logger() {
    if (logFile.is_open()) {
        logFile.close();
    }
}

long long Logger::getElapsedTime() const {

        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(now - start).count();
    }

void Logger::logEntry(const int& improvementIndex, const int& fitStep, const double& scatterFitFirst,
                      const double& writhePenalty, const double& overlapPenalty, const double& distanceConstraints,
                      const double& kmaxCurr, const std::string& scatterPath, const std::string& moleculePath,
                      const double& C2) {

    if (logFile.is_open()) {

        long long elapsed = getElapsedTime();

        logFile << "{";
        logFile << "\"ImprovementIndex\": " << improvementIndex << ", ";
        logFile << "\"FitStep\": " << fitStep << ", ";
        logFile << "\"ScatterFitFirst\": " << scatterFitFirst << ", ";
        logFile << "\"WrithePenalty\": " << writhePenalty << ", ";
        logFile << "\"OverlapPenalty\": " << overlapPenalty << ", ";
        logFile << "\"DistanceConstraints\": " << distanceConstraints << ", ";
        logFile << "\"HydrationDensity\": " << C2 << ", ";
        logFile << "\"ElapsedTime(µs)\": " << elapsed << ", ";
        logFile << "\"KmaxCurr\": " << kmaxCurr << ", ";
        logFile << "\"ScatterPath\": \"" << escapeJSON(scatterPath) << "\"" << ", ";
        logFile << "\"MoleculePath\": \"" << escapeJSON(moleculePath) << "\"";
        logFile << "}\n";

        logFile.flush();
    }
}

void Logger::logMetadata(const std::string& run, ModelParameters params) {

    if (logFile.is_open()) {

        logFile << "{";
        logFile << "\"Run\": \"" << escapeJSON(run) << "\", ";
        logFile << "\"qmin\": \"" << params.kmin << "\", ";
        logFile << "\"qmax\": \"" << params.kmax << "\", ";
        logFile << "\"qmax_curr\": \"" << params.kmaxCurr << "\", ";

        logFile << "\"HelixSolventRatios\": [";
        for (size_t i = 0; i < params.helRatList.size(); ++i) {
            logFile << params.helRatList[i];
            if (i < params.helRatList.size() - 1) {
                logFile << ", ";
            }
        }
        logFile << "], ";

        logFile << "\"lmin\":" << params.lmin << ", ";
        logFile << "\"rmin\":" << params.rmin << ", ";
        logFile << "\"rmax\":" << params.rmax << ", ";
        logFile << "\"closestApproachDist\":" << params.closestApproachDist << ", ";
        logFile << "\"lmin\":" << params.lmin << ", ";
        logFile << "\"Rin\":" << params.Rin << ", ";
        logFile << "\"Rout\":" << params.Rout << ", ";
        logFile << "\"RShell\":" << params.RShell << ", ";
        logFile << "\"ntrivs\":" << params.ntrivs << ", ";
        logFile << "\"solventsPerLink\":" << params.solventsPerLink << "}";
        logFile << "\n";

        logFile.flush();
    }
}

std::string Logger::escapeJSON(const std::string& s) {

    std::ostringstream o;

    for (auto c : s) {

        switch (c) {
            case '"': o << "\\\""; break;
            case '\\': o << "\\\\"; break;
            case '\b': o << "\\b"; break;
            case '\f': o << "\\f"; break;
            case '\n': o << "\\n"; break;
            case '\r': o << "\\r"; break;
            case '\t': o << "\\t"; break;
            default: o << c; break;
        }
    }

    return o.str();
}

void Logger::consoleInitial(const double& scatterFitFirst, const double& writhePenalty,
                            const double& overlapPenalty, const double& distanceConstraints) {

    // std::cout << std::left << std::setw(80) << "------------------------------ Initial Molecule ------------------------------- "
    std::cout << std::left << std::setw(80) << "------------------------------------------------------------------------------- "
                           << "\n";

    std::cout << std::left << std::setw(20) << "Scattering Fit"
                           << std::setw(20) << "Overlap Penalty"
                           << std::setw(20) << "Writhe Penalty"
                           << std::setw(20) << "Contact Penalty"
                           << "\n";

    std::cout << std::left << std::setw(20) << scatterFitFirst
                           << std::setw(20) << overlapPenalty
                           << std::setw(20) << writhePenalty
                           << std::setw(20) << distanceConstraints
                           << "\n"
                           << "------------------------------------------------------------------------------- "
                           << "\n";

}

void Logger::consoleCurrentStep(int step, int index, double currFit) {

    std::cout << std::left << std::setw(10) << "Step"
                            << std::setw(10) << step
                            << std::setw(10) << "Index"
                            << std::setw(10) << index
                            << std::setw(15) << "Current Fit"
                            << std::setw(10) << std::setprecision(5) << currFit
                            << "\n";

}

void Logger::consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond) {

    if (step==0){
        std::cout << std::left << std::setw(20) << "Improvement Index"
                               << std::setw(20) << "Current Fit Step"
                               << std::setw(20) << "Current Scatter"
                               << std::setw(20) << "Updated Scatter"
                               << "\n";

    }

    std::cout << std::left

                << std::setw(20)  << improveIndex
                << std::setw(20)  << step
                << std::setw(20) << scatterFitFirst
                << std::setw(20) << scatterFitSecond
                << "\n";
    }


void Logger::consoleChange(std::string updateType, ModelParameters& params) {

    if (updateType=="fitImprove") {
        std::cout << std::left << "                  --- Scatter Fit Improved! --- \n";
    }

    else if (updateType=="krangeIncrease") {
        std::cout << std::left << "                   --- K max increased to " << params.kmaxCurr << " --- \n";
    }

}
