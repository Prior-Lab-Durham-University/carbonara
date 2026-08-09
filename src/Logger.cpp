#include "Logger.h"
#include <unistd.h>
#include <cmath>
#include <algorithm>

bool Logger::colorEnabled() const {
    // Cheap, cached-per-call check rather than a static: this only runs a couple of
    // times per accepted move, nowhere near hot enough to matter, and avoids any
    // static-init-order question with iostream.
    return isatty(fileno(stdout));
}

std::string Logger::colorize(const std::string& text, const char* ansiCode) const {
    if (!colorEnabled()) { return text; }
    return std::string("\033[") + ansiCode + "m" + text + "\033[0m";
}

// Unicode block-element sparkline (▁▂▃▄▅▆▇█), scaled to the min/max of `values`.
// Flat (min==max, including a single-point history) renders as a mid-height line
// rather than dividing by zero.
std::string Logger::renderSparkline(const std::vector<double>& values) const {
    if (values.empty()) { return ""; }
    static const char* blocks[8] = {"▁","▂","▃","▄","▅","▆","▇","█"};
    double lo = *std::min_element(values.begin(), values.end());
    double hi = *std::max_element(values.begin(), values.end());
    std::string out;
    for (double v : values) {
        int level;
        if (hi - lo < 1e-12) {
            level = 3; // flat history -- neutral mid-height bar, not a divide-by-zero
        } else {
            double frac = (v - lo) / (hi - lo);
            level = (int)std::round(frac * 7.0);
            level = std::max(0, std::min(7, level));
        }
        out += blocks[level];
    }
    return out;
}

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
                      const double& C2, const double& writheDiffPenalty, const double& chi2,
                      const bool& hardSatisfied, const double& hardMaxViolation) {

    if (logFile.is_open()) {

        long long elapsed = getElapsedTime();

        logFile << "{";
        logFile << "\"ImprovementIndex\": " << improvementIndex << ", ";
        logFile << "\"FitStep\": " << fitStep << ", ";
        logFile << "\"ScatterFitFirst\": " << scatterFitFirst << ", ";
        logFile << "\"Chi2\": " << chi2 << ", ";
        logFile << "\"HardConstraintsSatisfied\": " << (hardSatisfied ? "true" : "false") << ", ";
        logFile << "\"HardConstraintsMaxViolation\": " << hardMaxViolation << ", ";
        logFile << "\"WrithePenalty\": " << writhePenalty << ", ";
        logFile << "\"OverlapPenalty\": " << overlapPenalty << ", ";
        logFile << "\"DistanceConstraints\": " << distanceConstraints << ", ";
        logFile << "\"WritheDiffPenalty\": " << writheDiffPenalty << ", ";
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
                            const double& overlapPenalty, const double& distanceConstraints,
                            const double& chi2, const bool& hardSatisfied, const double& hardMaxViolation) {

    chi2History.clear(); // fresh run -- forget any previous run's sparkline history

    std::ios_base::fmtflags savedFlags(std::cout.flags());
    std::streamsize savedPrecision = std::cout.precision();
    std::cout << std::fixed << std::setprecision(4);

    std::string title = colorize(" Initial Molecule ", "1;36");
    std::cout << "\n" << std::string(35, '-') << title << std::string(35, '-') << "\n";

    // Note: "Combined Fit" is chi2 * (1 + penaltyWeight * (overlap + distance + writhe +
    // writheDiff)) -- NOT the scattering fit alone. "Chi2" below it is the actual
    // unpenalized scattering agreement; watch that one to judge fit quality, the combined
    // column is what the search itself optimizes.
    std::cout << std::left << std::setw(15) << "Combined Fit"
                           << std::setw(15) << "Chi2"
                           << std::setw(15) << "Overlap Pen."
                           << std::setw(15) << "Writhe Pen."
                           << std::setw(15) << "Contact Pen."
                           << std::setw(15) << "Hard Constr."
                           << "\n";

    std::string hardCell = hardSatisfied
        ? colorize("OK", "1;32")
        : colorize("VIOLATED(" + std::to_string(hardMaxViolation) + ")", "1;31");

    std::cout << std::left << std::setw(15) << scatterFitFirst
                           << std::setw(15) << chi2
                           << std::setw(15) << overlapPenalty
                           << std::setw(15) << writhePenalty
                           << std::setw(15) << distanceConstraints
                           << hardCell
                           << "\n"
                           << std::string(89, '-') << "\n";

    std::cout.flags(savedFlags);
    std::cout.precision(savedPrecision);
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

void Logger::consoleFitAttempt(int step, int improveIndex, ModelParameters params, double scatterFitFirst, double scatterFitSecond,
                               const bool& hardSatisfied, const double& hardMaxViolation) {

    // Note: scatterFitFirst is the Combined Fit (chi2 * (1 + penaltyWeight * penalties)),
    // scatterFitSecond is the actual Chi2 -- despite the parameter names, matching the
    // logEntry/JSON convention. This is the *current accepted* state after this fit step,
    // not a proposed-vs-updated pair.
    std::ios_base::fmtflags savedFlags(std::cout.flags());
    std::streamsize savedPrecision = std::cout.precision();
    std::cout << std::fixed << std::setprecision(4);

    if (step==0){
        std::cout << std::left << std::setw(15) << "Improve Idx"
                               << std::setw(15) << "Fit Step"
                               << std::setw(15) << "Combined Fit"
                               << std::setw(15) << "Chi2"
                               << std::setw(15) << "Hard Constr."
                               << "\n";

    }

    std::string hardCell = hardSatisfied
        ? colorize("OK", "1;32")
        : colorize("VIOLATED(" + std::to_string(hardMaxViolation) + ")", "1;31");

    std::cout << std::left
                << std::setw(15)  << improveIndex
                << std::setw(15)  << step
                << std::setw(15) << scatterFitFirst
                << std::setw(15) << scatterFitSecond
                << hardCell
                << "\n";

    // Rolling Chi2 sparkline -- a quick "is this actually converging" glance without
    // having to eyeball a column of numbers scrolling past. Only drawn once there's
    // enough history to be worth looking at.
    chi2History.push_back(scatterFitSecond);
    if ((int)chi2History.size() > sparklineWidth) {
        chi2History.erase(chi2History.begin());
    }
    if (chi2History.size() >= 3) {
        double lo = *std::min_element(chi2History.begin(), chi2History.end());
        double hi = *std::max_element(chi2History.begin(), chi2History.end());
        std::cout << "  " << colorize("Chi2 trend:", "2") << " " << renderSparkline(chi2History)
                  << "  (" << lo << " .. " << hi << ")\n";
    }

    std::cout.flags(savedFlags);
    std::cout.precision(savedPrecision);
    }


void Logger::consoleChange(std::string updateType, ModelParameters& params) {

    if (updateType=="fitImprove") {
        std::cout << std::left << "                  --- Scatter Fit Improved! --- \n";
    }

    else if (updateType=="krangeIncrease") {
        std::cout << std::left << "                   --- K max increased to " << params.kmaxCurr << " --- \n";
    }

}
