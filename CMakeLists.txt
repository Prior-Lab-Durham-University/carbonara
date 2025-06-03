cmake_minimum_required(VERSION 3.10)
project(Carbonara)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Include directories
include_directories(src)

# Source files
set(SOURCES
    src/point.cpp
    src/polyHelix.cpp
    src/randomMolGen.cpp
    src/ktlMoleculeRandom.cpp
    src/experimentalData.cpp
    src/skmt.cpp
    src/writheFP.cpp
    src/moleculeFitAndState.cpp
    src/Logger.cpp
    src/helpers.cpp
    src/parameters.cpp
)

# Main executables
add_executable(predictStructureQvary ${SOURCES} src/mainPredictionFinalQvar.cpp)
add_executable(generate_structure ${SOURCES} src/Flexible_generator.cpp)
add_executable(single_fit ${SOURCES} src/single_fit.cpp)

# Set output directory for executables
set_target_properties(predictStructureQvary generate_structure single_fit
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
