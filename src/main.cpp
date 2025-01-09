#include <string_view>

#include "../include/fem.h"
#include "../include/grid.h"

static constexpr std::string_view GRID_INPUT_FILE_NAME = "../data/gridInfo.txt";
static constexpr std::string_view BOUNDARIES_INPUT_FILE_NAME =
    "../data/boundaries.txt";

int main() {
    FEM fem;

    const std::vector<Point> test_points = {{5E-1, 5E-1, 0E-1},
                                            {0E-1, 5E-1, 0E-1},
                                            {5E-1, 5E-1, 10E-1},
                                            {5E-1, 10E-1, 75E-2}};

    fem.generateLinearData(GRID_INPUT_FILE_NAME);
    fem.inputBoundaryConditions(BOUNDARIES_INPUT_FILE_NAME);
    fem.saveGridForVisualize();
    fem.generatePortrait();
    fem.assemblyGlobalComponents();
    fem.solveFEM();
    fem.printTestResults(test_points);
}