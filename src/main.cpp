#include <string_view>
#include <vector>

#include "../include/fem.h"
#include "../include/grid.h"

static constexpr std::string_view GRID_INPUT_FILE_NAME =
    "../data/input/gridInfo.txt";
static constexpr std::string_view BOUNDARIES_INPUT_FILE_NAME =
    "../data/input/boundaries.txt";

int main() {
    const std::vector<Point> test_points = {{5E-1, 5E-1, 0E-1},
                                            {0E-1, 5E-1, 0E-1},
                                            {5E-1, 5E-1, 10E-1},
                                            {5E-1, 10E-1, 75E-2}};

    FEM fem;

    fem.generateFEMData(GRID_INPUT_FILE_NAME);

    fem.inputBoundaryConditions(BOUNDARIES_INPUT_FILE_NAME);

    fem.solveFEM();

    fem.saveTestResults(test_points);
}