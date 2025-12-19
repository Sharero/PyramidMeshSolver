#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string_view>
#include <vector>

#include "../include/fem.h"
#include "../include/grid.h"

static constexpr std::string_view GRID_INPUT_FILE_NAME =
    "../data/input/gridInfo.txt";

static constexpr std::string_view BOUNDARIES_INPUT_FILE_NAME =
    "../data/input/boundaries.txt";

const BASIS_TYPE BASIS_FUNCTIONS_TYPE = BASIS_TYPE::Hierarhical;

const BASIS_ELEMENT_TYPE BASIS_FUNCTIONS_ELEMENTS_TYPE =
    BASIS_ELEMENT_TYPE::Quadratic;

int main() {
    const std::vector<Point> test_points = {
        {0.5, 0.5, 0.25}, {0.5, 0.5, 0.75}, {0.25, 0.5, 0.5}, {0.75, 0.5, 0.5}};

    FEM fem(BASIS_FUNCTIONS_TYPE, BASIS_FUNCTIONS_ELEMENTS_TYPE);

    fem.generateFEMData(GRID_INPUT_FILE_NAME);

    // fem.checkBasisFunctionsDeltaProperty();

    fem.inputBoundaryConditions(BOUNDARIES_INPUT_FILE_NAME);

    fem.solveFEM();

    fem.calculateResultAtPointsToFile(test_points);
}
