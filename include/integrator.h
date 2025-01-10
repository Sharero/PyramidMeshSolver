#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"

class Integrator {
   public:
    static double integrateBasisFunctions(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int, int>& integral_parameters);

    static double integrateDerivativeBasisFunctions(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int, int>& integral_parameters);

    static double integrateBasisFunctionForF(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int>& integral_parameters);
};
#endif
