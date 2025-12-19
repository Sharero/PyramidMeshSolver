#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"

using namespace boost;

class Integrator {
   public:
    static double integrateForMassMatrix(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int, int>& integral_parameters,
        BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_element_type);

    static double integrateForRightPartVector(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int>& integral_parameters,
        BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_element_type);

    static double integrateForStiffnessMatrix(
        const std::vector<Point>& nodes,
        const std::vector<Element>& finite_elements,
        const std::tuple<int, int, int>& integral_parameters,
        BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_element_type);
};
#endif
