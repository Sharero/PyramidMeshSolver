#include "../include/integrator.h"

#include <iostream>
#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"
#include "../include/mathUtils.h"

double Integrator::integrateForMassMatrix(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int, int>& integral_parameters,
    BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_element_type) {
    const int index_of_element = std::get<0>(integral_parameters);
    const int index_of_first_basis_function = std::get<1>(integral_parameters);
    const int index_of_second_basis_function = std::get<2>(integral_parameters);

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    Point node_0, node_1, node_2, node_3, vertex;

    MathUtils::getMainPyramideNodes(node_indexes, nodes, basis_functions_type,
                                    basis_functions_element_type, node_0,
                                    node_1, node_2, node_3, vertex);

    const std::vector<double> gauss_nodes = {
        0.025446, 0.129234, 0.297077, 0.500000, 0.702923, 0.870766, 0.974553};

    const std::vector<double> gauss_weights = {
        0.064742, 0.139853, 0.190915, 0.208980, 0.190915, 0.139853, 0.064742};

    double result = 0.0;

    for (int zi_idx = 0; zi_idx < gauss_nodes.size(); ++zi_idx) {
        double z = gauss_nodes[zi_idx];
        double w_z = gauss_weights[zi_idx];

        double L = 1.0 - z;
        double a = z / 2.0;

        for (int xi_idx = 0; xi_idx < gauss_nodes.size(); ++xi_idx) {
            double x = a + L * gauss_nodes[xi_idx];
            double w_x = L * gauss_weights[xi_idx];

            for (int yi_idx = 0; yi_idx < gauss_nodes.size(); ++yi_idx) {
                double y = a + L * gauss_nodes[yi_idx];
                double w_y = L * gauss_weights[yi_idx];

                result +=
                    MathUtils::getBasisFunction(
                        nodes, {x, y, z}, node_indexes, basis_functions_type,
                        basis_functions_element_type,
                        index_of_first_basis_function, false) *
                    MathUtils::getBasisFunction(
                        nodes, {x, y, z}, node_indexes, basis_functions_type,
                        basis_functions_element_type,
                        index_of_second_basis_function, false) *
                    w_x * w_y * w_z;
            }
        }
    }

    std::vector<std::vector<double>> J(3, std::vector<double>(3));

    MathUtils::calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex, J);

    const double determinant = MathUtils::calculateJacobian(J);

    return std::abs(determinant) * result;
}

double Integrator::integrateForRightPartVector(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int>& integral_parameters,
    BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_element_type) {
    int const index_of_element = std::get<0>(integral_parameters);
    int const index_of_basis_function = std::get<1>(integral_parameters);

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    Point node_0, node_1, node_2, node_3, vertex;

    MathUtils::getMainPyramideNodes(node_indexes, nodes, basis_functions_type,
                                    basis_functions_element_type, node_0,
                                    node_1, node_2, node_3, vertex);

    const std::vector<double> gauss_nodes = {
        0.025446, 0.129234, 0.297077, 0.500000, 0.702923, 0.870766, 0.974553};

    const std::vector<double> gauss_weights = {
        0.064742, 0.139853, 0.190915, 0.208980, 0.190915, 0.139853, 0.064742};

    double result = 0.0;

    for (int zi_idx = 0; zi_idx < gauss_nodes.size(); ++zi_idx) {
        double z = gauss_nodes[zi_idx];
        double w_z = gauss_weights[zi_idx];

        double L = 1.0 - z;
        double a = z / 2.0;

        for (int xi_idx = 0; xi_idx < gauss_nodes.size(); ++xi_idx) {
            double x = a + L * gauss_nodes[xi_idx];
            double w_x = L * gauss_weights[xi_idx];

            for (int yi_idx = 0; yi_idx < gauss_nodes.size(); ++yi_idx) {
                double y = a + L * gauss_nodes[yi_idx];
                double w_y = L * gauss_weights[yi_idx];

                double physical_x = 0.0;
                double physical_y = 0.0;
                double physical_z = 0.0;

                MathUtils::calculatePhysicalCoordinates(
                    {x, y, z}, node_0, node_1, node_2, node_3, vertex,
                    physical_x, physical_y, physical_z);

                result +=
                    MathUtils::calculateF(
                        {physical_x, physical_y, physical_z}) *
                    MathUtils::getBasisFunction(
                        nodes, {x, y, z}, node_indexes, basis_functions_type,
                        basis_functions_element_type, index_of_basis_function,
                        false) *
                    w_x * w_y * w_z;
            }
        }
    }

    std::vector<std::vector<double>> J(3, std::vector<double>(3));

    MathUtils::calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex, J);

    const double determinant = MathUtils::calculateJacobian(J);

    return std::abs(determinant) * result;
}

double Integrator::integrateForStiffnessMatrix(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int, int>& integral_parameters,
    BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_element_type) {
    int const index_of_element = std::get<0>(integral_parameters);
    int const index_of_first_basis_function = std::get<1>(integral_parameters);
    int const index_of_second_basis_function = std::get<2>(integral_parameters);

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    const std::vector<double> gauss_nodes = {
        0.025446, 0.129234, 0.297077, 0.500000, 0.702923, 0.870766, 0.974553};

    const std::vector<double> gauss_weights = {
        0.064742, 0.139853, 0.190915, 0.208980, 0.190915, 0.139853, 0.064742};

    double result = 0.0;

    for (int zi_idx = 0; zi_idx < gauss_nodes.size(); ++zi_idx) {
        double z = gauss_nodes[zi_idx];
        double w_z = gauss_weights[zi_idx];

        double L = 1.0 - z;
        double a = z / 2.0;

        for (int xi_idx = 0; xi_idx < gauss_nodes.size(); ++xi_idx) {
            double x = a + L * gauss_nodes[xi_idx];
            double w_x = L * gauss_weights[xi_idx];

            for (int yi_idx = 0; yi_idx < gauss_nodes.size(); ++yi_idx) {
                double y = a + L * gauss_nodes[yi_idx];
                double w_y = L * gauss_weights[yi_idx];

                result +=
                    MathUtils::getDerivativeBasisFunction(
                        nodes, {x, y, z}, node_indexes, basis_functions_type,
                        basis_functions_element_type,
                        index_of_first_basis_function,
                        index_of_second_basis_function, false) *
                    w_x * w_y * w_z;
            }
        }
    }

    return result;
}
