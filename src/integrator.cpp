#include "../include/integrator.h"

#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"
#include "../include/mathUtils.h"

const double INTEGRAL_NORMALIZATION_FACTOR = 0.5;

double Integrator::integrateForMassMatrix(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int, int>& integral_parameters) {
    int const count_of_segments = 20;

    int const index_of_element = std::get<0>(integral_parameters);
    int const index_of_first_basis_function = std::get<1>(integral_parameters);
    int const index_of_second_basis_function = std::get<2>(integral_parameters);

    double const x_integral_step = 1.0 / count_of_segments;
    double const y_integral_step = 1.0 / count_of_segments;
    double const z_integral_step = 1.0 / count_of_segments;

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    double integral = 0.0;

    for (int i = 0; i <= count_of_segments; i++) {
        for (int j = 0; j <= count_of_segments; j++) {
#pragma unroll 4
            for (int k = 0; k <= count_of_segments; k++) {
                Point const point = {nodes[0].x + k * x_integral_step,
                                     nodes[0].y + j * y_integral_step,
                                     nodes[0].z + i * z_integral_step};

                double weight = 1.0;

                if (i == 0 || i == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (j == 0 || j == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (k == 0 || k == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                integral +=
                    weight *
                    MathUtils::getLinearBasisFunction(
                        point, nodes[node_indexes.at(0)],
                        nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                        nodes[node_indexes.at(3)], nodes[nodes.size() - 1],
                        index_of_first_basis_function) *
                    MathUtils::getLinearBasisFunction(
                        point, nodes[node_indexes.at(0)],
                        nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                        nodes[node_indexes.at(3)], nodes[nodes.size() - 1],
                        index_of_second_basis_function);
            }
        }
    }

    integral *= x_integral_step * y_integral_step * z_integral_step;

    return integral;
}

double Integrator::integrateForRightPartVector(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int>& integral_parameters) {
    int const count_of_segments = 20;

    int const index_of_element = std::get<0>(integral_parameters);
    int const index_of_basis_function = std::get<1>(integral_parameters);

    double const x_integral_step = 1.0 / count_of_segments;
    double const y_integral_step = 1.0 / count_of_segments;
    double const z_integral_step = 1.0 / count_of_segments;

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    double integral = 0.0;

    for (int i = 0; i <= count_of_segments; i++) {
        for (int j = 0; j <= count_of_segments; j++) {
#pragma unroll 4
            for (int k = 0; k <= count_of_segments; k++) {
                Point const point = {nodes[0].x + k * x_integral_step,
                                     nodes[0].y + j * y_integral_step,
                                     nodes[0].z + i * z_integral_step};

                double weight = 1.0;

                if (i == 0 || i == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (j == 0 || j == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (k == 0 || k == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                integral +=
                    weight *
                    MathUtils::getLinearBasisFunction(
                        point, nodes[node_indexes.at(0)],
                        nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                        nodes[node_indexes.at(3)], nodes[nodes.size() - 1],
                        index_of_basis_function) *
                    MathUtils::calculateF(point);
            }
        }
    }

    integral *= x_integral_step * y_integral_step * z_integral_step;

    return integral;
}

double Integrator::integrateForStiffnessMatrix(
    const std::vector<Point>& nodes,
    const std::vector<Element>& finite_elements,
    const std::tuple<int, int, int>& integral_parameters) {
    int const count_of_segments = 20;

    int const index_of_element = std::get<0>(integral_parameters);
    int const index_of_first_basis_function = std::get<1>(integral_parameters);
    int const index_of_second_basis_function = std::get<2>(integral_parameters);

    double const x_integral_step = 1.0 / count_of_segments;
    double const y_integral_step = 1.0 / count_of_segments;
    double const z_integral_step = 1.0 / count_of_segments;

    double first_integral_component = 0.0;
    double second_integral_component = 0.0;
    double third_integral_component = 0.0;
    double fourth_integral_component = 0.0;
    double fifth_integral_component = 0.0;
    double sixth_integral_component = 0.0;

    const auto node_indexes =
        finite_elements[index_of_element].getNodeIndexes();

    double integral = 0.0;

    for (int i = 0; i <= count_of_segments; i++) {
        for (int j = 0; j <= count_of_segments; j++) {
#pragma unroll 4
            for (int k = 0; k <= count_of_segments; k++) {
                Point const point = {nodes[0].x + k * x_integral_step,
                                     nodes[0].y + j * y_integral_step,
                                     nodes[0].z + i * z_integral_step};

                double weight = 1.0;

                if (i == 0 || i == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (j == 0 || j == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                if (k == 0 || k == count_of_segments) {
                    weight *= INTEGRAL_NORMALIZATION_FACTOR;
                }

                second_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_second_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 0);

                fourth_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_second_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 1);

                sixth_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_second_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 2);

                first_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_first_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 0);

                third_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_first_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 1);

                fifth_integral_component =
                    MathUtils::getDerivativeLinearBasisFunction(
                        index_of_first_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        nodes[nodes.size() - 1], 2);

                integral +=
                    weight *
                    (first_integral_component * second_integral_component +
                     third_integral_component * fourth_integral_component +
                     fifth_integral_component * sixth_integral_component);
            }
        }
    }

    integral *= x_integral_step * y_integral_step * z_integral_step;

    return integral;
}
