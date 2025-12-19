#include "../include/mathUtils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"
#include "../include/phi.h"

const double TETRAHEDRON_VOLUME_DIVISOR = 6.0;

void MathUtils::filterRemainingIndexesForBaseCombinations(
    int vertex_index, const std::set<int>& forbidden_indexes,
    const std::set<int>& required_indexes,
    std::vector<int>& remaining_indexes) {
#pragma unroll 4
    for (int i = 0; i < vertex_index; ++i) {
        if (forbidden_indexes.count(i) == 0 && required_indexes.count(i) == 0) {
            remaining_indexes.push_back(i);
        }
    }
}

void MathUtils::calculateCombinations(
    int combination_size, const std::set<int>& required_indexes,
    const std::vector<int>& remaining_indexes,
    std::vector<Element>& base_nodes_combinations) {
    std::vector<int> sorted_required_indexes(required_indexes.begin(),
                                             required_indexes.end());

    std::sort(sorted_required_indexes.begin(), sorted_required_indexes.end());

    const int remaining_size =
        combination_size - static_cast<int>(sorted_required_indexes.size());

    std::vector<bool> indicators(remaining_indexes.size(), false);

    std::fill(indicators.begin(), indicators.begin() + remaining_size, true);

    while (true) {
        std::vector<int> combination(sorted_required_indexes.begin(),
                                     sorted_required_indexes.end());

        int last_index = combination.empty() ? -1 : combination.back();

#pragma unroll 4
        for (size_t i = 0; i < remaining_indexes.size(); ++i) {
            if (indicators[i] && remaining_indexes[i] > last_index) {
                combination.push_back(remaining_indexes[i]);
                last_index = remaining_indexes[i];
            }
        }

        if (static_cast<int>(combination.size()) == combination_size) {
            Element element;
            element.setNodeIndexes(combination);
            base_nodes_combinations.emplace_back(std::move(element));
        }

        if (!std::prev_permutation(indicators.begin(), indicators.end())) {
            break;
        }
    }
}

bool MathUtils::isCubicNodeIllegal(
    const std::tuple<int, int, int>& cubic_node_parameters) {
    const int grid_z_index = std::get<0>(cubic_node_parameters);
    const int grid_y_index = std::get<1>(cubic_node_parameters);
    const int grid_x_index = std::get<2>(cubic_node_parameters);

    const int last_grid_x_y_z_index = 5;

    if (grid_z_index == 1 || grid_z_index == 4) {
        return (grid_x_index != 1 && grid_x_index != 4) ||
               (grid_y_index != 1 && grid_y_index != 4);
    }

    if (grid_z_index == 0 || grid_z_index == last_grid_x_y_z_index ||
        grid_z_index == 2 || grid_z_index == 3) {
        return !(grid_y_index != 1 && grid_y_index != 4 && grid_x_index != 1 &&
                 grid_x_index != 4);
    }

    return true;
}

double MathUtils::calculateF(Point point) {
    // return 1;

    // return point.x + point.y + point.z;

    return -6.0 * (point.x + point.y + point.z) + point.x * point.x * point.x +
           point.y * point.y * point.y + point.z * point.z * point.z;

    // return point.x * point.x * point.x + point.y * point.y * point.y +
    //        point.z * point.z * point.z;

    // return -12 * (point.x * point.x + point.y * point.y + point.z * point.z)
    // +
    //        point.x * point.x * point.x * point.x +
    //        point.y * point.y * point.y * point.y +
    //        point.z * point.z * point.z * point.z;

    // return point.x * point.x * point.x * point.x +
    //        point.y * point.y * point.y * point.y +
    //        point.z * point.z * point.z * point.z;
}

double MathUtils::calculateU(Point point) {
    // return 1;

    // return point.x + point.y + point.z;

    // return point.x * point.x + point.y * point.y + point.z * point.z;

    return point.x * point.x * point.x + point.y * point.y * point.y +
           point.z * point.z * point.z;

    // return point.x * point.x * point.x * point.x +
    //        point.y * point.y * point.y * point.y +
    //        point.z * point.z * point.z * point.z;
}

bool MathUtils::isNumbersEqual(double first_number, double second_number,
                               double epsilon) {
    return std::fabs(first_number - second_number) < epsilon;
}

bool MathUtils::isPlane(std::vector<Point>& nodes, double average,
                        const Element& element, char axis) {
    const auto& node_indexes = element.getNodeIndexes();

    switch (axis) {
        case 'x':
            return std::all_of(
                node_indexes.begin(), node_indexes.end(), [&](int index) {
                    return isNumbersEqual(nodes[index].x, average,
                                          NUMBERS_EQUAL_EPSILON);
                });
        case 'y':
            return std::all_of(
                node_indexes.begin(), node_indexes.end(), [&](int index) {
                    return isNumbersEqual(nodes[index].y, average,
                                          NUMBERS_EQUAL_EPSILON);
                });
        case 'z':
            return std::all_of(
                node_indexes.begin(), node_indexes.end(), [&](int index) {
                    return isNumbersEqual(nodes[index].z, average,
                                          NUMBERS_EQUAL_EPSILON);
                });
        default:
            return false;
    }
}

Point MathUtils::calculateVectorCrossProduct(Point first_vector,
                                             Point second_vector) {
    Point cross_product_result = {
        first_vector.y * second_vector.z - first_vector.z * second_vector.y,
        first_vector.z * second_vector.x - first_vector.x * second_vector.z,
        first_vector.x * second_vector.y - first_vector.y * second_vector.x};

    return cross_product_result;
}

void MathUtils::calculatePhysicalCoordinates(Point local_point, Point node_0,
                                             Point node_1, Point node_2,
                                             Point node_3, Point vertex,
                                             double& x, double& y, double& z) {
    Point const center_of_base = {
        (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
        (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
        (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    Point v1 = {node_1.x - node_0.x, node_1.y - node_0.y, node_1.z - node_0.z};
    Point v2 = {node_2.x - node_0.x, node_2.y - node_0.y, node_2.z - node_0.z};
    Point v3 = {vertex.x - center_of_base.x, vertex.y - center_of_base.y,
                vertex.z - center_of_base.z};

    double P0v1 = node_0.x * v1.x + node_0.y * v1.y + node_0.z * v1.z;
    double P0v2 = node_0.x * v2.x + node_0.y * v2.y + node_0.z * v2.z;
    double Pcv3 = center_of_base.x * v3.x + center_of_base.y * v3.y +
                  center_of_base.z * v3.z;

    double v1v1 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
    double v2v2 = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;
    double v3v3 = v3.x * v3.x + v3.y * v3.y + v3.z * v3.z;

    Point d = {local_point.x * v1v1 + P0v1, local_point.y * v2v2 + P0v2,
               local_point.z * v3v3 + Pcv3};

    std::vector<std::vector<double>> A = {
        {v1.x, v1.y, v1.z}, {v2.x, v2.y, v2.z}, {v3.x, v3.y, v3.z}};
    std::vector<std::vector<double>> Ad1 = {
        {d.x, v1.y, v1.z}, {d.y, v2.y, v2.z}, {d.z, v3.y, v3.z}};
    std::vector<std::vector<double>> Ad2 = {
        {v1.x, d.x, v1.z}, {v2.x, d.y, v2.z}, {v3.x, d.z, v3.z}};
    std::vector<std::vector<double>> Ad3 = {
        {v1.x, v1.y, d.x}, {v2.x, v2.y, d.y}, {v3.x, v3.y, d.z}};

    double detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                  A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                  A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    double invA[3][3];
    invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / detA;
    invA[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / detA;
    invA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detA;

    invA[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / detA;
    invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / detA;
    invA[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / detA;

    invA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / detA;
    invA[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / detA;
    invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / detA;

    x = invA[0][0] * d.x + invA[0][1] * d.y + invA[0][2] * d.z;
    y = invA[1][0] * d.x + invA[1][1] * d.y + invA[1][2] * d.z;
    z = invA[2][0] * d.x + invA[2][1] * d.y + invA[2][2] * d.z;
}

void MathUtils::calculateLocalCoordinates(Point physical_point, Point node_0,
                                          Point node_1, Point node_2,
                                          Point node_3, Point vertex,
                                          double& ksi_parameter,
                                          double& nu_parameter,
                                          double& tetta_parameter) {
    Point const center_of_base = {
        (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
        (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
        (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    Point const vertice_0_minus_1 = {node_1.x - node_0.x, node_1.y - node_0.y,
                                     node_1.z - node_0.z};

    Point const vertice_0_minus_4 = {node_2.x - node_0.x, node_2.y - node_0.y,
                                     node_2.z - node_0.z};

    Point const vertice_0_minus_point = {physical_point.x - node_0.x,
                                         physical_point.y - node_0.y,
                                         physical_point.z - node_0.z};

    double const norm_of_vertice_0_minus_1 =
        vertice_0_minus_1.x * vertice_0_minus_1.x +
        vertice_0_minus_1.y * vertice_0_minus_1.y +
        vertice_0_minus_1.z * vertice_0_minus_1.z;

    double const norm_of_vertice_0_minus_4 =
        vertice_0_minus_4.x * vertice_0_minus_4.x +
        vertice_0_minus_4.y * vertice_0_minus_4.y +
        vertice_0_minus_4.z * vertice_0_minus_4.z;

    ksi_parameter = (vertice_0_minus_point.x * vertice_0_minus_1.x +
                     vertice_0_minus_point.y * vertice_0_minus_1.y +
                     vertice_0_minus_point.z * vertice_0_minus_1.z) /
                    norm_of_vertice_0_minus_1;

    nu_parameter = (vertice_0_minus_point.x * vertice_0_minus_4.x +
                    vertice_0_minus_point.y * vertice_0_minus_4.y +
                    vertice_0_minus_point.z * vertice_0_minus_4.z) /
                   norm_of_vertice_0_minus_4;

    Point const vetice_vertex_minus_center = {(vertex.x - center_of_base.x),
                                              (vertex.y - center_of_base.y),
                                              (vertex.z - center_of_base.z)};

    double const norm_of_vetice_vertex_minus_center =
        vetice_vertex_minus_center.x * vetice_vertex_minus_center.x +
        vetice_vertex_minus_center.y * vetice_vertex_minus_center.y +
        vetice_vertex_minus_center.z * vetice_vertex_minus_center.z;

    Point const vetice_point_minus_center = {
        physical_point.x - center_of_base.x,
        physical_point.y - center_of_base.y,
        physical_point.z - center_of_base.z};

    tetta_parameter =
        (vetice_point_minus_center.x * vetice_vertex_minus_center.x +
         vetice_point_minus_center.y * vetice_vertex_minus_center.y +
         vetice_point_minus_center.z * vetice_vertex_minus_center.z) /
        norm_of_vetice_vertex_minus_center;
}

double MathUtils::getBasisFunction(
    const std::vector<Point>& nodes, Point point,
    const std::vector<int>& node_indexes, BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_elements_type, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    Point node_0, node_1, node_2, node_3, vertex;

    getMainPyramideNodes(node_indexes, nodes, basis_functions_type,
                         basis_functions_elements_type, node_0, node_1, node_2,
                         node_3, vertex);

    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    return getLinearBasisFunction(
                        point, node_0, node_1, node_2, node_3, vertex,
                        number_basis_function,
                        need_coordinate_transform_to_local);
                    break;
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getQuadraticBasisFunction(
                        point, node_0, node_1, node_2, node_3, vertex,
                        number_basis_function,
                        need_coordinate_transform_to_local);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    return getCubicBasisFunction(
                        point, node_0, node_1, node_2, node_3, vertex,
                        number_basis_function,
                        need_coordinate_transform_to_local);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        case BASIS_TYPE::Hierarhical:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getHierarhicalQuadraticBasisFunction(
                        point, node_0, node_1, node_2, node_3, vertex,
                        number_basis_function,
                        need_coordinate_transform_to_local);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

double MathUtils::getDerivativeBasisFunction(
    const std::vector<Point>& nodes, Point point,
    const std::vector<int>& node_indexes, BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_elements_type,
    int number_of_first_basis_function, int number_of_second_basis_function,
    bool need_coordinate_transform_to_local) {
    Point node_0, node_1, node_2, node_3, vertex;

    getMainPyramideNodes(node_indexes, nodes, basis_functions_type,
                         basis_functions_elements_type, node_0, node_1, node_2,
                         node_3, vertex);

    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    return getDerivativeLinearBasisFunction(
                        number_of_first_basis_function,
                        number_of_second_basis_function, point, node_0, node_1,
                        node_2, node_3, vertex,
                        need_coordinate_transform_to_local);
                    break;
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getDerivativeQuadraticBasisFunction(
                        number_of_first_basis_function,
                        number_of_second_basis_function, point, node_0, node_1,
                        node_2, node_3, vertex,
                        need_coordinate_transform_to_local);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    return getDerivativeCubicBasisFunction(
                        number_of_first_basis_function,
                        number_of_second_basis_function, point, node_0, node_1,
                        node_2, node_3, vertex,
                        need_coordinate_transform_to_local);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        case BASIS_TYPE::Hierarhical:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getDerivativeHierarhicalQuadraticBasisFunction(
                        number_of_first_basis_function,
                        number_of_second_basis_function, point, node_0, node_1,
                        node_2, node_3, vertex,
                        need_coordinate_transform_to_local);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    return getDerivativeHierarhicalCubicBasisFunction(
                        number_of_first_basis_function,
                        number_of_second_basis_function, point, node_0, node_1,
                        node_2, node_3, vertex,
                        need_coordinate_transform_to_local);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

double MathUtils::getLinearBasisFunction(
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    double xi = 0.0;
    double eta = 0.0;
    double theta = 0.0;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  xi, eta, theta);
    } else {
        xi = point.x;
        eta = point.y;
        theta = point.z;
    }

    double xi_0 = 1 - xi;
    double xi_1 = xi;

    double eta_0 = 1 - eta;
    double eta_1 = eta;

    double theta_0 = 1 - theta;
    double theta_1 = theta;

    switch (number_basis_function) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_0 * eta_1 * theta_0;
            break;
        case 3:
            return xi_1 * eta_1 * theta_0;
            break;
        case 4:
            return theta_1;
            break;
        default:
            return 0.0;
            break;
    }
}

double MathUtils::getQuadraticBasisFunction(
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    double xi = 0.0;
    double eta = 0.0;
    double theta = 0.0;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  xi, eta, theta);
    } else {
        xi = point.x;
        eta = point.y;
        theta = point.z;
    }

    double xi_0 = (1 - 2 * xi) * (1 - xi);
    double xi_1 = 4 * xi * (1 - xi);
    double xi_2 = -xi * (1 - 2 * xi);

    double eta_0 = (1 - 2 * eta) * (1 - eta);
    double eta_1 = 4 * eta * (1 - eta);
    double eta_2 = -eta * (1 - 2 * eta);

    double theta_0 = (1 - 2 * theta) * (1 - theta);
    double theta_1 = 4 * theta * (1 - theta);
    double theta_2 = -theta * (1 - 2 * theta);

    switch (number_basis_function) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_2 * eta_0 * theta_0;
            break;
        case 3:
            return xi_0 * eta_1 * theta_0;
            break;
        case 4:
            return xi_1 * eta_1 * theta_0;
            break;
        case 5:
            return xi_2 * eta_1 * theta_0;
            break;
        case 6:
            return xi_0 * eta_2 * theta_0;
            break;
        case 7:
            return xi_1 * eta_2 * theta_0;
            break;
        case 8:
            return xi_2 * eta_2 * theta_0;
            break;
        case 9:
            return 4 * (0.75 - xi) * (0.75 - eta) * theta_1;
            break;
        case 10:
            return -4 * (0.25 - xi) * (0.75 - eta) * theta_1;
            break;
        case 11:
            return -4 * (0.75 - xi) * (0.25 - eta) * theta_1;
            break;
        case 12:
            return 4 * (0.25 - xi) * (0.25 - eta) * theta_1;
            break;
        case 13:
            return theta_2;
            break;
        default:
            return 0.0;
            break;
    }
}

double MathUtils::getCubicBasisFunction(
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    double xi = 0.0;
    double eta = 0.0;
    double theta = 0.0;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  xi, eta, theta);
    } else {
        xi = point.x;
        eta = point.y;
        theta = point.z;
    }

    double xi_0 = 1.0 / 2.0 * (1 - 3 * xi) * (2 - 3 * xi) * (1 - xi);
    double xi_1 = 9.0 / 2.0 * xi * (2 - 3 * xi) * (1 - xi);
    double xi_2 = -9.0 / 2.0 * xi * (1 - 3 * xi) * (1 - xi);
    double xi_3 = 1.0 / 2.0 * xi * (1 - 3 * xi) * (2 - 3 * xi);

    double eta_0 = 1.0 / 2.0 * (1 - 3 * eta) * (2 - 3 * eta) * (1 - eta);
    double eta_1 = 9.0 / 2.0 * eta * (2 - 3 * eta) * (1 - eta);
    double eta_2 = -9.0 / 2.0 * eta * (1 - 3 * eta) * (1 - eta);
    double eta_3 = 1.0 / 2.0 * eta * (1 - 3 * eta) * (2 - 3 * eta);

    double theta_0 =
        1.0 / 2.0 * (1 - 3 * theta) * (2 - 3 * theta) * (1 - theta);
    double theta_1 = 9.0 / 2.0 * theta * (2 - 3 * theta) * (1 - theta);
    double theta_2 = -9.0 / 2.0 * theta * (1 - 3 * theta) * (1 - theta);
    double theta_3 = 1.0 / 2.0 * theta * (1 - 3 * theta) * (2 - 3 * theta);

    switch (number_basis_function) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_2 * eta_0 * theta_0;
            break;
        case 3:
            return xi_3 * eta_0 * theta_0;
            break;
        case 4:
            return xi_0 * eta_1 * theta_0;
            break;
        case 5:
            return xi_1 * eta_1 * theta_0;
            break;
        case 6:
            return xi_2 * eta_1 * theta_0;
            break;
        case 7:
            return xi_3 * eta_1 * theta_0;
            break;
        case 8:
            return xi_0 * eta_2 * theta_0;
            break;
        case 9:
            return xi_1 * eta_2 * theta_0;
            break;
        case 10:
            return xi_2 * eta_2 * theta_0;
            break;
        case 11:
            return xi_3 * eta_2 * theta_0;
            break;
        case 12:
            return xi_0 * eta_3 * theta_0;
            break;
        case 13:
            return xi_1 * eta_3 * theta_0;
            break;
        case 14:
            return xi_2 * eta_3 * theta_0;
            break;
        case 15:
            return xi_3 * eta_3 * theta_0;
            break;
        case 16:
            return 1.0 / 16.0 * (5 - 6 * xi) * (5 - 6 * eta) * theta_1;
            break;
        case 17:
            return -1.0 / 16.0 * (1 - 6 * xi) * (5 - 6 * eta) * theta_1;
            break;
        case 18:
            return -1.0 / 16.0 * (5 - 6 * xi) * (1 - 6 * eta) * theta_1;
            break;
        case 19:
            return 1.0 / 16.0 * (1 - 6 * xi) * (1 - 6 * eta) * theta_1;
            break;
        case 20:
            return (2 - 3 * xi) * (2 - 3 * eta) * theta_2;
            break;
        case 21:
            return -(1 - 3 * xi) * (2 - 3 * eta) * theta_2;
            break;
        case 22:
            return -(2 - 3 * xi) * (1 - 3 * eta) * theta_2;
            break;
        case 23:
            return (1 - 3 * xi) * (1 - 3 * eta) * theta_2;
            break;
        case 24:
            return theta_3;
            break;
        default:
            return 0.0;
            break;
    }
}

double MathUtils::getHierarhicalQuadraticBasisFunction(
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    double xi = 0.0;
    double eta = 0.0;
    double theta = 0.0;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  xi, eta, theta);
    } else {
        xi = point.x;
        eta = point.y;
        theta = point.z;
    }

    double xi_0 = 1 - xi;
    double xi_1 = xi;
    double xi_2 = 1.0 / 8.0 * ((2 * xi - 1) * (2 * xi - 1) - 1);

    double eta_0 = 1 - eta;
    double eta_1 = eta;
    double eta_2 = 1.0 / 8.0 * ((2 * eta - 1) * (2 * eta - 1) - 1);

    double theta_0 = 1 - theta;
    double theta_1 = theta;
    double theta_2 = 1.0 / 8.0 * ((2 * theta - 1) * (2 * theta - 1) - 1);

    switch (number_basis_function) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_2 * eta_0 * theta_0;
            break;
        case 2:
            return xi_1 * eta_0 * theta_0;
            break;
        case 3:
            return xi_0 * eta_2 * theta_0;
            break;
        case 4:
            return xi_2 * eta_2 * theta_0;
            break;
        case 5:
            return xi_1 * eta_2 * theta_0;
            break;
        case 6:
            return xi_0 * eta_1 * theta_0;
            break;
        case 7:
            return xi_2 * eta_1 * theta_0;
            break;
        case 8:
            return xi_1 * eta_1 * theta_0;
            break;
        case 9:
            return xi_0 * eta_0 * theta_2;
            break;
        case 10:
            return xi_1 * eta_0 * theta_2;
            break;
        case 11:
            return xi_0 * eta_1 * theta_2;
            break;
        case 12:
            return xi_1 * eta_1 * theta_2;
            break;
        case 13:
            return theta_1;
            break;
        default:
            return 0.0;
            break;
    }
}

double MathUtils::getHierarhicalCubicBasisFunction(
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, int number_basis_function,
    bool need_coordinate_transform_to_local) {
    double xi = 0.0;
    double eta = 0.0;
    double theta = 0.0;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  xi, eta, theta);
    } else {
        xi = point.x;
        eta = point.y;
        theta = point.z;
    }

    const double xi_0 = 1 - xi;
    const double xi_1 = xi;
    const double xi_2 = ((2 * xi - 1) * (2 * xi - 1) - 1) / 8.0;

    const double eta_0 = 1 - eta;
    const double eta_1 = eta;
    const double eta_2 = ((2 * eta - 1) * (2 * eta - 1) - 1) / 8.0;

    const double theta_0 = 1 - theta;
    const double theta_1 = theta;
    const double theta_2 = ((2 * theta - 1) * (2 * theta - 1) - 1) / 8.0;

    switch (number_basis_function) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_2 * eta_0 * theta_0;
            break;
        case 2:
            return xi_1 * eta_0 * theta_0;
            break;
        case 3:
            return xi_0 * eta_2 * theta_0;
            break;
        case 4:
            return xi_1 * eta_2 * theta_0;
            break;
        case 5:
            return xi_0 * eta_1 * theta_0;
            break;
        case 6:
            return xi_2 * eta_1 * theta_0;
            break;
        case 7:
            return xi_1 * eta_1 * theta_0;
            break;
        case 8:
            return xi_0 * eta_0 * theta_2;
            break;
        case 9:
            return xi_1 * eta_0 * theta_2;
            break;
        case 10:
            return xi_0 * eta_1 * theta_2;
            break;
        case 11:
            return xi_1 * eta_1 * theta_2;
            break;
        case 12:
            return theta_1;
            break;
        default:
            return 0.0;
            break;
    }
}

double MathUtils::getDerivativeLinearBasisFunction(
    int number_of_first_basis_function, int number_of_second_basis_function,
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, bool need_coordinate_transform_to_local) {
    Point basis_point;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  basis_point.x, basis_point.y, basis_point.z);
    } else {
        basis_point.x = point.x;
        basis_point.y = point.y;
        basis_point.z = point.z;
    }

    double first_d_psi_d_xi = 0.0;
    double first_d_psi_d_eta = 0.0;
    double first_d_psi_d_theta = 0.0;

    dual xi = basis_point.x;
    dual eta = basis_point.y;
    dual theta = basis_point.z;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Linear,
                      number_of_first_basis_function, xi, eta, theta,
                      first_d_psi_d_xi, first_d_psi_d_eta, first_d_psi_d_theta);

    double second_d_psi_d_xi = 0.0;
    double second_d_psi_d_eta = 0.0;
    double second_d_psi_d_theta = 0.0;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Linear,
                      number_of_second_basis_function, xi, eta, theta,
                      second_d_psi_d_xi, second_d_psi_d_eta,
                      second_d_psi_d_theta);

    std::vector<std::vector<double>> jacobi_matrix(3, std::vector<double>(3));
    std::vector<std::vector<double>> transposed_jacobi_matrix(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> multiplied_matrixes_jacobi_and_transposed(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> inverted_multiplied_matrixes(
        3, std::vector<double>(3));

    calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex,
                          jacobi_matrix);

    const double determinant = calculateJacobian(jacobi_matrix);

    transponateMatrix(jacobi_matrix, transposed_jacobi_matrix);

    multiplyMatrixByMatrix(jacobi_matrix, transposed_jacobi_matrix,
                           multiplied_matrixes_jacobi_and_transposed);

    calculateInvertedMatrix(multiplied_matrixes_jacobi_and_transposed,
                            inverted_multiplied_matrixes);

    const Point first = {
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][0] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][0] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][0],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][1] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][1] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][1],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][2] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][2] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][2]};

    return (first.x * second_d_psi_d_xi + first.y * second_d_psi_d_eta +
            first.z * second_d_psi_d_theta) *
           std::abs(determinant);
}

double MathUtils::getDerivativeQuadraticBasisFunction(
    int number_of_first_basis_function, int number_of_second_basis_function,
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, bool need_coordinate_transform_to_local) {
    Point basis_point;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  basis_point.x, basis_point.y, basis_point.z);
    } else {
        basis_point.x = point.x;
        basis_point.y = point.y;
        basis_point.z = point.z;
    }

    double first_d_psi_d_xi = 0.0;
    double first_d_psi_d_eta = 0.0;
    double first_d_psi_d_theta = 0.0;

    dual xi = basis_point.x;
    dual eta = basis_point.y;
    dual theta = basis_point.z;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Quadratic,
                      number_of_first_basis_function, xi, eta, theta,
                      first_d_psi_d_xi, first_d_psi_d_eta, first_d_psi_d_theta);

    double second_d_psi_d_xi = 0.0;
    double second_d_psi_d_eta = 0.0;
    double second_d_psi_d_theta = 0.0;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Quadratic,
                      number_of_second_basis_function, xi, eta, theta,
                      second_d_psi_d_xi, second_d_psi_d_eta,
                      second_d_psi_d_theta);

    std::vector<std::vector<double>> jacobi_matrix(3, std::vector<double>(3));
    std::vector<std::vector<double>> transposed_jacobi_matrix(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> multiplied_matrixes_jacobi_and_transposed(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> inverted_multiplied_matrixes(
        3, std::vector<double>(3));

    calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex,
                          jacobi_matrix);

    const double determinant = calculateJacobian(jacobi_matrix);

    transponateMatrix(jacobi_matrix, transposed_jacobi_matrix);

    multiplyMatrixByMatrix(jacobi_matrix, transposed_jacobi_matrix,
                           multiplied_matrixes_jacobi_and_transposed);

    calculateInvertedMatrix(multiplied_matrixes_jacobi_and_transposed,
                            inverted_multiplied_matrixes);

    const Point first = {
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][0] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][0] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][0],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][1] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][1] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][1],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][2] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][2] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][2]};

    return (first.x * second_d_psi_d_xi + first.y * second_d_psi_d_eta +
            first.z * second_d_psi_d_theta) *
           std::abs(determinant);
}

double MathUtils::getDerivativeCubicBasisFunction(
    int number_of_first_basis_function, int number_of_second_basis_function,
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, bool need_coordinate_transform_to_local) {
    Point basis_point;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  basis_point.x, basis_point.y, basis_point.z);
    } else {
        basis_point.x = point.x;
        basis_point.y = point.y;
        basis_point.z = point.z;
    }

    double first_d_psi_d_xi = 0.0;
    double first_d_psi_d_eta = 0.0;
    double first_d_psi_d_theta = 0.0;

    dual xi = basis_point.x;
    dual eta = basis_point.y;
    dual theta = basis_point.z;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Cubic,
                      number_of_first_basis_function, xi, eta, theta,
                      first_d_psi_d_xi, first_d_psi_d_eta, first_d_psi_d_theta);

    double second_d_psi_d_xi = 0.0;
    double second_d_psi_d_eta = 0.0;
    double second_d_psi_d_theta = 0.0;

    computeDerivative(BASIS_TYPE::Lagrange, BASIS_ELEMENT_TYPE::Cubic,
                      number_of_second_basis_function, xi, eta, theta,
                      second_d_psi_d_xi, second_d_psi_d_eta,
                      second_d_psi_d_theta);

    std::vector<std::vector<double>> jacobi_matrix(3, std::vector<double>(3));
    std::vector<std::vector<double>> transposed_jacobi_matrix(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> multiplied_matrixes_jacobi_and_transposed(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> inverted_multiplied_matrixes(
        3, std::vector<double>(3));

    calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex,
                          jacobi_matrix);

    const double determinant = calculateJacobian(jacobi_matrix);

    transponateMatrix(jacobi_matrix, transposed_jacobi_matrix);

    multiplyMatrixByMatrix(jacobi_matrix, transposed_jacobi_matrix,
                           multiplied_matrixes_jacobi_and_transposed);

    calculateInvertedMatrix(multiplied_matrixes_jacobi_and_transposed,
                            inverted_multiplied_matrixes);

    const Point first = {
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][0] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][0] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][0],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][1] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][1] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][1],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][2] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][2] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][2]};

    return (first.x * second_d_psi_d_xi + first.y * second_d_psi_d_eta +
            first.z * second_d_psi_d_theta) *
           std::abs(determinant);
}

double MathUtils::getDerivativeHierarhicalQuadraticBasisFunction(
    int number_of_first_basis_function, int number_of_second_basis_function,
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, bool need_coordinate_transform_to_local) {
    Point basis_point;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  basis_point.x, basis_point.y, basis_point.z);
    } else {
        basis_point.x = point.x;
        basis_point.y = point.y;
        basis_point.z = point.z;
    }

    double first_d_psi_d_xi = 0.0;
    double first_d_psi_d_eta = 0.0;
    double first_d_psi_d_theta = 0.0;

    dual xi = basis_point.x;
    dual eta = basis_point.y;
    dual theta = basis_point.z;

    computeDerivative(BASIS_TYPE::Hierarhical, BASIS_ELEMENT_TYPE::Quadratic,
                      number_of_first_basis_function, xi, eta, theta,
                      first_d_psi_d_xi, first_d_psi_d_eta, first_d_psi_d_theta);

    double second_d_psi_d_xi = 0.0;
    double second_d_psi_d_eta = 0.0;
    double second_d_psi_d_theta = 0.0;

    computeDerivative(BASIS_TYPE::Hierarhical, BASIS_ELEMENT_TYPE::Quadratic,
                      number_of_second_basis_function, xi, eta, theta,
                      second_d_psi_d_xi, second_d_psi_d_eta,
                      second_d_psi_d_theta);

    std::vector<std::vector<double>> jacobi_matrix(3, std::vector<double>(3));
    std::vector<std::vector<double>> transposed_jacobi_matrix(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> multiplied_matrixes_jacobi_and_transposed(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> inverted_multiplied_matrixes(
        3, std::vector<double>(3));

    calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex,
                          jacobi_matrix);

    const double determinant = calculateJacobian(jacobi_matrix);

    transponateMatrix(jacobi_matrix, transposed_jacobi_matrix);

    multiplyMatrixByMatrix(jacobi_matrix, transposed_jacobi_matrix,
                           multiplied_matrixes_jacobi_and_transposed);

    calculateInvertedMatrix(multiplied_matrixes_jacobi_and_transposed,
                            inverted_multiplied_matrixes);

    const Point first = {
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][0] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][0] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][0],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][1] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][1] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][1],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][2] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][2] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][2]};

    return (first.x * second_d_psi_d_xi + first.y * second_d_psi_d_eta +
            first.z * second_d_psi_d_theta) *
           std::abs(determinant);
}

double MathUtils::getDerivativeHierarhicalCubicBasisFunction(
    int number_of_first_basis_function, int number_of_second_basis_function,
    Point point, Point node_0, Point node_1, Point node_2, Point node_3,
    Point vertex, bool need_coordinate_transform_to_local) {
    Point basis_point;

    if (need_coordinate_transform_to_local) {
        calculateLocalCoordinates(point, node_0, node_1, node_2, node_3, vertex,
                                  basis_point.x, basis_point.y, basis_point.z);
    } else {
        basis_point.x = point.x;
        basis_point.y = point.y;
        basis_point.z = point.z;
    }

    double first_d_psi_d_xi = 0.0;
    double first_d_psi_d_eta = 0.0;
    double first_d_psi_d_theta = 0.0;

    dual xi = basis_point.x;
    dual eta = basis_point.y;
    dual theta = basis_point.z;

    computeDerivative(BASIS_TYPE::Hierarhical, BASIS_ELEMENT_TYPE::Cubic,
                      number_of_first_basis_function, xi, eta, theta,
                      first_d_psi_d_xi, first_d_psi_d_eta, first_d_psi_d_theta);

    double second_d_psi_d_xi = 0.0;
    double second_d_psi_d_eta = 0.0;
    double second_d_psi_d_theta = 0.0;

    computeDerivative(BASIS_TYPE::Hierarhical, BASIS_ELEMENT_TYPE::Cubic,
                      number_of_second_basis_function, xi, eta, theta,
                      second_d_psi_d_xi, second_d_psi_d_eta,
                      second_d_psi_d_theta);

    std::vector<std::vector<double>> jacobi_matrix(3, std::vector<double>(3));
    std::vector<std::vector<double>> transposed_jacobi_matrix(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> multiplied_matrixes_jacobi_and_transposed(
        3, std::vector<double>(3));
    std::vector<std::vector<double>> inverted_multiplied_matrixes(
        3, std::vector<double>(3));

    calculateJacobiMatrix(node_0, node_1, node_2, node_3, vertex,
                          jacobi_matrix);

    const double determinant = calculateJacobian(jacobi_matrix);

    transponateMatrix(jacobi_matrix, transposed_jacobi_matrix);

    multiplyMatrixByMatrix(jacobi_matrix, transposed_jacobi_matrix,
                           multiplied_matrixes_jacobi_and_transposed);

    calculateInvertedMatrix(multiplied_matrixes_jacobi_and_transposed,
                            inverted_multiplied_matrixes);

    const Point first = {
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][0] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][0] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][0],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][1] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][1] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][1],
        first_d_psi_d_xi * inverted_multiplied_matrixes[0][2] +
            first_d_psi_d_eta * inverted_multiplied_matrixes[1][2] +
            first_d_psi_d_theta * inverted_multiplied_matrixes[2][2]};

    return (first.x * second_d_psi_d_xi + first.y * second_d_psi_d_eta +
            first.z * second_d_psi_d_theta) *
           std::abs(determinant);
}

void MathUtils::getMainPyramideNodes(
    const std::vector<int>& node_indexes, const std::vector<Point>& nodes,
    BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_elements_type, Point& node_0,
    Point& node_1, Point& node_2, Point& node_3, Point& node_4) {
    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    node_0 = nodes[node_indexes.at(0)];
                    node_1 = nodes[node_indexes.at(1)];
                    node_2 = nodes[node_indexes.at(2)];
                    node_3 = nodes[node_indexes.at(3)];
                    node_4 = nodes[nodes.size() - 1];
                    break;
                case BASIS_ELEMENT_TYPE::Quadratic:
                    node_0 = nodes[node_indexes.at(0)];
                    node_1 = nodes[node_indexes.at(2)];
                    node_2 = nodes[node_indexes.at(6)];
                    node_3 = nodes[node_indexes.at(8)];
                    node_4 = nodes[nodes.size() - 1];
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    node_0 = nodes[node_indexes.at(0)];
                    node_1 = nodes[node_indexes.at(3)];
                    node_2 = nodes[node_indexes.at(12)];
                    node_3 = nodes[node_indexes.at(15)];
                    node_4 = nodes[nodes.size() - 1];
                    break;
                default:
                    break;
            }
            break;
        case BASIS_TYPE::Hierarhical:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Quadratic:
                    node_0 = nodes[node_indexes.at(0)];
                    node_1 = nodes[node_indexes.at(2)];
                    node_2 = nodes[node_indexes.at(6)];
                    node_3 = nodes[node_indexes.at(8)];
                    node_4 = nodes[nodes.size() - 1];
                    break;
            }
            break;
        default:
            break;
    }
}

double MathUtils::calculateJacobian(std::vector<std::vector<double>>& J) {
    return J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2]) -
           J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2]) +
           J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
}

void MathUtils::calculateJacobiMatrix(Point node_0, Point node_1, Point node_2,
                                      Point node_3, Point node_4,
                                      std::vector<std::vector<double>>& J) {
    Point node_c = {(node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
                    (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
                    (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    double D10 = pow(node_1.x - node_0.x, 2) + pow(node_1.y - node_0.y, 2) +
                 pow(node_1.z - node_0.z, 2);

    double D20 = pow(node_2.x - node_0.x, 2) + pow(node_2.y - node_0.y, 2) +
                 pow(node_2.z - node_0.z, 2);

    double D4c = pow(node_4.x - node_c.x, 2) + pow(node_4.y - node_c.y, 2) +
                 pow(node_4.z - node_c.z, 2);

    J[0][0] = D10 / (node_1.x - node_0.x);  // d_x_d_ksi
    J[0][1] = D10 / (node_1.y - node_0.y);  // d_y_d_ksi
    J[0][2] = D10 / (node_1.z - node_0.z);  // d_z_d_ksi

    J[1][0] = D20 / (node_2.x - node_0.x);  // d_x_d_nu
    J[1][1] = D20 / (node_2.y - node_0.y);  // d_y_d_nu
    J[1][2] = D20 / (node_2.z - node_0.z);  // d_z_d_nu

    J[2][0] = D4c / (node_4.x - node_c.x);  // d_x_d_tetta
    J[2][1] = D4c / (node_4.y - node_c.y);  // d_y_d_tetta
    J[2][2] = D4c / (node_4.z - node_c.z);  // d_z_d_tetta

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (J[i][j] == INFINITY) {
                J[i][j] = 0.0;
            }
        }
    }
}
void MathUtils::transponateMatrix(std::vector<std::vector<double>>& J,
                                  std::vector<std::vector<double>>& JT) {
    JT[0][0] = J[0][0];
    JT[0][1] = J[1][0];
    JT[0][2] = J[2][0];

    JT[1][0] = J[0][1];
    JT[1][1] = J[1][1];
    JT[1][2] = J[2][1];

    JT[2][0] = J[0][2];
    JT[2][1] = J[1][2];
    JT[2][2] = J[2][2];
}

void MathUtils::calculateInvertedMatrix(
    std::vector<std::vector<double>>& J,
    std::vector<std::vector<double>>& JInverted) {
    const double determinant = calculateJacobian(J);

    std::vector<std::vector<double>> J_Cofactor(3, std::vector<double>(3));

    J_Cofactor[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
    J_Cofactor[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]);
    J_Cofactor[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]);

    J_Cofactor[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]);
    J_Cofactor[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
    J_Cofactor[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]);

    J_Cofactor[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    J_Cofactor[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]);
    J_Cofactor[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);

    JInverted[0][0] = J_Cofactor[0][0] / determinant;
    JInverted[0][1] = J_Cofactor[1][0] / determinant;
    JInverted[0][2] = J_Cofactor[2][0] / determinant;

    JInverted[1][0] = J_Cofactor[0][1] / determinant;
    JInverted[1][1] = J_Cofactor[1][1] / determinant;
    JInverted[1][2] = J_Cofactor[2][1] / determinant;

    JInverted[2][0] = J_Cofactor[0][2] / determinant;
    JInverted[2][1] = J_Cofactor[1][2] / determinant;
    JInverted[2][2] = J_Cofactor[2][2] / determinant;
}

void MathUtils::multiplyMatrixByMatrix(
    std::vector<std::vector<double>>& first_matrix,
    std::vector<std::vector<double>>& second_matrix,
    std::vector<std::vector<double>>& result) {
    for (int i = 0; i < first_matrix.size(); ++i) {
        for (int j = 0; j < first_matrix.size(); ++j) {
            result[i][j] = 0;

            for (int k = 0; k < first_matrix.size(); ++k) {
                result[i][j] += first_matrix[i][k] * second_matrix[k][j];
            }
        }
    }
}