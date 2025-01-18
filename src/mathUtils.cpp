#include "../include/mathUtils.h"

#include <cmath>
#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"

const double NUMBERS_EQUAL_EPSILON = 1e-15;

const double TETRAHEDRON_VOLUME_DIVISOR = 6.0;

bool MathUtils::isQudraticNodeIllegal(
    const std::tuple<int, double, double, double>& qudratic_node_parameters) {
    const int grid_z_index = std::get<0>(qudratic_node_parameters);
    const double grid_x_val = std::get<1>(qudratic_node_parameters);
    const double grid_y_val = std::get<2>(qudratic_node_parameters);
    const double center_of_line = std::get<3>(qudratic_node_parameters);

    if (grid_z_index == 0 || grid_z_index == 4) {
        return grid_x_val == grid_y_val && grid_y_val == center_of_line;
    }

    if (grid_z_index == 2) {
        return grid_x_val == center_of_line || grid_y_val == center_of_line;
    }

    return false;
}

bool MathUtils::isCubicNodeIllegal(
    const std::tuple<int, int, int>& cubic_node_parameters) {
    const int grid_z_index = std::get<0>(cubic_node_parameters);
    const int grid_y_index = std::get<1>(cubic_node_parameters);
    const int grid_x_index = std::get<2>(cubic_node_parameters);

    const int last_grid_x_y_z_index = 5;

    if (grid_z_index == 0 || grid_z_index == last_grid_x_y_z_index) {
        const bool first_condition = (grid_y_index != 1 && grid_y_index != 4 &&
                                      grid_x_index != 1 && grid_x_index != 4);

        const bool second_condition =
            (grid_y_index != 2 || grid_x_index != 2) &&
            (grid_y_index != 2 || grid_x_index != 3) &&
            (grid_y_index != 3 || grid_x_index != 2) &&
            (grid_y_index != 3 || grid_x_index != 3);

        return !(first_condition && second_condition);
    }

    if (grid_z_index == 2 || grid_z_index == 3) {
        return (grid_y_index != 0 && grid_y_index != last_grid_x_y_z_index) ||
               (grid_x_index != 0 && grid_x_index != last_grid_x_y_z_index);
    }

    return true;
}

double MathUtils::calculateF(Point point) {
    const double const_right_part = 5.0;
    const double quadratic_div_grad = -6.0;

    // return const_right_part + 0 * (point.x + point.y + point.z);
    return point.x + point.y + point.z;
    // return quadratic_div_grad + pow(point.x, 2) + pow(point.y, 2) +
    //        pow(point.z, 2);
}

double MathUtils::calculateU(Point point) {
    const double const_right_part = 5.0;

    // return const_right_part + 0 * (point.x + point.y + point.z);
    return point.x + point.y + point.z;
    // return pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2);
    // return pow(point.x, 3) + pow(point.y, 3) + pow(point.z, 3);
    // return sin(point.x * point.y * point.z);
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
            return isNumbersEqual(nodes[node_indexes[0]].x, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[1]].x, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[2]].x, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[3]].x, average,
                                  NUMBERS_EQUAL_EPSILON);
        case 'y':
            return isNumbersEqual(nodes[node_indexes[0]].y, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[1]].y, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[2]].y, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[3]].y, average,
                                  NUMBERS_EQUAL_EPSILON);
        case 'z':
            return isNumbersEqual(nodes[node_indexes[0]].z, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[1]].z, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[2]].z, average,
                                  NUMBERS_EQUAL_EPSILON) &&
                   isNumbersEqual(nodes[node_indexes[3]].z, average,
                                  NUMBERS_EQUAL_EPSILON);
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

bool MathUtils::isPointInPyramid(std::vector<Point>& nodes, const Point& point,
                                 const Element& element) {
    const auto& node_indexes = element.getNodeIndexes();

    Point const first_point_of_pyramid = nodes[node_indexes.at(0)];
    Point const second_point_of_pyramid = nodes[node_indexes.at(1)];
    Point const third_point_of_pyramid = nodes[node_indexes.at(2)];
    Point const fourth_point_of_pyramid = nodes[node_indexes.at(3)];
    Point const fifth_point_of_pyramid = nodes[node_indexes.at(4)];

    if (isPointInsideTetrahedron(
            point, first_point_of_pyramid, second_point_of_pyramid,
            third_point_of_pyramid, fifth_point_of_pyramid)) {
        return true;
    }

    if (isPointInsideTetrahedron(
            point, first_point_of_pyramid, third_point_of_pyramid,
            fourth_point_of_pyramid, fifth_point_of_pyramid)) {
        return true;
    }

    return false;
}

double MathUtils::calculateTetrahedronVolume(Point first_tetrahedron_point,
                                             Point second_tetrahedron_point,
                                             Point third_tetrahedron_point,
                                             Point fourth_tetrahedron_point) {
    Point const vertice_1 = {
        second_tetrahedron_point.x - first_tetrahedron_point.x,
        second_tetrahedron_point.y - first_tetrahedron_point.y,
        second_tetrahedron_point.z - first_tetrahedron_point.z};

    Point const vertice_2 = {
        third_tetrahedron_point.x - first_tetrahedron_point.x,
        third_tetrahedron_point.y - first_tetrahedron_point.y,
        third_tetrahedron_point.z - first_tetrahedron_point.z};

    Point const vertice_3 = {
        fourth_tetrahedron_point.x - first_tetrahedron_point.x,
        fourth_tetrahedron_point.y - first_tetrahedron_point.y,
        fourth_tetrahedron_point.z - first_tetrahedron_point.z};

    Point const cross = calculateVectorCrossProduct(vertice_2, vertice_3);

    return std::abs(vertice_1.x * cross.x + vertice_1.y * cross.y +
                    vertice_1.z * cross.z) /
           TETRAHEDRON_VOLUME_DIVISOR;
}

bool MathUtils::isPointInsideTetrahedron(Point point,
                                         Point first_tetrahedron_point,
                                         Point second_tetrahedron_point,
                                         Point third_tetrahedron_point,
                                         Point fourth_tetrahedron_point) {
    const double epsilon = 0.00001;

    double const total_tetrahedron_volume = calculateTetrahedronVolume(
        first_tetrahedron_point, second_tetrahedron_point,
        third_tetrahedron_point, fourth_tetrahedron_point);

    double const first_tetrahedron_volume = calculateTetrahedronVolume(
        point, second_tetrahedron_point, third_tetrahedron_point,
        fourth_tetrahedron_point);

    double const second_tetrahedron_volume = calculateTetrahedronVolume(
        first_tetrahedron_point, point, third_tetrahedron_point,
        fourth_tetrahedron_point);

    double const third_tetrahedron_volume = calculateTetrahedronVolume(
        first_tetrahedron_point, second_tetrahedron_point, point,
        fourth_tetrahedron_point);

    double const fourth_tetrahedron_volume = calculateTetrahedronVolume(
        first_tetrahedron_point, second_tetrahedron_point,
        third_tetrahedron_point, point);

    double const difference =
        std::abs(first_tetrahedron_volume + second_tetrahedron_volume +
                 third_tetrahedron_volume + fourth_tetrahedron_volume -
                 total_tetrahedron_volume);

    return difference < epsilon;
}

double MathUtils::getBasisFunction(
    const std::vector<Point>& nodes, Point point,
    const std::vector<int>& node_indexes, BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_elements_type,
    int number_basis_function) {
    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    return getLinearBasisFunction(
                        point, nodes[node_indexes.at(0)],
                        nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                        nodes[node_indexes.at(3)], nodes[node_indexes.at(4)],
                        number_basis_function);
                    break;

                case BASIS_ELEMENT_TYPE::Quadratic:
                    break;

                default:
                    return 0.0;
                    break;
            }

        case BASIS_TYPE::Hermite:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Cubic:
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

double MathUtils::getLinearBasisFunction(Point point, Point node_0,
                                         Point node_1, Point node_2,
                                         Point node_3, Point vertex,
                                         int number_basis_function) {
    Point const center_of_base = {
        (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
        (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
        (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    Point const vertice_0_minus_1 = {node_1.x - node_0.x, node_1.y - node_0.y,
                                     node_1.z - node_0.z};

    Point const vertice_0_minus_4 = {node_2.x - node_0.x, node_2.y - node_0.y,
                                     node_2.z - node_0.z};

    Point const vertice_0_minus_point = {point.x, point.y, point.z};

    double const norm_of_vertice_0_minus_1 =
        vertice_0_minus_1.x * vertice_0_minus_1.x +
        vertice_0_minus_1.y * vertice_0_minus_1.y +
        vertice_0_minus_1.z * vertice_0_minus_1.z;

    double const norm_of_vertice_0_minus_4 =
        vertice_0_minus_4.x * vertice_0_minus_4.x +
        vertice_0_minus_4.y * vertice_0_minus_4.y +
        vertice_0_minus_4.z * vertice_0_minus_4.z;

    double const ksi_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_1.x +
         vertice_0_minus_point.y * vertice_0_minus_1.y +
         vertice_0_minus_point.z * vertice_0_minus_1.z) /
        norm_of_vertice_0_minus_1;

    double const nu_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_4.x +
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

    Point const vetice_point_minus_center = {point.x - center_of_base.x,
                                             point.y - center_of_base.y,
                                             point.z - center_of_base.z};

    double const tetta_parameter =
        (vetice_point_minus_center.x * vetice_vertex_minus_center.x +
         vetice_point_minus_center.y * vetice_vertex_minus_center.y +
         vetice_point_minus_center.z * vetice_vertex_minus_center.z) /
        norm_of_vetice_vertex_minus_center;

    switch (number_basis_function) {
        case 0:
            return (1 - ksi_parameter) * (1 - nu_parameter) *
                   (1 - tetta_parameter);
            break;
        case 1:
            return ksi_parameter * (1 - nu_parameter) * (1 - tetta_parameter);
            break;
        case 2:
            return (1 - ksi_parameter) * nu_parameter * (1 - tetta_parameter);
            break;
        case 3:
            return ksi_parameter * nu_parameter * (1 - tetta_parameter);
            break;
        case 4:
            return tetta_parameter;
            break;
        default:
            return 0.0;
    }
}

double MathUtils::getQuadraticBasisFunction(Point point, Point node_0,
                                            Point node_1, Point node_2,
                                            Point node_3, Point vertex,
                                            int number_basis_function) {
    Point const center_of_base = {
        (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
        (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
        (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    Point const vertice_0_minus_1 = {node_1.x - node_0.x, node_1.y - node_0.y,
                                     node_1.z - node_0.z};

    Point const vertice_0_minus_4 = {node_2.x - node_0.x, node_2.y - node_0.y,
                                     node_2.z - node_0.z};

    Point const vertice_0_minus_point = {point.x, point.y, point.z};

    double const norm_of_vertice_0_minus_1 =
        vertice_0_minus_1.x * vertice_0_minus_1.x +
        vertice_0_minus_1.y * vertice_0_minus_1.y +
        vertice_0_minus_1.z * vertice_0_minus_1.z;

    double const norm_of_vertice_0_minus_4 =
        vertice_0_minus_4.x * vertice_0_minus_4.x +
        vertice_0_minus_4.y * vertice_0_minus_4.y +
        vertice_0_minus_4.z * vertice_0_minus_4.z;

    double const ksi_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_1.x +
         vertice_0_minus_point.y * vertice_0_minus_1.y +
         vertice_0_minus_point.z * vertice_0_minus_1.z) /
        norm_of_vertice_0_minus_1;

    double const nu_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_4.x +
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

    Point const vetice_point_minus_center = {point.x - center_of_base.x,
                                             point.y - center_of_base.y,
                                             point.z - center_of_base.z};

    double const tetta_parameter =
        (vetice_point_minus_center.x * vetice_vertex_minus_center.x +
         vetice_point_minus_center.y * vetice_vertex_minus_center.y +
         vetice_point_minus_center.z * vetice_vertex_minus_center.z) /
        norm_of_vetice_vertex_minus_center;

    switch (number_basis_function) {
        case 0:
            return (1 - ksi_parameter) * (1 - nu_parameter) *
                   (1 - tetta_parameter);
            break;
        case 1:
            return ksi_parameter * (1 - nu_parameter) * (1 - tetta_parameter);
            break;
        case 2:
            return (1 - ksi_parameter) * nu_parameter * (1 - tetta_parameter);
            break;
        case 3:
            return ksi_parameter * nu_parameter * (1 - tetta_parameter);
            break;
        case 4:
            return tetta_parameter;
            break;
        default:
            return 0.0;
    }
}

double MathUtils::getDerivativeLinearBasisFunction(
    int number_basis_function, Point point, Point node_0, Point node_1,
    Point node_2, Point node_3, Point vertex, int number_derivative_parameter) {
    Point const center_of_base = {
        (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
        (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
        (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

    Point const vertice_0_minus_1 = {node_1.x - node_0.x, node_1.y - node_0.y,
                                     node_1.z - node_0.z};

    Point const vertice_0_minus_4 = {node_2.x - node_0.x, node_2.y - node_0.y,
                                     node_2.z - node_0.z};

    Point const vertice_0_minus_point = {point.x, point.y, point.z};

    double const norm_of_vertice_0_minus_1 =
        vertice_0_minus_1.x * vertice_0_minus_1.x +
        vertice_0_minus_1.y * vertice_0_minus_1.y +
        vertice_0_minus_1.z * vertice_0_minus_1.z;

    double const norm_of_vertice_0_minus_4 =
        vertice_0_minus_4.x * vertice_0_minus_4.x +
        vertice_0_minus_4.y * vertice_0_minus_4.y +
        vertice_0_minus_4.z * vertice_0_minus_4.z;

    double const ksi_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_1.x +
         vertice_0_minus_point.y * vertice_0_minus_1.y +
         vertice_0_minus_point.z * vertice_0_minus_1.z) /
        norm_of_vertice_0_minus_1;

    double const nu_parameter =
        (vertice_0_minus_point.x * vertice_0_minus_4.x +
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

    Point const vetice_point_minus_center = {point.x - center_of_base.x,
                                             point.y - center_of_base.y,
                                             point.z - center_of_base.z};

    double const tetta_parameter =
        (vetice_point_minus_center.x * vetice_vertex_minus_center.x +
         vetice_point_minus_center.y * vetice_vertex_minus_center.y +
         vetice_point_minus_center.z * vetice_vertex_minus_center.z) /
        norm_of_vetice_vertex_minus_center;

    std::vector<std::vector<double>> jacobian(3, std::vector<double>(3));

    std::vector<std::vector<double>> inversed_jacobian(3,
                                                       std::vector<double>(3));

    calculateJacobian(node_0, node_1, node_2, vertex, jacobian);

    inverseMatrix(jacobian, inversed_jacobian);

    double d_phi_d_ksi = 0.0;
    double d_phi_d_nu = 0.0;
    double d_phi_d_tetta = 0.0;

    switch (number_basis_function) {
        case 0:
            d_phi_d_ksi = -(1 - nu_parameter) * (1 - tetta_parameter);
            d_phi_d_nu = -(1 - ksi_parameter) * (1 - tetta_parameter);
            d_phi_d_tetta = -(1 - ksi_parameter) * (1 - nu_parameter);
            break;
        case 1:
            d_phi_d_ksi = (1 - nu_parameter) * (1 - tetta_parameter);
            d_phi_d_nu = -ksi_parameter * (1 - tetta_parameter);
            d_phi_d_tetta = -ksi_parameter * (1 - nu_parameter);
            break;
        case 2:
            d_phi_d_ksi = -nu_parameter * (1 - tetta_parameter);
            d_phi_d_nu = (1 - ksi_parameter) * (1 - tetta_parameter);
            d_phi_d_tetta = -(1 - ksi_parameter) * nu_parameter;
            break;
        case 3:
            d_phi_d_ksi = nu_parameter * (1 - tetta_parameter);
            d_phi_d_nu = ksi_parameter * (1 - tetta_parameter);
            d_phi_d_tetta = -ksi_parameter * nu_parameter;
            break;
        case 4:
            d_phi_d_ksi = 0;
            d_phi_d_nu = 0;
            d_phi_d_tetta = 1;
            break;
        default:
            break;
    }

    Point const derivatives = {d_phi_d_ksi * inversed_jacobian[0][0] +
                                   d_phi_d_nu * inversed_jacobian[1][0] +
                                   d_phi_d_tetta * inversed_jacobian[2][0],
                               d_phi_d_ksi * inversed_jacobian[0][1] +
                                   d_phi_d_nu * inversed_jacobian[1][1] +
                                   d_phi_d_tetta * inversed_jacobian[2][1],
                               d_phi_d_ksi * inversed_jacobian[0][2] +
                                   d_phi_d_nu * inversed_jacobian[1][2] +
                                   d_phi_d_tetta * inversed_jacobian[2][2]};

    switch (number_derivative_parameter) {
        case 0:
            return derivatives.x;
            break;
        case 1:
            return derivatives.y;
            break;
        case 2:
            return derivatives.z;
            break;
        default:
            return 0.0;
            break;
    }
}

void MathUtils::calculateJacobian(Point node_0, Point node_1, Point node_2,
                                  Point node_4,
                                  std::vector<std::vector<double>>& jacobian) {
    jacobian[0][0] = node_1.x - node_0.x;
    jacobian[0][1] = node_2.x - node_0.x;
    jacobian[0][2] = node_4.x - node_0.x;

    jacobian[1][0] = node_1.y - node_0.y;
    jacobian[1][1] = node_2.y - node_0.y;
    jacobian[1][2] = node_4.y - node_0.y;

    jacobian[2][0] = node_1.z - node_0.z;
    jacobian[2][1] = node_2.z - node_0.z;
    jacobian[2][2] = node_4.z - node_0.z;
}

double MathUtils::calculateDeterminant(
    std::vector<std::vector<double>>& matrix) {
    return matrix[0][0] *
               (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] *
               (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] *
               (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

void MathUtils::inverseMatrix(
    std::vector<std::vector<double>>& matrix,
    std::vector<std::vector<double>>& inversed_matrix) {
    double const determinant = calculateDeterminant(matrix);

    inversed_matrix[0][0] =
        (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) /
        determinant;
    inversed_matrix[0][1] =
        (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) /
        determinant;
    inversed_matrix[0][2] =
        (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) /
        determinant;

    inversed_matrix[1][0] =
        (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) /
        determinant;
    inversed_matrix[1][1] =
        (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) /
        determinant;
    inversed_matrix[1][2] =
        (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) /
        determinant;

    inversed_matrix[2][0] =
        (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) /
        determinant;
    inversed_matrix[2][1] =
        (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) /
        determinant;
    inversed_matrix[2][2] =
        (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) /
        determinant;
}
