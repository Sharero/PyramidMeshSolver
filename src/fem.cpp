#include "../include/fem.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string_view>
#include <tuple>
#include <vector>

#include "../include/grid.h"

static constexpr std::string_view RESULT_OUTPUT_FILE_NAME =
    "../data/output/result.txt";

const double MIDPOINT_DIVISOR = 2.0;

const double NUMBERS_EQUAL_EPSILON = 1e-15;

const double TETRAHEDRON_VOLUME_DIVISOR = 6.0;

const double INTEGRAL_NORMALIZATION_FACTOR = 0.5;

double FEM::calculateF(Point point) {
    const double const_right_part = 5.0;
    const double quadratic_div_grad = -6.0;

    // return const_right_part + 0 * (point.x + point.y + point.z);
    return point.x + point.y + point.z;
    // return quadratic_div_grad + pow(point.x, 2) + pow(point.y, 2) +
    //        pow(point.z, 2);
}

double FEM::calculateU(Point point) {
    const double const_right_part = 5.0;

    // return const_right_part + 0 * (point.x + point.y + point.z);
    return point.x + point.y + point.z;
    // return pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2);
    // return pow(point.x, 3) + pow(point.y, 3) + pow(point.z, 3);
    // return sin(point.x * point.y * point.z);
}

void FEM::saveTestResults(const std::vector<Point>& test_points) {
    int const first_column_step = 12;
    int const second_column_step = 15;
    int const third_column_step = 18;

    int const numbers_in_float = 6;

    {
        std::filesystem::path const file_name = RESULT_OUTPUT_FILE_NAME;
        std::ofstream result_out(file_name);

        result_out << std::setw(first_column_step) << "x"
                   << std::setw(second_column_step) << "y"
                   << std::setw(second_column_step) << "z"
                   << std::setw(second_column_step) << "u"
                   << std::setw(second_column_step) << "u*"
                   << std::setw(third_column_step) << "|u - u*|" << '\n';

        result_out << std::fixed << std::setprecision(numbers_in_float);

#pragma unroll 4
        for (const auto& point : test_points) {
            double const true_result = calculateU(point);
            double const result = getResultAtPoint(point);

            result_out << std::setw(second_column_step) << point.x
                       << std::setw(second_column_step) << point.y
                       << std::setw(second_column_step) << point.z
                       << std::setw(second_column_step) << result
                       << std::setw(second_column_step) << true_result
                       << std::setw(second_column_step)
                       << std::abs(true_result - result) << '\n';
        }
    }
}

bool FEM::isNumbersEqual(double first_number, double second_number,
                         double epsilon) {
    return std::fabs(first_number - second_number) < epsilon;
}

bool FEM::isPlane(double average, const Element& element, char axis) {
    const auto node_indexes = element.getNodeIndexes();

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

void FEM::generateFEMData(std::string_view const& input_file_name) {
    std::filesystem::path const file_name = input_file_name;

    mesh.generateGrid(file_name);

    const std::vector<double> grid_x = mesh.getGridData('x');
    const std::vector<double> grid_y = mesh.getGridData('y');
    const std::vector<double> grid_z = mesh.getGridData('z');

    const Point pyramid_height = {
        (grid_x[grid_x.size() - 1] + grid_x[0]) / MIDPOINT_DIVISOR,
        (grid_y[grid_y.size() - 1] + grid_y[0]) / MIDPOINT_DIVISOR,
        (grid_z[grid_z.size() - 1] + grid_z[0]) / MIDPOINT_DIVISOR};

    for (const double& z_element : grid_z) {
        for (const double& y_element : grid_y) {
#pragma unroll 4
            for (const double& x_element : grid_x) {
                nodes.push_back({x_element, y_element, z_element});
            }
        }
    }

    nodes.push_back(pyramid_height);

    int const vertex_index = static_cast<int>(nodes.size()) - 1;
    int const current_element_index = 0;

    for (int node_1 = 0; node_1 < vertex_index; ++node_1) {
        for (int node_2 = node_1 + 1; node_2 < vertex_index; ++node_2) {
            for (int node_3 = node_2 + 1; node_3 < vertex_index; ++node_3) {
#pragma unroll 4
                for (int node_4 = node_3 + 1; node_4 < vertex_index; ++node_4) {
                    double const average_z_coordinate =
                        (nodes[node_1].z + nodes[node_2].z + nodes[node_3].z +
                         nodes[node_4].z) /
                        4.0;

                    bool const is_z_plane = isPlane(
                        average_z_coordinate,
                        Element(node_1, node_2, node_3, node_4, node_1), 'z');

                    double const average_x_coordinate =
                        (nodes[node_1].x + nodes[node_2].x + nodes[node_3].x +
                         nodes[node_4].x) /
                        4.0;

                    bool const is_x_plane = isPlane(
                        average_x_coordinate,
                        Element(node_1, node_2, node_3, node_4, node_1), 'x');

                    double const average_y_coordinate =
                        (nodes[node_1].y + nodes[node_2].y + nodes[node_3].y +
                         nodes[node_4].y) /
                        4.0;

                    bool const is_y_plane = isPlane(
                        average_y_coordinate,
                        Element(node_1, node_2, node_3, node_4, node_1), 'y');

                    if (is_z_plane || is_x_plane || is_y_plane) {
                        finite_elements.emplace_back(node_1, node_2, node_3,
                                                     node_4, vertex_index);
                    }
                }
            }
        }
    }

    nodes_count = nodes.size();

    finite_elements_count = finite_elements.size();

    stiffness_matrix_local.resize(
        NUMBER_OF_VERTICES_OF_PYRAMID,
        std::vector<double>(NUMBER_OF_VERTICES_OF_PYRAMID));

    mass_matrix_local.resize(
        NUMBER_OF_VERTICES_OF_PYRAMID,
        std::vector<double>(NUMBER_OF_VERTICES_OF_PYRAMID));

    right_part_local.resize(NUMBER_OF_VERTICES_OF_PYRAMID);
}

void FEM::inputBoundaryConditions(std::string_view const& input_file_name) {
    int count_first_condition_nodes = 0;
    std::filesystem::path const file_name = input_file_name;

    {
        std::ifstream input_boundaries(file_name);

        input_boundaries >> count_first_condition_nodes;

        first_boundary_condition_nodes.resize(count_first_condition_nodes);

#pragma unroll 4
        for (int i = 0; i < count_first_condition_nodes; i++) {
            input_boundaries >> first_boundary_condition_nodes[i].first;
            first_boundary_condition_nodes[i].second =
                calculateU(nodes[first_boundary_condition_nodes[i].first]);
        }
    }
}

void FEM::generatePortrait() {
    std::size_t lower_triangle_size = 0;

    std::vector<std::set<int>> connections(nodes_count);

    for (int i = 0; i < finite_elements_count; i++) {
        for (int j = 0; j < NUMBER_OF_VERTICES_OF_PYRAMID; j++) {
#pragma unroll 4
            for (int k = 0; k < j; k++) {
                const auto node_indexes = finite_elements[i].getNodeIndexes();
                connections[node_indexes.at(k)].insert(node_indexes.at(j));
            }
        }
    }

#pragma unroll 4
    for (int i = 0; i < nodes_count; i++) {
        lower_triangle_size += connections[i].size();
    }

    slae.setSlaeSize(nodes_count);

    slae.allocateMemory(lower_triangle_size);

    slae.setIndexArrays(connections);
}

void FEM::calculateLocalComponents(int index_of_finite_element) {
    for (int i = 0; i < NUMBER_OF_VERTICES_OF_PYRAMID; i++) {
#pragma unroll 4
        for (int j = 0; j < NUMBER_OF_VERTICES_OF_PYRAMID; j++) {
            stiffness_matrix_local[i][j] =
                lambda * integrateDerivativeBasisFunctions(
                             std::make_tuple(index_of_finite_element, i, j));

            mass_matrix_local[i][j] =
                gamma * integrateBasisFunctions(
                            std::make_tuple(index_of_finite_element, i, j));
        }
    }

#pragma unroll 4
    for (int i = 0; i < NUMBER_OF_VERTICES_OF_PYRAMID; i++) {
        right_part_local[i] = integrateBasisFunctionForF(
            std::make_tuple(index_of_finite_element, i));
    }
}

void FEM::assemblyGlobalComponents() {
    for (int i = 0; i < finite_elements_count; i++) {
        calculateLocalComponents(i);

        auto const node_indexes = finite_elements[i].getNodeIndexes();

        for (int j = 0; j < NUMBER_OF_VERTICES_OF_PYRAMID; j++) {
#pragma unroll 4
            for (int k = 0; k < j; k++) {
                double const lower_triangle_element =
                    mass_matrix_local[j][k] + stiffness_matrix_local[j][k];

                if (node_indexes.at(j) > node_indexes.at(k)) {
                    slae.addLowerTriangleElement(
                        std::make_tuple(node_indexes.at(j), node_indexes.at(k),
                                        lower_triangle_element));
                } else {
                    slae.addLowerTriangleElement(
                        std::make_tuple(node_indexes.at(k), node_indexes.at(j),
                                        lower_triangle_element));
                }
            }

            slae.addDiagonalElement(
                node_indexes.at(j),
                mass_matrix_local[j][j] + stiffness_matrix_local[j][j]);

            slae.addRightPartElement(node_indexes.at(j), right_part_local[j]);
        }
    }
}

void FEM::solveFEM() {
    saveGridForVisualize();

    generatePortrait();

    assemblyGlobalComponents();

    slae.applyFirstBoundaryConditions(first_boundary_condition_nodes);

    slae.solveSLAE();
}

void FEM::saveGridForVisualize() {
    {
        std::ofstream out_nodes("..data/output/nodes.txt");

#pragma unroll 4
        for (auto& node : nodes) {
            out_nodes << node.x << " " << node.y << " " << node.z << '\n';
        }
    }

    {
        std::ofstream out_elements("..data/output/elements.txt");

#pragma unroll 4
        for (auto& element : finite_elements) {
            const auto node_indexes = element.getNodeIndexes();
            out_elements << node_indexes.at(0) << " " << node_indexes.at(1)
                         << " " << node_indexes.at(2) << " "
                         << node_indexes.at(3) << " " << node_indexes.at(4)
                         << '\n';
        }
    }
}

int FEM::getFiniteElementIndex(Point point) {
    int position = -1;

#pragma unroll 4
    for (int i = 0; i < finite_elements.size(); ++i) {
        if (isPointInPyramid(point, finite_elements[i])) {
            position = i;
        }
    }

    return position;
}

Point FEM::calculateVectorCrossProduct(Point first_vector,
                                       Point second_vector) {
    Point cross_product_result = {
        first_vector.y * second_vector.z - first_vector.z * second_vector.y,
        first_vector.z * second_vector.x - first_vector.x * second_vector.z,
        first_vector.x * second_vector.y - first_vector.y * second_vector.x};

    return cross_product_result;
}

bool FEM::isPointInPyramid(const Point& point, const Element& element) {
    const auto node_indexes = element.getNodeIndexes();

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

double FEM::calculateTetrahedronVolume(Point first_tetrahedron_point,
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

bool FEM::isPointInsideTetrahedron(Point point, Point first_tetrahedron_point,
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

double FEM::getResultAtPoint(Point point) {
    double result = 0.0;

    int const finite_element_index = getFiniteElementIndex(point);

    const auto slae_result = slae.getResultVector();

    const auto node_indexes =
        finite_elements[finite_element_index].getNodeIndexes();

#pragma unroll 4
    for (int i = 0; i < NUMBER_OF_VERTICES_OF_PYRAMID; i++) {
        result += slae_result[node_indexes.at(i)] *
                  getLinearBasisFunctionTest(point, nodes[node_indexes.at(0)],
                                             nodes[node_indexes.at(1)],
                                             nodes[node_indexes.at(2)],
                                             nodes[node_indexes.at(3)], i);
    }

    return result;
}

double FEM::getLinearBasisFunctionTest(Point point, Point node_0, Point node_1,
                                       Point node_2, Point node_3,
                                       int number_basis_function) {
    Point const vertex = {0.5, 0.5, 0.5};

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

double FEM::getDerivativeLinearBasisFunctionTest(
    int number_basis_function, Point point, Point node_0, Point node_1,
    Point node_2, Point node_3, int number_derivative_parameter) {
    Point const vertex = {0.5, 0.5, 0.5};

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

void FEM::calculateJacobian(Point node_0, Point node_1, Point node_2,
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

double FEM::calculateDeterminant(std::vector<std::vector<double>>& matrix) {
    return matrix[0][0] *
               (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] *
               (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] *
               (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

void FEM::inverseMatrix(std::vector<std::vector<double>>& matrix,
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

double FEM::integrateBasisFunctions(
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
                    getLinearBasisFunctionTest(point, nodes[node_indexes.at(0)],
                                               nodes[node_indexes.at(1)],
                                               nodes[node_indexes.at(2)],
                                               nodes[node_indexes.at(3)],
                                               index_of_first_basis_function) *
                    getLinearBasisFunctionTest(point, nodes[node_indexes.at(0)],
                                               nodes[node_indexes.at(1)],
                                               nodes[node_indexes.at(2)],
                                               nodes[node_indexes.at(3)],
                                               index_of_second_basis_function);
            }
        }
    }

    integral *= x_integral_step * y_integral_step * z_integral_step;

    return integral;
}

double FEM::integrateBasisFunctionForF(
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
                    getLinearBasisFunctionTest(
                        point, nodes[node_indexes.at(0)],
                        nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                        nodes[node_indexes.at(3)], index_of_basis_function) *
                    calculateF(point);
            }
        }
    }

    integral *= x_integral_step * y_integral_step * z_integral_step;

    return integral;
}

double FEM::integrateDerivativeBasisFunctions(
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
                    getDerivativeLinearBasisFunctionTest(
                        index_of_second_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        0);

                fourth_integral_component =
                    getDerivativeLinearBasisFunctionTest(
                        index_of_second_basis_function, point,
                        nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                        nodes[node_indexes.at(2)], nodes[node_indexes.at(3)],
                        1);

                sixth_integral_component = getDerivativeLinearBasisFunctionTest(
                    index_of_second_basis_function, point,
                    nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                    nodes[node_indexes.at(2)], nodes[node_indexes.at(3)], 2);

                first_integral_component = getDerivativeLinearBasisFunctionTest(
                    index_of_first_basis_function, point,
                    nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                    nodes[node_indexes.at(2)], nodes[node_indexes.at(3)], 0);

                third_integral_component = getDerivativeLinearBasisFunctionTest(
                    index_of_first_basis_function, point,
                    nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                    nodes[node_indexes.at(2)], nodes[node_indexes.at(3)], 1);

                fifth_integral_component = getDerivativeLinearBasisFunctionTest(
                    index_of_first_basis_function, point,
                    nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                    nodes[node_indexes.at(2)], nodes[node_indexes.at(3)], 2);

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
