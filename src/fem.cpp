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

#include "../include/element.h"
#include "../include/grid.h"
#include "../include/integrator.h"
#include "../include/mathUtils.h"

static constexpr std::string_view RESULT_OUTPUT_FILE_NAME =
    "../data/output/result.txt";

const double MIDPOINT_DIVISOR = 2.0;

int FEM::getFiniteElementIndex(Point point) {
    int position = -1;

#pragma unroll 4
    for (int i = 0; i < finite_elements.size(); ++i) {
        if (MathUtils::isPointInPyramid(nodes, point, finite_elements[i])) {
            position = i;
        }
    }

    return position;
}

double FEM::getResultAtPoint(Point point) {
    double result = 0.0;

    int const finite_element_index = getFiniteElementIndex(point);

    const auto slae_result = slae.getResultVector();

    const auto node_indexes =
        finite_elements[finite_element_index].getNodeIndexes();

#pragma unroll 4
    for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
        result += slae_result[node_indexes.at(i)] *
                  MathUtils::getLinearBasisFunction(
                      point, nodes[node_indexes.at(0)],
                      nodes[node_indexes.at(1)], nodes[node_indexes.at(2)],
                      nodes[node_indexes.at(3)], nodes[nodes.size() - 1], i);
    }

    return result;
}

void FEM::generateLinearNodes() {
    const std::vector<double> grid_x = mesh.getGridData('x');
    const std::vector<double> grid_y = mesh.getGridData('y');
    const std::vector<double> grid_z = mesh.getGridData('z');

    const Point pyramid_height = {(grid_x[1] + grid_x[0]) / MIDPOINT_DIVISOR,
                                  (grid_y[1] + grid_y[0]) / MIDPOINT_DIVISOR,
                                  (grid_z[1] + grid_z[0]) / MIDPOINT_DIVISOR};

    for (const double& z_element : grid_z) {
        for (const double& y_element : grid_y) {
#pragma unroll 4
            for (const double& x_element : grid_x) {
                nodes.push_back({x_element, y_element, z_element});
            }
        }
    }

    nodes.push_back(pyramid_height);
}

void FEM::generateQuadraticNodes() {
    const std::vector<double> grid_x = mesh.getGridData('x');
    const std::vector<double> grid_y = mesh.getGridData('y');
    const std::vector<double> grid_z = mesh.getGridData('z');

    const Point pyramid_height = {(grid_x[4] + grid_x[0]) / MIDPOINT_DIVISOR,
                                  (grid_y[4] + grid_y[0]) / MIDPOINT_DIVISOR,
                                  (grid_z[4] + grid_z[0]) / MIDPOINT_DIVISOR};

    const double center_of_line = 0.5;

    for (int i = 0; i < grid_z.size(); ++i) {
        for (int j = (i % 2); j < grid_y.size(); j += 2) {
#pragma unroll 4
            for (int k = (i % 2); k < grid_x.size(); k += 2) {
                if (!MathUtils::isQudraticNodeIllegal(std::make_tuple(
                        i, grid_x[k], grid_y[j], center_of_line))) {
                    nodes.push_back({grid_x[k], grid_y[j], grid_z[i]});
                }
            }
        }
    }

    nodes.push_back(pyramid_height);
}

void FEM::generateCubicNodes() {
    const std::vector<double> grid_x = mesh.getGridData('x');
    const std::vector<double> grid_y = mesh.getGridData('y');
    const std::vector<double> grid_z = mesh.getGridData('z');

    const Point pyramid_height = {(grid_x[5] + grid_x[0]) / MIDPOINT_DIVISOR,
                                  (grid_y[5] + grid_y[0]) / MIDPOINT_DIVISOR,
                                  (grid_z[5] + grid_z[0]) / MIDPOINT_DIVISOR};

    for (int i = 0; i < grid_z.size(); ++i) {
        for (int j = 0; j < grid_y.size(); j++) {
#pragma unroll 4
            for (int k = 0; k < grid_x.size(); k++) {
                if (k == 0 && (i == 1 || i == 4)) {
                    if (j == 1 || j == 4) {
                        nodes.push_back({grid_x[1], grid_y[j], grid_z[i]});
                        nodes.push_back({grid_x[4], grid_y[j], grid_z[i]});
                    } else if (j != 0 && j != grid_z.size() - 1) {
                        nodes.push_back({grid_x[2], grid_y[j], grid_z[i]});
                        nodes.push_back({grid_x[3], grid_y[j], grid_z[i]});
                    }
                } else if (!MathUtils::isCubicNodeIllegal(
                               std::make_tuple(i, j, k))) {
                    nodes.push_back({grid_x[k], grid_y[j], grid_z[i]});
                }
            }
        }
    }

    nodes.push_back(pyramid_height);
}

void FEM::generateFEMData(std::string_view const& input_file_name,
                          BASIS_TYPE basis_functions_type,
                          BASIS_ELEMENT_TYPE basis_functions_elements_type) {
    std::filesystem::path const file_name = input_file_name;

    mesh.generateGrid(file_name, basis_functions_type,
                      basis_functions_elements_type);

    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    generateLinearNodes();
                    break;

                case BASIS_ELEMENT_TYPE::Quadratic:
                    generateQuadraticNodes();
                    break;

                case BASIS_ELEMENT_TYPE::Cubic:
                    generateCubicNodes();
                    break;

                default:
                    break;
            }
            break;

        default:
            generateLinearNodes();
            break;
    }

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

                    bool const is_z_plane = MathUtils::isPlane(
                        nodes, average_z_coordinate,
                        {node_1, node_2, node_3, node_4}, 'z');

                    double const average_x_coordinate =
                        (nodes[node_1].x + nodes[node_2].x + nodes[node_3].x +
                         nodes[node_4].x) /
                        4.0;

                    bool const is_x_plane = MathUtils::isPlane(
                        nodes, average_x_coordinate,
                        {node_1, node_2, node_3, node_4}, 'x');

                    double const average_y_coordinate =
                        (nodes[node_1].y + nodes[node_2].y + nodes[node_3].y +
                         nodes[node_4].y) /
                        4.0;

                    bool const is_y_plane = MathUtils::isPlane(
                        nodes, average_y_coordinate,
                        {node_1, node_2, node_3, node_4}, 'y');

                    if (is_z_plane || is_x_plane || is_y_plane) {
                        finite_elements.emplace_back(Element(
                            {node_1, node_2, node_3, node_4, vertex_index}));
                    }
                }
            }
        }
    }

    number_of_vertices_of_pyramid = finite_elements[0].getNodeIndexes().size();

    nodes_count = nodes.size();

    finite_elements_count = finite_elements.size();

    stiffness_matrix_local.resize(
        number_of_vertices_of_pyramid,
        std::vector<double>(number_of_vertices_of_pyramid));

    mass_matrix_local.resize(
        number_of_vertices_of_pyramid,
        std::vector<double>(number_of_vertices_of_pyramid));

    right_part_local.resize(number_of_vertices_of_pyramid);

    saveGridForVisualize();
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
            first_boundary_condition_nodes[i].second = MathUtils::calculateU(
                nodes[first_boundary_condition_nodes[i].first]);
        }
    }
}

void FEM::generatePortrait() {
    std::size_t lower_triangle_size = 0;

    std::vector<std::set<int>> connections(nodes_count);

    for (int i = 0; i < finite_elements_count; i++) {
        for (int j = 0; j < number_of_vertices_of_pyramid; j++) {
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
    for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
#pragma unroll 4
        for (int j = 0; j < number_of_vertices_of_pyramid; j++) {
            stiffness_matrix_local[i][j] =
                lambda * Integrator::integrateForStiffnessMatrix(
                             nodes, finite_elements,
                             std::make_tuple(index_of_finite_element, i, j));

            mass_matrix_local[i][j] =
                gamma * Integrator::integrateForMassMatrix(
                            nodes, finite_elements,
                            std::make_tuple(index_of_finite_element, i, j));
        }
    }

#pragma unroll 4
    for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
        right_part_local[i] = Integrator::integrateForRightPartVector(
            nodes, finite_elements,
            std::make_tuple(index_of_finite_element, i));
    }
}

void FEM::assemblyGlobalComponents() {
    for (int i = 0; i < finite_elements_count; i++) {
        calculateLocalComponents(i);

        auto const node_indexes = finite_elements[i].getNodeIndexes();

        for (int j = 0; j < number_of_vertices_of_pyramid; j++) {
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

void FEM::checkBasisFunctionsToEqualsOne(
    BASIS_TYPE basis_functions_type,
    BASIS_ELEMENT_TYPE basis_functions_elements_type) {
    bool is_all_basis_functions_good = true;

    std::vector<int> wrong_basis_function_indexes;

    for (const auto& element : finite_elements) {
        const auto node_indexes = element.getNodeIndexes();

#pragma unroll 4
        for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
            const int current_pyramid_point = node_indexes.at(i);
            if (MathUtils::getBasisFunction(nodes, nodes[node_indexes.at(i)],
                                            node_indexes, basis_functions_type,
                                            basis_functions_elements_type,
                                            i) != 1.0) {
                is_all_basis_functions_good = false;

                wrong_basis_function_indexes.push_back(i);
            }
        }

        if (!wrong_basis_function_indexes.empty()) {
            std::cout << "Detected wrong basis function in element = {"
                      << node_indexes.at(0) << ", " << node_indexes.at(1)
                      << ", " << node_indexes.at(2) << ", "
                      << node_indexes.at(3) << ", " << node_indexes.at(4)
                      << "}\nat points = (";

#pragma unroll 4
            for (int i = 0; i < wrong_basis_function_indexes.size(); i++) {
                std::cout << nodes[wrong_basis_function_indexes[i]].x << ", "
                          << nodes[wrong_basis_function_indexes[i]].y << ", "
                          << nodes[wrong_basis_function_indexes[i]].z << ")";

                if (i != wrong_basis_function_indexes.size() - 1) {
                    std::cout << ", ";
                }
            }

            std::cout << '\n';
        }
    }

    if (is_all_basis_functions_good) {
        std::cout << "All basis functions are correct for all elements!\n";
    }
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
            double const true_result = MathUtils::calculateU(point);
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

void FEM::saveGridForVisualize() {
    int const nodes_out_step = 10;
    int const elements_out_step = 5;

    {
        std::ofstream out_nodes("../data/output/nodes.txt");

#pragma unroll 4
        for (auto& node : nodes) {
            out_nodes << std::setw(nodes_out_step) << node.x
                      << std::setw(nodes_out_step) << node.y
                      << std::setw(nodes_out_step) << node.z << '\n';
        }
    }

    {
        std::ofstream out_elements("../data/output/elements.txt");

#pragma unroll 4
        for (auto& element : finite_elements) {
            const auto node_indexes = element.getNodeIndexes();
            out_elements << std::setw(elements_out_step) << node_indexes.at(0)
                         << std::setw(elements_out_step) << node_indexes.at(1)
                         << std::setw(elements_out_step) << node_indexes.at(2)
                         << std::setw(elements_out_step) << node_indexes.at(3)
                         << std::setw(elements_out_step) << node_indexes.at(4)
                         << '\n';
        }
    }
}