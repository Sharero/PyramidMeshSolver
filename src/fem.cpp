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

                    bool const is_z_plane = MathUtils::isPlane(
                        nodes, average_z_coordinate,
                        Element(node_1, node_2, node_3, node_4, node_1), 'z');

                    double const average_x_coordinate =
                        (nodes[node_1].x + nodes[node_2].x + nodes[node_3].x +
                         nodes[node_4].x) /
                        4.0;

                    bool const is_x_plane = MathUtils::isPlane(
                        nodes, average_x_coordinate,
                        Element(node_1, node_2, node_3, node_4, node_1), 'x');

                    double const average_y_coordinate =
                        (nodes[node_1].y + nodes[node_2].y + nodes[node_3].y +
                         nodes[node_4].y) /
                        4.0;

                    bool const is_y_plane = MathUtils::isPlane(
                        nodes, average_y_coordinate,
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
            first_boundary_condition_nodes[i].second = MathUtils::calculateU(
                nodes[first_boundary_condition_nodes[i].first]);
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
                lambda * Integrator::integrateDerivativeBasisFunctions(
                             nodes, finite_elements,
                             std::make_tuple(index_of_finite_element, i, j));

            mass_matrix_local[i][j] =
                gamma * Integrator::integrateBasisFunctions(
                            nodes, finite_elements,
                            std::make_tuple(index_of_finite_element, i, j));
        }
    }

#pragma unroll 4
    for (int i = 0; i < NUMBER_OF_VERTICES_OF_PYRAMID; i++) {
        right_part_local[i] = Integrator::integrateBasisFunctionForF(
            nodes, finite_elements,
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
    for (int i = 0; i < NUMBER_OF_VERTICES_OF_PYRAMID; i++) {
        result +=
            slae_result[node_indexes.at(i)] *
            MathUtils::getLinearBasisFunctionTest(
                point, nodes[node_indexes.at(0)], nodes[node_indexes.at(1)],
                nodes[node_indexes.at(2)], nodes[node_indexes.at(3)], i);
    }

    return result;
}
