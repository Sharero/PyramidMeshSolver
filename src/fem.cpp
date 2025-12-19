#include "../include/fem.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
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
const double THIRD_DIVISOR = 3.0;

const int COUNT_BASE_NODES_AT_LINEAR_COMBINATION = 4;
const int COUNT_BASE_NODES_AT_QUADRATIC_COMBINATION = 9;
const int COUNT_BASE_NODES_AT_CUBIC_COMBINATION = 16;

int FEM::getFiniteElementIndex(Point point) {
    int position = 0;

    double d[6] = {point.z,     point.y,     point.x,
                   1 - point.y, 1 - point.x, 1 - point.z};

    double dmin = d[0];

    for (int i = 1; i < 6; ++i) {
        if (d[i] < dmin) {
            dmin = d[i];
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
                  MathUtils::getBasisFunction(
                      nodes, point, node_indexes, basis_functions_type,
                      basis_functions_element_type, i, true);
    }

    return result;
}

int FEM::findFaceNodeIndex(int node_1, int node_2, double divisor) const {
    double alpha = 0.0;
    double beta = 0.0;

    if (divisor == MIDPOINT_DIVISOR) {
        alpha = 1.0;
        beta = 1.0;
    } else if (divisor == THIRD_DIVISOR) {
        alpha = 1.0;
        beta = 2.0;
    } else if (divisor == MIDPOINT_DIVISOR * THIRD_DIVISOR) {
        alpha = 1.0;
        beta = 5.0;
    }

    const Point desired_point = {
        (alpha * nodes[node_1].x + beta * nodes[node_2].x) / divisor,
        (alpha * nodes[node_1].y + beta * nodes[node_2].y) / divisor,
        (alpha * nodes[node_1].z + beta * nodes[node_2].z) / divisor};

#pragma unroll 4
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i] == desired_point) {
            return static_cast<int>(i);
        }
    }

    return -1;
}

void FEM::generateBaseNodesCombinations(
    int combination_size, std::vector<Element>& base_nodes_combinations) {
    int const vertex_index = static_cast<int>(nodes.size()) - 1;

    std::set<int> forbidden_indexes;
    std::vector<std::set<int>> required_indexes_sets;

    switch (combination_size) {
        case COUNT_BASE_NODES_AT_LINEAR_COMBINATION: {
            const int right_top_index = 5;

            forbidden_indexes = {};

            required_indexes_sets = {
                {0, 1}, {0, 2}, {2, 3}, {1, 3}, {4, right_top_index}};

            break;
        }
        case COUNT_BASE_NODES_AT_QUADRATIC_COMBINATION: {
            const std::vector<int> first_forbidden_layer_indexes = {9, 10, 11,
                                                                    12};
            const std::vector<int> second_forbidden_layer_indexes = {21, 22, 23,
                                                                     24};
            const std::set<int> first_plane_start_indexes = {0, 1, 2};
            const std::set<int> second_plane_start_indexes = {0, 3, 6};
            const std::set<int> third_plane_start_indexes = {6, 7, 8};
            const std::set<int> fourth_plane_start_indexes = {2, 5, 8};
            const std::set<int> fifth_plane_start_indexes = {25, 26, 27};

#pragma unroll 4
            for (const auto index : first_forbidden_layer_indexes) {
                forbidden_indexes.insert(index);
            }

#pragma unroll 4
            for (const auto index : second_forbidden_layer_indexes) {
                forbidden_indexes.insert(index);
            }

            required_indexes_sets = {
                first_plane_start_indexes, second_plane_start_indexes,
                third_plane_start_indexes, fourth_plane_start_indexes,
                fifth_plane_start_indexes};

            break;
        }
        default: {
            forbidden_indexes = {};
            required_indexes_sets = {};
            break;
        }
    }

#pragma unroll 4
    for (const auto& required_indexes : required_indexes_sets) {
        std::vector<int> remaining_indexes;

        MathUtils::filterRemainingIndexesForBaseCombinations(
            vertex_index, forbidden_indexes, required_indexes,
            remaining_indexes);

        MathUtils::calculateCombinations(combination_size, required_indexes,
                                         remaining_indexes,
                                         base_nodes_combinations);
    }
}

void FEM::generateElements(int p) {
    int const vertex_index = (int)nodes.size() - 1;
    int const baseN = p + 1;
    int const baseCount = baseN * baseN;

    std::array<int, 4> corner = {0, baseN - 1, baseN * (baseN - 1),
                                 baseCount - 1};

    std::vector<Element> base_combs;

    generateBaseNodesCombinations(baseCount, base_combs);

    for (auto& comb : base_combs) {
        auto node_idxs = comb.getNodeIndexes();

        double avgZ = 0;
        double avgY = 0;
        double avgX = 0;

        for (int idx : node_idxs) {
            avgX += nodes[idx].x;
            avgY += nodes[idx].y;
            avgZ += nodes[idx].z;
        }

        avgX /= baseCount;
        avgY /= baseCount;
        avgZ /= baseCount;

        if (!(MathUtils::isPlane(nodes, avgZ, comb, 'z') ||
              MathUtils::isPlane(nodes, avgY, comb, 'y') ||
              MathUtils::isPlane(nodes, avgX, comb, 'x')))
            continue;

        std::vector<int> elem_idxs;

        for (auto index : node_idxs) {
            elem_idxs.push_back(index);
        }

        for (int n = 1; n < p; ++n) {
            for (int c : corner) {
                elem_idxs.push_back(findFaceNodeIndex(
                    node_idxs[c], vertex_index, n * MIDPOINT_DIVISOR));
            }
        }

        elem_idxs.push_back(vertex_index);

        Element E;
        E.setNodeIndexes(elem_idxs);
        finite_elements.push_back(std::move(E));
    }
}

void FEM::generateNodes(int basis_type) {
    const auto grid_x = mesh.getGridData('x');
    const auto grid_y = mesh.getGridData('y');
    const auto grid_z = mesh.getGridData('z');

    int stride = 0;
    int lastGridIndex = 0;

    switch (basis_type) {
        case 1:
            stride = 1;
            lastGridIndex = 1;
            break;
        case 2:
            stride = 2;
            lastGridIndex = 4;
            break;
        case 3:
            stride = 1;
            lastGridIndex = 5;
            break;
        default:
            break;
    }

    Point pyramid_height = {
        (grid_x[lastGridIndex] + grid_x[0]) / MIDPOINT_DIVISOR,
        (grid_y[lastGridIndex] + grid_y[0]) / MIDPOINT_DIVISOR,
        (grid_z[lastGridIndex] + grid_z[0]) / MIDPOINT_DIVISOR};

    auto illegal = [&](int i, int j, int k) {
        if (basis_type == 3) {
            return MathUtils::isCubicNodeIllegal(std::make_tuple(i, j, k));
        } else if (basis_type == 2) {
            return grid_x[k] == 0.5 && grid_y[j] == 0.5 && grid_z[i] == 0.5;
        }
        return false;
    };

    for (int i = 0; i < grid_z.size(); ++i) {
        int jStart = (basis_type == 2 ? (i % stride) : 0);

        for (int j = jStart; j < grid_y.size(); j += stride) {
            int kStart = (basis_type == 2 ? (i % stride) : 0);

            for (int k = kStart; k < grid_x.size(); k += stride) {
                if (!illegal(i, j, k)) {
                    nodes.push_back({grid_x[k], grid_y[j], grid_z[i]});
                }
            }
        }
    }

    nodes.push_back(pyramid_height);
}

void FEM::generateFEMData(std::string_view const& input_file_name) {
    std::filesystem::path const file_name = input_file_name;

    mesh.generateGrid(file_name, basis_functions_type,
                      basis_functions_element_type);

    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_element_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    generateNodes(1);
                    generateElements(1);
                    break;
                case BASIS_ELEMENT_TYPE::Quadratic:
                    generateNodes(2);
                    generateElements(2);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    generateNodes(3);
                    finite_elements = {
                        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                         13, 14, 15, 16, 17, 18, 19, 25, 26, 29, 30, 72},
                        {0,  1,  2,  3,  20, 21, 22, 23, 36, 37, 38, 39, 56,
                         57, 58, 59, 16, 17, 52, 53, 25, 26, 41, 42, 72},
                        {0,  4,  8,  12, 20, 24, 28, 32, 36, 40, 44, 48, 56,
                         60, 64, 68, 16, 18, 52, 54, 25, 29, 41, 45, 72},
                        {12, 13, 14, 15, 32, 33, 34, 35, 48, 49, 50, 51, 68,
                         69, 70, 71, 18, 19, 54, 55, 29, 30, 45, 46, 72},
                        {3,  7,  11, 15, 23, 27, 31, 35, 39, 43, 47, 51, 59,
                         63, 67, 71, 17, 19, 53, 55, 26, 30, 42, 46, 72},
                        {56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                         69, 70, 71, 52, 53, 54, 55, 41, 42, 45, 46, 72}};
                    break;
                default:
                    break;
            }
            break;
        case BASIS_TYPE::Hierarhical:
            switch (basis_functions_element_type) {
                case BASIS_ELEMENT_TYPE::Quadratic:
                    generateNodes(2);
                    generateElements(2);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    generateNodes(1);
                    generateElements(1);
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }

    writeGridInformationToFile();

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
}

void FEM::inputBoundaryConditions(std::string_view const& input_file_name) {
    int count_first_condition_nodes = 0;
    std::filesystem::path const file_name = input_file_name;

    {
        std::ifstream input_boundaries(file_name);

        input_boundaries >> count_first_condition_nodes;

        first_boundary_condition_nodes.resize(count_first_condition_nodes);

        std::vector<int> global_vertexes = {0, 2, 6, 8, 25, 27, 31, 33};
        std::vector<std::array<int, 3>> global_edges = {
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {13, 14, 15},
            {18, 19, 20}, {25, 26, 27}, {28, 29, 30}, {31, 32, 33},
            {0, 3, 6},    {2, 5, 8},    {13, 16, 18}, {15, 17, 20},
            {25, 28, 31}, {27, 30, 33}};

#pragma unroll 4
        for (int i = 0; i < count_first_condition_nodes; i++) {
            input_boundaries >> first_boundary_condition_nodes[i].first;
            if (find(global_vertexes.begin(), global_vertexes.end(),
                     first_boundary_condition_nodes[i].first) !=
                global_vertexes.end()) {
                first_boundary_condition_nodes[i].second =
                    MathUtils::calculateU(
                        nodes[first_boundary_condition_nodes[i].first]);
            } else {
                for (auto edge : global_edges) {
                    if (edge.at(1) == first_boundary_condition_nodes[i].first) {
                        double value_0 =
                            MathUtils::calculateU(nodes[edge.at(0)]);
                        double value_1 =
                            MathUtils::calculateU(nodes[edge.at(2)]);
                        double value_2 =
                            MathUtils::calculateU(nodes[edge.at(1)]);
                        first_boundary_condition_nodes[i].second =
                            value_2 - (value_0 + value_1) / 2.0;
                    }
                }
            }
        }
    }
}

void FEM::generatePortrait() {
    std::size_t lower_triangle_size = 0;

    std::vector<std::set<int>> connections(nodes_count);

    for (int i = 0; i < finite_elements_count; i++) {
        const auto node_indexes = finite_elements[i].getNodeIndexes();
        for (int j = 0; j < number_of_vertices_of_pyramid; j++) {
            const int index_1 = node_indexes.at(j);
            for (int k = 0; k < number_of_vertices_of_pyramid; k++) {
                const int index_2 = node_indexes.at(k);
                if (index_1 < index_2) {
                    connections[index_1].insert(index_2);
                }
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
        for (int j = i; j < number_of_vertices_of_pyramid; j++) {
            stiffness_matrix_local[i][j] =
                lambda * Integrator::integrateForStiffnessMatrix(
                             nodes, finite_elements,
                             std::make_tuple(index_of_finite_element, i, j),
                             basis_functions_type,
                             basis_functions_element_type);

            mass_matrix_local[i][j] =
                gamma * Integrator::integrateForMassMatrix(
                            nodes, finite_elements,
                            std::make_tuple(index_of_finite_element, i, j),
                            basis_functions_type, basis_functions_element_type);

            stiffness_matrix_local[j][i] = stiffness_matrix_local[i][j];
            mass_matrix_local[j][i] = mass_matrix_local[i][j];
        }
    }

    const auto node_indexes =
        finite_elements[index_of_finite_element].getNodeIndexes();

    for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
        right_part_local[i] = 0.0;
    }

    for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
        right_part_local[i] = Integrator::integrateForRightPartVector(
            nodes, finite_elements, std::make_tuple(index_of_finite_element, i),
            basis_functions_type, basis_functions_element_type);
    }

    writeLocalComponentsToFiles(index_of_finite_element, stiffness_matrix_local,
                                mass_matrix_local, right_part_local);
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
    generatePortrait();

    assemblyGlobalComponents();

    slae.applyFirstBoundaryConditions(first_boundary_condition_nodes);

    slae.solveSLAE();

    const auto slae_result = slae.getResultVector();

    for (const auto result : slae_result) {
        std::cout << result << '\n';
    }
}

void FEM::checkBasisFunctionsDeltaProperty() {
    bool all_ok = true;

    struct Error {
        std::size_t basis_function_index;
        std::size_t node_index;

        double calculated_value;
        double expected_value;
    };

    for (std::size_t i = 0; i < finite_elements_count; ++i) {
        const auto node_indexes = finite_elements[i].getNodeIndexes();

        std::vector<Error> errors;

        for (std::size_t j = 0; j < number_of_vertices_of_pyramid; ++j) {
            const Point& test_node = nodes[node_indexes[j]];

            double val_self = MathUtils::getBasisFunction(
                nodes, test_node, node_indexes, basis_functions_type,
                basis_functions_element_type, static_cast<int>(j), true);

            if (std::abs(val_self - 1.0) > 1e-12) {
                errors.push_back({j, j, val_self, 1.0});
            }

            for (std::size_t k = 0; k < number_of_vertices_of_pyramid; ++k) {
                if (k == j) {
                    continue;
                }

                const Point& other_node = nodes[node_indexes[k]];

                double val_other = MathUtils::getBasisFunction(
                    nodes, other_node, node_indexes, basis_functions_type,
                    basis_functions_element_type, static_cast<int>(j), true);

                if (std::abs(val_other) > 1e-12) {
                    errors.push_back({j, k, val_other, 0.0});
                }
            }
        }

        if (!errors.empty()) {
            all_ok = false;

            std::cout << "Element #" << i << " (nodes: ";

            for (auto idx : node_indexes) {
                std::cout << idx << ' ';
            }

            std::cout << ")\n";

            for (auto& err : errors) {
                std::cout << "  Basis #" << err.basis_function_index
                          << " at node idx " << node_indexes[err.node_index]
                          << ": got " << err.calculated_value << ", expected "
                          << err.expected_value << '\n';
            }
            std::cout << '\n';
        }
    }

    if (all_ok) {
        std::cout << "All basis functions passed the delta-property on all "
                     "elements.\n";
    }
}

void FEM::calculateResultAtPointsToFile(const std::vector<Point>& test_points) {
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

void FEM::writeGridInformationToFile() {
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

#pragma unroll 4
            for (const auto index : node_indexes) {
                out_elements << std::setw(elements_out_step) << index << " ";
            }

            out_elements << '\n';
        }
    }
}

void FEM::writeLocalComponentsToFiles(
    int index_of_finite_element,
    const std::vector<std::vector<double>>& stiffness_matrix,
    const std::vector<std::vector<double>>& mass_matrix,
    const std::vector<double>& right_part_vector) {
    std::ofstream output_file("../data/output/test" +
                              std::to_string(index_of_finite_element) + ".txt");

    for (int i = 0; i < number_of_vertices_of_pyramid; ++i) {
        for (int j = 0; j < number_of_vertices_of_pyramid; ++j) {
            output_file << std::setw(20) << stiffness_matrix[i][j];
        }
        output_file << "\n";
    }
    output_file << "\n";

    for (int i = 0; i < number_of_vertices_of_pyramid; ++i) {
        for (int j = 0; j < number_of_vertices_of_pyramid; ++j) {
            output_file << std::setw(20) << mass_matrix[i][j];
        }
        output_file << "\n";
    }
    output_file << "\n";

    output_file << std::setprecision(12);

    for (int i = 0; i < number_of_vertices_of_pyramid; ++i) {
        output_file << std::setw(20) << right_part_vector[i];
    }

    output_file << "\n";
}