#include "../include/slae.h"

#include <algorithm>
#include <cstddef>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "../include/msg.h"

void SLAE::allocateMemory(std::size_t size) {
    ig.resize(slae_size + 1);
    jg.resize(size);
    di.resize(slae_size);
    gg.resize(size);
    f.resize(slae_size);
    q.resize(slae_size);
}

void SLAE::addLowerTriangleElement(
    const std::tuple<int, int, double>& element_parameters) {
    int const row_index = std::get<0>(element_parameters);
    int const column_index = std::get<1>(element_parameters);
    double const value = std::get<2>(element_parameters);

    bool flag = false;
    int index = 0;

#pragma unroll 4
    for (int k = ig[row_index]; k < ig[row_index + 1] && !flag; k++) {
        if (jg[k] == column_index) {
            index = k;
            flag = true;
        }
    }

    gg[index] += value;
}

void SLAE::addDiagonalElement(int diagonal_element_index,
                              double diagonal_element_value) {
    di[diagonal_element_index] += diagonal_element_value;
}

void SLAE::addRightPartElement(int right_part_element_index,
                               double right_part_element_value) {
    f[right_part_element_index] += right_part_element_value;
}

void SLAE::solveSLAE() {
    MSG msg(ig, di, jg, gg, slae_size, f);
    msg.calculateSLAE(q);
}

void SLAE::applyFirstBoundaryConditions(
    const std::vector<std::pair<int, double>>& first_boundary_condition_nodes) {
    for (const auto& node : first_boundary_condition_nodes) {
        di[node.first] = 1.;
        f[node.first] = node.second;

#pragma unroll 4
        for (int i = ig[node.first]; i < ig[node.first + 1]; i++) {
            f[jg[i]] -= gg[i] * node.second;
            gg[i] = 0.;
        }

        for (int j = node.first + 1; j < slae_size; j++) {
#pragma unroll 4
            for (int i = ig[j]; i < ig[j + 1]; i++) {
                if (jg[i] == node.first) {
                    f[j] -= gg[i] * node.second;
                    gg[i] = 0.;
                }
            }
        }
    }
}

void SLAE::setIndexArrays(std::vector<std::set<int>>& connections) {
    int current_column_index = 0;

    ig[0] = 0;

    for (int i = 0; i < slae_size; i++) {
        int current_count_of_elements = 0;

#pragma unroll 4
        for (int j = 0; j <= i; j++) {
            if (count(connections[j].begin(), connections[j].end(), i) != 0) {
                jg[current_column_index] = j;
                current_column_index++;
                current_count_of_elements++;
            }
        }

        ig[i + 1] = ig[i] + current_count_of_elements;
    }
}