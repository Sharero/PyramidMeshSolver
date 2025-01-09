#include "../include/msg.h"

#include <cmath>
#include <cstddef>
#include <vector>

void MSG::multiplyMatrixAndVector(const std::vector<double>& f_vector,
                                  std::vector<double>& x_vector) {
    for (int i = 0; i < slae_size; i++) {
        x_vector[i] = di[i] * f_vector[i];

#pragma unroll 4
        for (int k = gi[i], k1 = gi[i + 1]; k < k1; k++) {
            x_vector[i] += gg[k] * f_vector[gj[k]];
            x_vector[gj[k]] += gg[k] * f_vector[i];
        }
    }
}

void MSG::calculateLLT(const std::vector<double>& f_vector,
                       std::vector<double>& x_vector) {
    calculateL(f_vector, x_vector);
    calculateLT(x_vector, x_vector);
}

void MSG::calculateLT(const std::vector<double>& f_vector,
                      std::vector<double>& x_vector) {
    std::vector<double> f_vector_copy = f_vector;

    for (std::size_t k = slae_size, k1 = slae_size - 1; k > 0; k--, k1--) {
        x_vector[k1] = f_vector_copy[k1] / di_of_lower_triangle[k1];

#pragma unroll 4
        for (int i = gi[k1]; i < gi[k]; i++) {
            f_vector_copy[gj[i]] -= gg_of_lower_triangle[i] * x_vector[k1];
        }
    }
}

void MSG::calculateL(const std::vector<double>& f_vector,
                     std::vector<double>& x_vector) {
    for (int k = 1, k1 = 0; k <= slae_size; k++, k1++) {
        double sum = 0.0;

#pragma unroll 4
        for (int i = gi[k1]; i < gi[k]; i++) {
            sum += gg_of_lower_triangle[i] * x_vector[gj[i]];
        }

        x_vector[k1] = (f_vector[k1] - sum) / di_of_lower_triangle[k1];
    }
}

void MSG::makeLLTDecomposition() {
    double diagonal_sum = 0.0;
    double lower_triangle_sum = 0.0;

    for (int k = 0; k < slae_size; k++) {
        diagonal_sum = 0;

        const int row_start_index = gi[k];
        const int row_end_index = gi[k + 1];

        for (int i = row_start_index; i < row_end_index; i++) {
            lower_triangle_sum = 0;

            int column_start_index = gi[gj[i]];
            const int column_end_index = gi[gj[i] + 1];

            for (int current_row_index = row_start_index; current_row_index < i;
                 current_row_index++) {
#pragma unroll 4
                for (int j = column_start_index; j < column_end_index; j++) {
                    if (gj[current_row_index] == gj[j]) {
                        lower_triangle_sum +=
                            gg_of_lower_triangle[current_row_index] *
                            gg_of_lower_triangle[j];
                        column_start_index++;
                    }
                }
            }

            gg_of_lower_triangle[i] =
                (gg_of_lower_triangle[i] - lower_triangle_sum) /
                di_of_lower_triangle[gj[i]];

            diagonal_sum += gg_of_lower_triangle[i] * gg_of_lower_triangle[i];
        }

        di_of_lower_triangle[k] = sqrt(di_of_lower_triangle[k] - diagonal_sum);
    }
}

void MSG::calculateSLAE(std::vector<double>& solution) {
    int const maximum_iterations = 1000;
    double const epsilon = 1E-15;

    bool end = false;

    double discrepansy = 0.0;
    double rp_norm = 0.0;
    double scalar1 = 0.0;
    double scalar2 = 0.0;
    double betta = 0.0;
    double alpha = 0.0;

    multiplyMatrixAndVector(x0, r);

#pragma unroll 4
    for (int i = 0; i < slae_size; i++) {
        r[i] = rp[i] - r[i];
    }

    makeLLTDecomposition();
    calculateLLT(r, z);

#pragma unroll 4
    for (int i = 0; i < slae_size; i++) {
        p[i] = z[i];
    }

    rp_norm = sqrt(calculateScalarProduct(rp, rp));
    scalar1 = calculateScalarProduct(p, r);

    for (int current_iteration = 0;
         current_iteration < maximum_iterations && !end; current_iteration++) {
        discrepansy = sqrt(calculateScalarProduct(r, r));

        if (epsilon < discrepansy / rp_norm) {
            multiplyMatrixAndVector(z, s);

            alpha = scalar1 / calculateScalarProduct(s, z);

#pragma unroll 4
            for (int i = 0; i < slae_size; i++) {
                x0[i] += alpha * z[i];
                r[i] -= alpha * s[i];
            }

            calculateLLT(r, p);
            scalar2 = calculateScalarProduct(p, r);

            betta = scalar2 / scalar1;

            scalar1 = scalar2;

#pragma unroll 4
            for (int i = 0; i < slae_size; i++) {
                z[i] = p[i] + betta * z[i];
            }
        } else {
            end = true;
        }
    }

#pragma unroll 4
    for (int i = 0; i < slae_size; i++) {
        solution[i] = x0[i];
    }
}

double MSG::calculateScalarProduct(
    const std::vector<double>& first_vector,
    const std::vector<double>& second_vector) const {
    double scalar = 0.0;

#pragma unroll 4
    for (int i = 0; i < slae_size; i++) {
        scalar += first_vector[i] * second_vector[i];
    }

    return scalar;
}