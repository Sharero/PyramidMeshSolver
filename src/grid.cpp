#include "../include/grid.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

void Grid::generateGrid(const std::filesystem::path& input_file_name) {
    {
        std::ifstream input_file(input_file_name, std::ios::in);

        input_file >> minimum_x_coordinate >> maximum_x_coordinate >>
            minimum_y_coordinate >> maximum_y_coordinate >>
            minimum_z_coordinate >> maximum_z_coordinate;

        input_file >> grid_x_step >> grid_y_step >> grid_z_step;

        input_file >> grid_x_scale >> grid_y_scale >> grid_z_scale;
    }

    count_x_points = static_cast<int>(
        (maximum_x_coordinate - minimum_x_coordinate) / grid_x_step + 1);

    count_y_points = static_cast<int>(
        (maximum_y_coordinate - minimum_y_coordinate) / grid_y_step + 1);

    count_z_points = static_cast<int>(
        (maximum_z_coordinate - minimum_z_coordinate) / grid_z_step + 1);

    grid_x.resize(count_x_points);
    grid_y.resize(count_y_points);
    grid_z.resize(count_z_points);

    const double step_x = grid_x_step;
    const double step_y = grid_y_step;
    const double step_z = grid_z_step;

    const double minimum_x = minimum_x_coordinate;
    const double minimum_y = minimum_y_coordinate;
    const double minimum_z = minimum_z_coordinate;

    std::transform(
        grid_x.begin(), grid_x.end(), grid_x.begin(),
        [step_x, minimum_x, current_index = 0](double& grid_x_value) mutable {
            grid_x_value += current_index * step_x + minimum_x;
            ++current_index;
            return grid_x_value;
        });

    std::transform(
        grid_y.begin(), grid_y.end(), grid_y.begin(),
        [step_y, minimum_y, current_index = 0](double& grid_y_value) mutable {
            grid_y_value += current_index * step_y + minimum_y;
            ++current_index;
            return grid_y_value;
        });

    std::transform(
        grid_z.begin(), grid_z.end(), grid_z.begin(),
        [step_z, minimum_z, current_index = 0](double& grid_z_value) mutable {
            grid_z_value += current_index * step_z + minimum_z;
            ++current_index;
            return grid_z_value;
        });
}

void Grid::printGridData() {
    std::cout << "stepX: " << grid_x_step << " ";
    std::cout << "scaleX: " << grid_x_scale << '\n';
    std::cout << "grid_x: ";

#pragma unroll 4
    for (double const grid_x_element : grid_x) {
        std::cout << grid_x_element << " ";
    }

    std::cout << '\n' << '\n';
    std::cout << "stepY: " << grid_y_step << " ";
    std::cout << "scaleY: " << grid_y_scale << '\n';
    std::cout << "grid_y: ";

#pragma unroll 4
    for (double const grid_y_element : grid_y) {
        std::cout << grid_y_element << " ";
    }

    std::cout << '\n' << '\n';
    std::cout << "stepZ: " << grid_z_step << " ";
    std::cout << "scaleZ: " << grid_z_scale << '\n';
    std::cout << "grid_z: ";

#pragma unroll 4
    for (double const grid_z_element : grid_z) {
        std::cout << grid_z_element << " ";
    }

    std::cout << '\n';
}
