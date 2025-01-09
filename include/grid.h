#ifndef GRID_H
#define GRID_H

#include <filesystem>
#include <vector>

struct Point {
    double x;
    double y;
    double z;
};

class Grid {
   private:
    std::vector<double> grid_x;
    std::vector<double> grid_y;
    std::vector<double> grid_z;

    int count_x_points{0}, count_y_points{0}, count_z_points{0};
    static int count_x_points1;
    double grid_x_step{0.0}, grid_y_step{0.0}, grid_z_step{0.0};

    double grid_x_scale{0.0}, grid_y_scale{0.0}, grid_z_scale{0.0};

    double minimum_x_coordinate{0.0}, minimum_y_coordinate{0.0},
        minimum_z_coordinate{0.0};

    double maximum_x_coordinate{0.0}, maximum_y_coordinate{0.0},
        maximum_z_coordinate{0.0};

    Point pyramidHeight{0.0, 0.0, 0.0};

   public:
    void generateGrid(const std::filesystem::path& input_file_name);
    void printGridInfo();
};
#endif
