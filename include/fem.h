#ifndef FEM_H
#define FEM_H

#include <string_view>
#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"
#include "../include/integrator.h"
#include "../include/mathUtils.h"
#include "../include/slae.h"

class FEM {
   private:
    SLAE slae;

    Integrator integrator;

    Grid mesh;

    std::vector<Point> nodes;

    std::vector<Element> finite_elements;

    std::vector<std::vector<double>> mass_matrix_local;
    std::vector<std::vector<double>> stiffness_matrix_local;

    std::vector<double> right_part_local;

    std::vector<std::pair<int, double>> first_boundary_condition_nodes;

    std::size_t nodes_count{0};
    std::size_t finite_elements_count{0};
    std::size_t number_of_vertices_of_pyramid{0};

    double lambda{1.0};
    double gamma{1.0};

    void generateLinearNodes();

    void generateQuadraticNodes();

    void generateCubicNodes();

    double getResultAtPoint(Point point);

    void saveGridForVisualize();

    int getFiniteElementIndex(Point point);

    void generatePortrait();

    void calculateLocalComponents(int index_of_finite_element);

    void assemblyGlobalComponents();

   public:
    void checkBasisFunctionsToEqualsOne(
        BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_elements_type);

    void saveTestResults(const std::vector<Point>& test_points);

    void generateFEMData(std::string_view const& input_file_name,
                         BASIS_TYPE basis_functions_type,
                         BASIS_ELEMENT_TYPE basis_functions_elements_type);

    void inputBoundaryConditions(std::string_view const& input_file_name);

    void solveFEM();
};
#endif
