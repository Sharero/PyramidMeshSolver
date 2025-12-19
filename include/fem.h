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

    BASIS_TYPE basis_functions_type;
    BASIS_ELEMENT_TYPE basis_functions_element_type;

    double lambda{1.0};
    double gamma{1.0};

    [[nodiscard]]
    int findFaceNodeIndex(int node_1, int node_2, double divisor) const;

    void generateBaseNodesCombinations(
        int combination_size, std::vector<Element>& base_nodes_combinations);

    void generateNodes(int basis_type);

    void generateElements(int p);

    double getResultAtPoint(Point point);

    void writeGridInformationToFile();

    int getFiniteElementIndex(Point point);

    void generatePortrait();

    void calculateLocalComponents(int index_of_finite_element);

    void writeLocalComponentsToFiles(
        int index_of_finite_element,
        const std::vector<std::vector<double>>& stiffness_matrix,
        const std::vector<std::vector<double>>& mass_matrix,
        const std::vector<double>& right_part_vector);

    void assemblyGlobalComponents();

   public:
    FEM(BASIS_TYPE basis_type, BASIS_ELEMENT_TYPE basis_element_type)
        : basis_functions_type(basis_type),
          basis_functions_element_type(basis_element_type) {}

    void checkBasisFunctionsDeltaProperty();

    void calculateResultAtPointsToFile(const std::vector<Point>& test_points);

    void generateFEMData(std::string_view const& input_file_name);

    void inputBoundaryConditions(std::string_view const& input_file_name);

    void solveFEM();
};
#endif
