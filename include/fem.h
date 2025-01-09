#ifndef FEM_H
#define FEM_H

const int NUMBER_OF_VERTICES_OF_PYRAMID = 5;

#include <array>
#include <string_view>
#include <tuple>
#include <vector>

#include "../include/grid.h"
#include "../include/slae.h"

struct Element final {
   private:
    std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID> node_indexes;

   public:
    Element() = default;
    Element(int vertice_index_1, int vertice_index_2, int vertice_index_3,
            int vertice_index_4, int vertice_index_5)
        : node_indexes{vertice_index_1, vertice_index_2, vertice_index_3,
                       vertice_index_4, vertice_index_5} {}

    [[nodiscard]]
    const std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID>& getNodeIndexes()
        const {
        return node_indexes;
    }

    void setNodeIndexes(
        const std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID>& new_indexes) {
        node_indexes = new_indexes;
    }
};

class FEM {
   private:
    SLAE slae;

    Grid mesh;

    std::vector<Point> nodes;

    std::vector<Element> finite_elements;

    std::vector<std::vector<double>> mass_matrix_local;
    std::vector<std::vector<double>> stiffness_matrix_local;

    std::vector<double> right_part_local;

    std::vector<std::pair<int, double>> first_boundary_condition_nodes;

    std::size_t nodes_count{0};
    std::size_t finite_elements_count{0};

    double lambda{1.0};
    double gamma{1.0};

   public:
    static double calculateF(Point point);

    static double calculateU(Point point);

    void saveTestResults(const std::vector<Point>& test_points);

    static double getLinearBasisFunctionTest(Point point, Point node_0,
                                             Point node_1, Point node_2,
                                             Point node_3,
                                             int number_basis_function);

    static double getDerivativeLinearBasisFunctionTest(
        int number_basis_function, Point point, Point node_0, Point node_1,
        Point node_2, Point node_3, int number_derivative_parameter);

    static void calculateJacobian(Point node_0, Point node_1, Point node_2,
                                  Point node_4,
                                  std::vector<std::vector<double>>& jacobian);

    static double calculateDeterminant(
        std::vector<std::vector<double>>& matrix);

    static void inverseMatrix(
        std::vector<std::vector<double>>& matrix,
        std::vector<std::vector<double>>& inversed_matrix);

    double integrateBasisFunctions(
        const std::tuple<int, int, int>& integral_parameters);

    double integrateDerivativeBasisFunctions(
        const std::tuple<int, int, int>& integral_parameters);

    double integrateBasisFunctionForF(
        const std::tuple<int, int>& integral_parameters);

    void generateLinearData(std::string_view const& input_file_name);

    void inputBoundaryConditions(std::string_view const& input_file_name);

    void generatePortrait();

    void calculateLocalComponents(int index_of_finite_element);

    void assemblyGlobalComponents();

    void solveFEM();

    double getResultAtPoint(Point point);

    static Point calculateVectorCrossProduct(Point first_vector,
                                             Point second_vector);

    void saveGridForVisualize();

    int getFiniteElementIndex(Point point);

    bool isPointInPyramid(const Point& point, const Element& element);

    static bool isPointInsideTetrahedron(Point point,
                                         Point first_tetrahedron_point,
                                         Point second_tetrahedron_point,
                                         Point third_tetrahedron_point,
                                         Point fourth_tetrahedron_point);

    static double calculateTetrahedronVolume(Point first_tetrahedron_point,
                                             Point second_tetrahedron_point,
                                             Point third_tetrahedron_point,
                                             Point fourth_tetrahedron_point);
};
#endif
