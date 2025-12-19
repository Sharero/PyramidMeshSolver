#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <set>
#include <tuple>
#include <vector>

#include "../include/element.h"
#include "../include/grid.h"

class MathUtils {
   public:
    static void calculateCombinations(
        int combination_size, const std::set<int>& required_indexes,
        const std::vector<int>& remaining_indexes,
        std::vector<Element>& base_nodes_combinations);

    static void filterRemainingIndexesForBaseCombinations(
        int vertex_index, const std::set<int>& forbidden_indexes,
        const std::set<int>& required_indexes,
        std::vector<int>& remaining_indexes);

    static bool isCubicNodeIllegal(
        const std::tuple<int, int, int>& cubic_node_parameters);

    static double calculateF(Point point);

    static double calculateU(Point point);

    static bool isNumbersEqual(double first_number, double second_number,
                               double epsilon);

    static bool isPlane(std::vector<Point>& nodes, double average,
                        const Element& element, char axis);

    static Point calculateVectorCrossProduct(Point first_vector,
                                             Point second_vector);

    static void calculateLocalCoordinates(Point physical_point, Point node_0,
                                          Point node_1, Point node_2,
                                          Point node_3, Point vertex,
                                          double& ksi_parameter,
                                          double& nu_parameter,
                                          double& tetta_parameter);

    static void calculatePhysicalCoordinates(Point local_point, Point node_0,
                                             Point node_1, Point node_2,
                                             Point node_3, Point vertex,
                                             double& x, double& y, double& z);

    static double getBasisFunction(
        const std::vector<Point>& nodes, Point point,
        const std::vector<int>& node_indexes, BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_elements_type,
        int number_basis_function, bool need_coordinate_transform_to_local);

    static double getDerivativeBasisFunction(
        const std::vector<Point>& nodes, Point point,
        const std::vector<int>& node_indexes, BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_elements_type,
        int number_of_first_basis_function, int number_of_second_basis_function,
        bool need_coordinate_transform_to_local);

    static double getLinearBasisFunction(
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, int number_basis_function,
        bool need_coordinate_transform_to_local);

    static double getQuadraticBasisFunction(
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, int number_basis_function,
        bool need_coordinate_transform_to_local);

    static double getCubicBasisFunction(
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, int number_basis_function,
        bool need_coordinate_transform_to_local);

    static double getHierarhicalQuadraticBasisFunction(
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, int number_basis_function,
        bool need_coordinate_transform_to_local);

    static double getHierarhicalCubicBasisFunction(
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, int number_basis_function,
        bool need_coordinate_transform_to_local);

    static double getDerivativeLinearBasisFunction(
        int number_of_first_basis_function, int number_of_second_basis_function,
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, bool need_coordinate_transform_to_local);

    static double getDerivativeQuadraticBasisFunction(
        int number_of_first_basis_function, int number_of_second_basis_function,
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, bool need_coordinate_transform_to_local);

    static double getDerivativeCubicBasisFunction(
        int number_of_first_basis_function, int number_of_second_basis_function,
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, bool need_coordinate_transform_to_local);

    static double getDerivativeHierarhicalQuadraticBasisFunction(
        int number_of_first_basis_function, int number_of_second_basis_function,
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, bool need_coordinate_transform_to_local);

    static double getDerivativeHierarhicalCubicBasisFunction(
        int number_of_first_basis_function, int number_of_second_basis_function,
        Point point, Point node_0, Point node_1, Point node_2, Point node_3,
        Point vertex, bool need_coordinate_transform_to_local);

    static void getMainPyramideNodes(
        const std::vector<int>& node_indexes, const std::vector<Point>& nodes,
        BASIS_TYPE basis_functions_type,
        BASIS_ELEMENT_TYPE basis_functions_elements_type, Point& node_0,
        Point& node_1, Point& node_2, Point& node_3, Point& node_4);

    static void multiplyMatrixByMatrix(
        std::vector<std::vector<double>>& first_matrix,
        std::vector<std::vector<double>>& second_matrix,
        std::vector<std::vector<double>>& result);

    static double calculateJacobian(std::vector<std::vector<double>>& J);

    static void calculateJacobiMatrix(Point node_0, Point node_1, Point node_2,
                                      Point node_3, Point node_4,
                                      std::vector<std::vector<double>>& J);

    static void transponateMatrix(std::vector<std::vector<double>>& J,
                                  std::vector<std::vector<double>>& JT);

    static void calculateInvertedMatrix(
        std::vector<std::vector<double>>& J,
        std::vector<std::vector<double>>& JInverted);
};
#endif
