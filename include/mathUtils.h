#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "../include/element.h"
#include "../include/grid.h"

class MathUtils {
   public:
    static double calculateF(Point point);

    static double calculateU(Point point);

    static bool isNumbersEqual(double first_number, double second_number,
                               double epsilon);

    static bool isPlane(std::vector<Point>& nodes, double average,
                        const Element& element, char axis);

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

    static Point calculateVectorCrossProduct(Point first_vector,
                                             Point second_vector);

    static bool isPointInPyramid(std::vector<Point>& nodes, const Point& point,
                                 const Element& element);

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
