#pragma once
#include "grid.h"
#include "slae.h"

#include <array>
#include <list>
#include <algorithm>
#include <set>

struct Element {

    array<int, 5> nodeIndexes;
};

class FEM {

public:
    SLAE slae;
    Grid grid;
    vector<Point> nodes;
    vector<Element> finiteElements;

    vector<vector<double>> massMatrixLocal, stiffnessMatrixLocal;
    vector<double> rightPartLocal;

    vector<int> firstBoundaryConditionNodes;

    int nodesCount, finiteElementsCount;
    
    double lambda = 1, gamma = 1;

    double UFunction(Point point);

    double GetLinearBasisFunctionTest(Point point, Point node0, Point node1, Point node2, Point node3, int numberBasisFunction);
    double GetDerivativeLinearBasisFunctionTest(Point point, Point node0, Point node1, Point node2, Point node3, int numberBasisFunction, int numberDerivativeParameter);

    void CalculateJacobian(Point node0,Point node1,Point node2,vector<vector<double>>& J);
    double CalculateDeterminant(vector<vector<double>>& matrix);
    void InverseMatrix(vector<vector<double>>& matrix, vector<vector<double>>& inv);

    double IntegrateBasisFunctions(int elementIndex, int firstBasis, int secondBasis);
    double IntegrateDerivativeBasisFunctions(int elementIndex, int firstBasis, int secondBasis);
    double IntegrateBasisFunctionForF(int elementIndex, int firstBasis);

    void GenerateLinearData(string inputFileName);
    void InputBoundaryConditions(string inputFileName);

    void GenerateMatrixPortrait();
    void GenerateLocalMatrix(int finiteElementNumber);
    void AssemblyGlobalMatrix();
    void ApplyFirstBoundaryConditions();
    void SolveFEM();
    double GetResultAtPoint(Point point);    

    void CrossProduct(Point v1, Point v2, Point v3);
    void SaveGridForVisualize();
    int GetFiniteElement(Point point);
    bool IsPointInPyramid(const Point& point, const Element& element);
    bool IsPointInsideTetrahedron(Point point, Point a, Point b, Point c, Point d);
    double CalculateTetrahedronVolume(Point p1, Point p2, Point p3, Point p4);
};