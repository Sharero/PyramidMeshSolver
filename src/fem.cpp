#include "fem.h"

double fFunction(double x, double y, double z) {

    // return 5;
    // return x + y + z;
    return -6 + pow(x, 2) + pow(y, 2) + pow(z, 2);
    // return -6 * (x + y + z) + pow(x, 3) + pow(y, 3) + pow(z, 3);
    // return sin(x * y * z) * (1 + pow(x * y, 2) + pow(x * z, 2) + pow(y * z, 2));
}
double FEM::UFunction(Point point) {

    // return 5;
    // return point.x + point.y + point.z;
    return pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2);
    // return pow(point.x, 3) + pow(point.y, 3) + pow(point.z, 3);
    // return sin(point.x * point.y * point.z);
}

void FEM::GenerateLinearData(string inputFileName) {

    nodes.resize(9);
    finiteElements.resize(6);

    nodes[0] = {0,0,0};
    nodes[1] = {1,0,0};
    nodes[2] = {0,1,0};
    nodes[3] = {1,1,0};
    nodes[4] = {0,0,1};
    nodes[5] = {1,0,1};
    nodes[6] = {0,1,1};
    nodes[7] = {1,1,1};
    nodes[8] = {0.5,0.5,0.5};

    finiteElements[0] = {0, 1, 4, 5, 8};
    finiteElements[1] = {0, 2, 4, 6, 8};
    finiteElements[2] = {1, 3, 5, 7, 8};
    finiteElements[3] = {2, 3, 6, 7, 8};
    finiteElements[4] = {0, 1, 2, 3, 8};
    finiteElements[5] = {4, 5, 6, 7, 8};

    // grid.GenerateGrid(inputFileName);

    // for (size_t j = 0; j < grid.countY; j++) {

    //     for (size_t k = 0; k < grid.countX; k++) {

    //         nodes.push_back({grid.gridX[k], grid.gridY[j], grid.gridZ[0]});
    //     }
    // }

    // for (size_t i = 1; i < grid.countZ; i++) {

    //     double scale = 1 - (grid.gridZ[i] / grid.pyramidHeight.z);

    //     for (size_t j = 0; j < grid.countY * grid.countX; j++) {
            
    //         double x = nodes[j].x + (grid.pyramidHeight.x - nodes[j].x) * (grid.gridZ[i] / grid.pyramidHeight.z);
    //         double y = nodes[j].y + (grid.pyramidHeight.y - nodes[j].y) * (grid.gridZ[i] / grid.pyramidHeight.z);

    //         nodes.push_back({x, y, grid.gridZ[i]});
    //     }   
    // }
    
    // nodes.push_back(grid.pyramidHeight);

    // int apexIndex = nodes.size() - 1;

    // for (size_t k = 0; k < grid.countZ; k++) {

    //     for (int i = 0; i < grid.countX - 1; ++i) {

    //         for (int j = 0; j < grid.countY - 1; ++j) {

    //             int idx1 = i * grid.countX + j + k * grid.countX * grid.countY;
    //             int idx2 = idx1 + 1;
    //             int idx3 = (i + 1) * grid.countX + j + k * grid.countX * grid.countY;
    //             int idx4 = idx3 + 1;

    //             finiteElements.push_back({idx1, idx2, idx3, idx4, apexIndex});
    //         }
    //     }
    // }

    nodesCount = nodes.size();
    finiteElementsCount = finiteElements.size();
    stiffnessMatrixLocal.resize(5,vector<double>(5));
    massMatrixLocal.resize(5,vector<double>(5));
    rightPartLocal.resize(5);
}

void FEM::InputBoundaryConditions(string inputFileName) {

    int countFirstConditionNodes;

    {
        ifstream inputBoundaries(inputFileName);

        inputBoundaries >> countFirstConditionNodes;

        firstBoundaryConditionNodes.resize(countFirstConditionNodes);

        for (int i = 0; i < countFirstConditionNodes; i++) {
            
            inputBoundaries >> firstBoundaryConditionNodes[i];
        }
    }
}

void FEM::GenerateMatrixPortrait() {

    int ggSize = 0;
    int tmp = 0;

    slae.n = nodesCount;

    vector<set<int>> connections(nodesCount);

    for (int i = 0; i < finiteElementsCount; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < j; k++)
                connections[finiteElements[i].nodeIndexes[k]].insert(finiteElements[i].nodeIndexes[j]);
    
    for (int i = 0; i < nodesCount; i++) {

        ggSize += connections[i].size();
    }

    slae.AllocateMemory(ggSize);

    slae.ig[0] = 0;

    for (int i = 0; i < slae.n; i++) {

        int k = 0;
        for (int j = 0; j <= i; j++) {

            if (count(connections[j].begin(), connections[j].end(), i)) {

                slae.jg[tmp] = j;
                tmp++;
                k++;
            }
        }

        slae.ig[i + 1] = slae.ig[i] + k;
    }    
}

void FEM::GenerateLocalMatrix(int finiteElementNumber) {

    for (size_t i = 0; i < 5; i++) {

        for (size_t j = 0; j < 5; j++) {

            stiffnessMatrixLocal[i][j] = lambda * IntegrateDerivativeBasisFunctions(finiteElementNumber,i,j);
            massMatrixLocal[i][j] = gamma * IntegrateBasisFunctions(finiteElementNumber,i,j);
        }    
    }

    for (size_t i = 0; i < 5; i++) {

        rightPartLocal[i] = IntegrateBasisFunctionForF(finiteElementNumber,i);
    }
}

void FEM::AssemblyGlobalMatrix() {

    for (int i = 0; i < finiteElementsCount; i++) {

        GenerateLocalMatrix(i);

        for (int j = 0; j < 5; j++) {

            for (int k = 0; k < j; k++) {

                double a = massMatrixLocal[j][k] + stiffnessMatrixLocal[j][k];

                if (finiteElements[i].nodeIndexes[j] > finiteElements[i].nodeIndexes[k])
                    slae.AddElement(finiteElements[i].nodeIndexes[j], finiteElements[i].nodeIndexes[k], a);
                else
                    slae.AddElement(finiteElements[i].nodeIndexes[k], finiteElements[i].nodeIndexes[j], a);
            }

            slae.di[finiteElements[i].nodeIndexes[j]] += massMatrixLocal[j][j] + stiffnessMatrixLocal[j][j];
            slae.f[finiteElements[i].nodeIndexes[j]] += rightPartLocal[j];
        }
    }
}

void FEM::ApplyFirstBoundaryConditions() {

    for (const auto& node : firstBoundaryConditionNodes) {

        slae.di[node] = 1.;
        slae.f[node] = UFunction(nodes[node]);

        int startIndex = slae.ig[node];
        int endIndex = slae.ig[node + 1];

        for (int i = startIndex; i < endIndex; i++) {

            slae.f[slae.jg[i]] -= slae.gg[i] * UFunction(nodes[node]);
            slae.gg[i] = 0.;
        }

        for (int p = node + 1; p < nodesCount; p++) {

            int startIndex = slae.ig[p];
            int endIndex = slae.ig[p + 1];

            for (int i = startIndex; i < endIndex; i++) {

                if (slae.jg[i] == node) {

                    slae.f[p] -= slae.gg[i] * UFunction(nodes[node]);
                    slae.gg[i] = 0.;
                }
            }
        }
    }
}

void FEM::SolveFEM() {

    ApplyFirstBoundaryConditions();
    slae.Solve();        
}

void FEM::SaveGridForVisualize() {

    ofstream outNodes("data/nodes.txt");

    for (auto& node : nodes) {

        outNodes << node.x << " " << node.y << " " << node.z << endl;
    }

    ofstream outElements("data/elements.txt");

    for (auto& element : finiteElements) {

        outElements << element.nodeIndexes[0] << " " 
                    << element.nodeIndexes[1] << " " 
                    << element.nodeIndexes[2] << " "
                    << element.nodeIndexes[3] << " "
                    << element.nodeIndexes[4] << endl;
    }
}

int FEM::GetFiniteElement(Point point) {

    int position = -1;

    for (size_t i = 0; i < finiteElements.size(); ++i) {

        if (IsPointInPyramid(point, finiteElements[i])) {

            position = i;
        }
    }

    return position;
}

void FEM::CrossProduct(Point v1, Point v2, Point v3) {
    v3.x = v1.y * v2.z - v1.z * v2.y;
    v3.y = v1.z * v2.x - v1.x * v2.z;
    v3.z = v1.x * v2.y - v1.y * v2.x;
}
bool FEM::IsPointInPyramid(const Point& point, const Element& element) {

    Point a = nodes[element.nodeIndexes[0]], 
    b = nodes[element.nodeIndexes[1]], 
    c = nodes[element.nodeIndexes[2]], 
    d = nodes[element.nodeIndexes[3]], 
    e = nodes[element.nodeIndexes[4]];

    if (IsPointInsideTetrahedron(point, a, b, c, e)) return true;

    if (IsPointInsideTetrahedron(point, a, c, d, e)) return true;

    return false;
}
double FEM::CalculateTetrahedronVolume(Point p1, Point p2, Point p3, Point p4) {
    Point v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    Point v2 = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
    Point v3 = {p4.x - p1.x, p4.y - p1.y, p4.z - p1.z};

    Point cross;
    CrossProduct(v2,v3,cross);

    double scalar = v1.x * cross.x + v1.y * cross.y + v1.z * cross.z;
    return std::abs(scalar) / 6.0;
}
bool FEM::IsPointInsideTetrahedron(Point point, Point a, Point b, Point c, Point d) {
    double volumeTotal = CalculateTetrahedronVolume(a, b, c, d);
    double v1 = CalculateTetrahedronVolume(point, b, c, d);
    double v2 = CalculateTetrahedronVolume(a, point, c, d);
    double v3 = CalculateTetrahedronVolume(a, b, point, d);
    double v4 = CalculateTetrahedronVolume(a, b, c, point);
    double sumOfVolumes = v1 + v2 + v3 + v4;
    double difference = std::abs(sumOfVolumes - volumeTotal);
    return difference < 0.00001;

}

double FEM::GetResultAtPoint(Point point) {

    double result = 0.;
    int finiteElementIndex = GetFiniteElement(point);

    for (int i = 0; i < 5; i++)
        result += slae.q[finiteElements[finiteElementIndex].nodeIndexes[i]] * 
        GetLinearBasisFunctionTest(point, nodes[finiteElements[finiteElementIndex].nodeIndexes[0]]
        , nodes[finiteElements[finiteElementIndex].nodeIndexes[1]]
        , nodes[finiteElements[finiteElementIndex].nodeIndexes[2]]
        , nodes[finiteElements[finiteElementIndex].nodeIndexes[3]], i);

    return result;
}

double FEM::GetLinearBasisFunctionTest(Point point, Point node0, Point node1, Point node2, Point node3, int numberBasisFunction) {
    
    Point vertex = {0.5,0.5,0.5};

    Point baseCenter;
    baseCenter.x = (node0.x + node1.x + node2.x + node3.x) / 4.0;
    baseCenter.y = (node0.y + node1.y + node2.y + node3.y) / 4.0;
    baseCenter.z = (node0.z + node1.z + node2.z + node3.z) / 4.0;

    Point v01 = {node1.x - node0.x, node1.y - node0.y, node1.z - node0.z};
    Point v04 = {node2.x - node0.x, node2.y - node0.y, node2.z - node0.z};
    Point v0p = {point.x,point.y,point.z};
    
    double v01norm = v01.x * v01.x + v01.y * v01.y + v01.z * v01.z;
    double v04norm = v04.x * v04.x + v04.y * v04.y + v04.z * v04.z;

    double ksi = (v0p.x * v01.x + v0p.y * v01.y + v0p.z * v01.z) / v01norm;
    double nu = (v0p.x * v04.x + v0p.y * v04.y + v0p.z * v04.z) / v04norm;

    double tx = (vertex.x - baseCenter.x), 
           ty = (vertex.y - baseCenter.y), 
           tz = (vertex.z - baseCenter.z);

    double tnorm = tx * tx + ty * ty + tz * tz;

    Point pointVector = {point.x - baseCenter.x, point.y - baseCenter.y, point.z - baseCenter.z};
    double tetta = (pointVector.x * tx + pointVector.y * ty + pointVector.z * tz) / tnorm;   

    switch(numberBasisFunction) {

        case 0:
            return (1 - ksi) * (1 - nu) * (1 - tetta);
            break;
        case 1:
            return ksi * (1 - nu) * (1 - tetta);
            break;
        case 2:
            return (1 - ksi) * nu * (1 - tetta);
            break;
        case 3:
            return ksi * nu * (1 - tetta);
            break;
        case 4:
            return tetta;
            break;
    }

    return 0.0;
}
double FEM::GetDerivativeLinearBasisFunctionTest(Point point, Point node0, Point node1, Point node2, Point node3, int numberBasisFunction, int numberDerivativeParameter) {

    Point vertex = {0.5,0.5,0.5};

    Point baseCenter;
    baseCenter.x = (node0.x + node1.x + node2.x + node3.x) / 4.0;
    baseCenter.y = (node0.y + node1.y + node2.y + node3.y) / 4.0;
    baseCenter.z = (node0.z + node1.z + node2.z + node3.z) / 4.0;

    Point v01 = {node1.x - node0.x, node1.y - node0.y, node1.z - node0.z};
    Point v04 = {node2.x - node0.x, node2.y - node0.y, node2.z - node0.z};
    Point v0p = {point.x,point.y,point.z};
    
    double v01norm = v01.x * v01.x + v01.y * v01.y + v01.z * v01.z;
    double v04norm = v04.x * v04.x + v04.y * v04.y + v04.z * v04.z;

    double ksi = (v0p.x * v01.x + v0p.y * v01.y + v0p.z * v01.z) / v01norm;
    double nu = (v0p.x * v04.x + v0p.y * v04.y + v0p.z * v04.z) / v04norm;

    double tx = (vertex.x - baseCenter.x), 
           ty = (vertex.y - baseCenter.y), 
           tz = (vertex.z - baseCenter.z);

    double tnorm = tx * tx + ty * ty + tz * tz;

    Point pointVector = {point.x - baseCenter.x, point.y - baseCenter.y, point.z - baseCenter.z};
    double tetta = (pointVector.x * tx + pointVector.y * ty + pointVector.z * tz) / tnorm;   

    vector<vector<double>> J(3, vector<double>(3));
    vector<vector<double>> JInv(3, vector<double>(3));

    CalculateJacobian(node0, node1, node2, J);
    InverseMatrix(J, JInv);
        
    double dN_dxi, dN_dnu, dN_dtetta;

    switch (numberBasisFunction) {
        case 0:
            dN_dxi = - (1 - nu) * (1 - tetta);
            dN_dnu = - (1 - ksi) * (1 - tetta);
            dN_dtetta = - (1 - ksi) * (1 - nu);
            break;
        case 1:
            dN_dxi = (1 - nu) * (1 - tetta);
            dN_dnu = - ksi * (1 - tetta);
            dN_dtetta = - ksi * (1 - nu);
            break;
        case 2:
            dN_dxi = - nu * (1 - tetta);
            dN_dnu = (1 - ksi) * (1 - tetta);
            dN_dtetta = - (1 - ksi) * nu;
            break;
        case 3:
            dN_dxi = nu * (1 - tetta);
            dN_dnu = ksi * (1 - tetta);
            dN_dtetta = - ksi * nu;
            break;
        case 4:
            dN_dxi = 0;
            dN_dnu = 0;
            dN_dtetta = 1;
            break;
        default:
            dN_dxi = 0;
            dN_dnu = 0;
            dN_dtetta = 0;
            break;
    }

    Point derivatives;

    derivatives.x =   dN_dxi * JInv[0][0] +
                      dN_dnu * JInv[1][0] +
                      dN_dtetta * JInv[2][0];
    derivatives.y =   dN_dxi * JInv[0][1] +
                      dN_dnu * JInv[1][1] +
                      dN_dtetta * JInv[2][1];
    derivatives.z =   dN_dxi * JInv[0][2] +
                      dN_dnu * JInv[1][2] +
                      dN_dtetta * JInv[2][2];

    switch (numberDerivativeParameter) {

    case 0:
        return derivatives.x;
        break;
    case 1:
        return derivatives.y;
        break;
    case 2:
        return derivatives.z;
        break;
    default:
        break;
    }

    return 0.0;
}

void FEM::CalculateJacobian(Point node0, Point node1, Point node2, vector<vector<double>>& J) {

    J[0][0] = node1.x - node0.x;
    J[0][1] = node2.x - node0.x;
    J[0][2] = nodes[8].x - node0.x;

    J[1][0] = node1.y - node0.y;
    J[1][1] = node2.y - node0.y;
    J[1][2] = nodes[8].y - node0.y;

    J[2][0] = node1.z - node0.z;
    J[2][1] = node2.z - node0.z;
    J[2][2] = nodes[8].z - node0.z;
}
double FEM::CalculateDeterminant(vector<vector<double>>& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}
void FEM::InverseMatrix(vector<vector<double>>& matrix, vector<vector<double>>& inv) {

    double det = CalculateDeterminant(matrix);

    inv[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
    inv[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det;
    inv[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;

    inv[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det;
    inv[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
    inv[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det;

    inv[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
    inv[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det;
    inv[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;
}

double FEM::IntegrateBasisFunctions(int elementIndex, int firstBasis, int secondBasis) {

    int N = 20;

    double dx = 1.0 / N;
    double dy = 1.0 / N;
    double dz = 1.0 / N;

    double integral = 0.0;

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {

                Point point;

                point.x = nodes[0].x + k * dx;
                point.y = nodes[0].y + j * dy;
                point.z = nodes[0].z + i * dz;
                double weight = 1.0;

                if (i == 0 || i == N) weight *= 0.5;
                if (j == 0 || j == N) weight *= 0.5;
                if (k == 0 || k == N) weight *= 0.5;

                integral += weight * GetLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], firstBasis) 
                                   * GetLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], secondBasis);
            }
        }
    }

    integral *= dx * dy * dz;

    return integral;
}
double FEM::IntegrateBasisFunctionForF(int elementIndex, int firstBasis) {

    int N = 20;

    double dx = 1.0 / N;
    double dy = 1.0 / N;
    double dz = 1.0 / N;

    double integral = 0.0;

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {

                Point point;

                point.x = nodes[0].x + k * dx;
                point.y = nodes[0].y + j * dy;
                point.z = nodes[0].z + i * dz;
                double weight = 1.0;

                if (i == 0 || i == N) weight *= 0.5;
                if (j == 0 || j == N) weight *= 0.5;
                if (k == 0 || k == N) weight *= 0.5;

                integral += weight * GetLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], firstBasis) 
                                   * fFunction(point.x, point.y, point.z);
            }
        }
    }

    integral *= dx * dy * dz;

    return integral;
}
double FEM::IntegrateDerivativeBasisFunctions(int elementIndex, int firstBasis, int secondBasis) {

    int N = 20;

    double dx = 1.0 / N;
    double dy = 1.0 / N;
    double dz = 1.0 / N;
    double integral = 0.0;

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {

                Point point;

                point.x = nodes[0].x + k * dx;
                point.y = nodes[0].y + j * dy;
                point.z = nodes[0].z + i * dz;
                double weight = 1.0;

                if (i == 0 || i == N) weight *= 0.5;
                if (j == 0 || j == N) weight *= 0.5;
                if (k == 0 || k == N) weight *= 0.5;

                integral += weight * (GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], firstBasis, 0)  *
                                      GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], secondBasis, 0) +  
                                      GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], firstBasis, 1)  *
                                      GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], secondBasis, 1) +  
                                      GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], firstBasis, 2)  *
                                      GetDerivativeLinearBasisFunctionTest(point, nodes[finiteElements[elementIndex].nodeIndexes[0]], nodes[finiteElements[elementIndex].nodeIndexes[1]], nodes[finiteElements[elementIndex].nodeIndexes[2]], nodes[finiteElements[elementIndex].nodeIndexes[3]], secondBasis, 2))
                                    ;
            }
        }
    }

    integral *= dx * dy * dz;

    return integral;
}
