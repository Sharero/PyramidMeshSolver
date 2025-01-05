#include "fem.h"

using namespace std;

int main(void) {
    FEM fem;

    int dota_jopa = 2;

    Point point1 = {5E-1, 5E-1, 0E-1};
    Point point2 = {0E-1, 5E-1, 0E-1};
    Point point3 = {5E-1, 5E-1, 10E-1};
    Point point4 = {5E-1, 10E-1, 75E-2};

    fem.GenerateLinearData("data/gridInfo.txt");
    fem.InputBoundaryConditions("data/boundaries.txt");
    fem.SaveGridForVisualize();
    fem.GenerateMatrixPortrait();
    fem.AssemblyGlobalMatrix();
    fem.SolveFEM();

    // cout << scientific << fem.GetResultAtPoint(point1, 4) << " " <<
    // fem.UFunction(point1) << " " << abs(fem.GetResultAtPoint(point1, 4) -
    // fem.UFunction(point1)) << endl; cout << scientific <<
    // fem.GetResultAtPoint(point2, 1) << " " << fem.UFunction(point2) << " " <<
    // abs(fem.GetResultAtPoint(point2, 1) - fem.UFunction(point2)) << endl;
    // cout
    // << scientific << fem.GetResultAtPoint(point3, 5) << " " <<
    // fem.UFunction(point3) << " " << abs(fem.GetResultAtPoint(point3, 5) -
    // fem.UFunction(point3)) << endl; cout << scientific <<
    // fem.GetResultAtPoint(point4, 3) << " " << fem.UFunction(point4) << " " <<
    // abs(fem.GetResultAtPoint(point4, 3) - fem.UFunction(point4)) << endl;
}