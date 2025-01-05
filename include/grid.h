#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

struct Point {

    double x, y, z;
};

class Grid {

public:
    vector<double> gridX, gridY, gridZ;

    int countX, countY, countZ;

    double stepX, stepY, stepZ;
    double scaleX, scaleY, scaleZ;
    double minimumX, minimumY, minimumZ;
    double maximumX, maximumY, maximumZ;

    Point pyramidHeight;

    void GenerateGrid(string inputFileName);
    void PrintGridInfo();
};