#include "grid.h"

void Grid::GenerateGrid(string inputFileName) {

    {
        ifstream inputFile(inputFileName);

        inputFile >> minimumX >> maximumX >>
                     minimumY >> maximumY >>
                     minimumZ >> maximumZ;

        inputFile >> stepX >> stepY >> stepZ;
        inputFile >> scaleX >> scaleY >> scaleZ;
    }

    countX = (maximumX - minimumX) / stepX + 1;
    countY = (maximumY - minimumY) / stepY + 1;
    countZ = (maximumZ - minimumZ) / stepZ;

    pyramidHeight = {(maximumX + minimumX) / 2.0, (maximumY + minimumY) / 2.0, maximumZ};

    gridX.resize(countX);
    gridY.resize(countY);
    gridZ.resize(countZ);

    for (size_t i = 0; i < countX; i++) {

        gridX[i] = minimumX + i * stepX;
    }

    for (size_t i = 0; i < countY; i++) {

        gridY[i] = minimumY + i * stepY;
    }

    for (size_t i = 0; i < countZ; i++) {

        gridZ[i] = minimumZ + i * stepZ;
    }
}

void Grid::PrintGridInfo() {


    cout << "stepX: " << stepX << " ";
    cout << "scaleX: " << scaleX << endl;
    cout << "gridX: ";
    for (size_t i = 0; i < gridX.size(); i++) {

        cout << gridX[i] << " ";
    }

    cout << endl << endl;

    cout << "stepY: " << stepY << " ";
    cout << "scaleY: " << scaleY << endl;
    cout << "gridY: ";
    for (size_t i = 0; i < gridY.size(); i++) {

        cout << gridY[i] << " ";
    }

    cout << endl << endl;

    cout << "stepZ: " << stepZ << " ";
    cout << "scaleZ: " << scaleZ << endl;
    cout << "gridZ: ";
    for (size_t i = 0; i < gridZ.size(); i++) {

        cout << gridZ[i] << " ";
    }

    cout << endl;
}
