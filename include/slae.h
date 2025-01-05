#include <vector>
#include "msg.h"

using namespace std;

class SLAE {

public:
    int n;

    vector<int> ig, jg;
    vector<double> di, gg, f, q;

    MSG msg;

    void AllocateMemory(int size);
    void AddElement(int i, int j, double element);
    void Solve();
};