#include "slae.h"

void SLAE::AllocateMemory(int size) {

    ig.resize(n + 1);
    jg.resize(size);
    di.resize(n);
    gg.resize(size);
    f.resize(n);
    q.resize(n);
}

void SLAE::AddElement(int i, int j, double element) {

    bool flag = false;
    int index = 0;

    for (int k = ig[i]; k < ig[i + 1] && !flag; k++) {

        if (jg[k] == j) {

            index = k;
            flag = true;
        }
    }

    gg[index] += element;
}

void SLAE::Solve() {

    msg.Initialize(ig, jg, di, gg, f, n);
    msg.Calculate(q);
}
