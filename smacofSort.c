#include "smacofSS.h"

struct quintuple {
    double value;
    double weight;
    double dist;
    int index;
    int jndex;
};

int myComp(const void* px, const void* py) {
    double x = ((struct quintuple*)px)->value;
    double y = ((struct quintuple*)py)->value;
    return (int)copysign(1.0, x - y);
}

void mySort(double* x, double* w, double* d, int* iind, int* jind,
            const int* n) {
    int nn = *n;
    struct quintuple* xi = xmalloc(nn * sizeof(struct quintuple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].weight = w[i];
        xi[i].dist = d[i];
        xi[i].index = iind[i];
        xi[i].jndex = jind[i];
    }
    (void)qsort(xi, nn, sizeof(struct quintuple), myComp);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        w[i] = xi[i].weight;
        d[i] = xi[i].dist;
        iind[i] = xi[i].index;
        jind[i] = xi[i].jndex;
    }
    xfree(xi);
    return;
}