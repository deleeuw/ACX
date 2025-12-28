#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double x[10] = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2};  // dissimilarities
double y[10] = {1, 2, 5, 4, 3, 4, 3, 2, 1, 5};  // distances
double z[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};  // indices
int n = 10;

#define EPS 1e-15

struct tuple {
    double value;
    double dist;
    int index;
};

int myComp(const void* px, const void* py) {
    double x = ((struct tuple*)px)->value;
    double y = ((struct tuple*)py)->value;
    double a = ((struct tuple*)px)->dist;
    double b = ((struct tuple*)py)->dist;
    double d = 0.0;
    if (fabs(x - y) < DBL_EPSILON * fabs(x + y)) {
        d = a - b;
    } else {
        d = x - y;
    }
    return (int)copysign(1.0, d);
}

int main() {
    struct tuple* xi = malloc(n * sizeof(struct tuple));
    for (int i = 0; i < n; i++) {
        xi[i].value = x[i];
        xi[i].dist = y[i];
        xi[i].index = z[i];
    }
    (void)qsort(xi, n, sizeof(struct tuple), myComp);
    for (int i = 0; i < n; i++) {
        x[i] = xi[i].value;
        y[i] = xi[i].dist;
        z[i] = xi[i].index;
    }
    free(xi);
    for (int i = 0; i < n; i++) {
        printf("%f %f %f\n", x[i], y[i], z[i]);
    }
    printf("\n\n");
    return EXIT_SUCCESS;
}
