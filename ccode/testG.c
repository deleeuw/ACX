#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SQUARE(x) ((x) * (x))

void smacofSSGuttmanTransform(const int* nobj, const int* ndim, const int* ndat,
                              int* iind, int* jind, const int* weighted,
                              double* wght, double* vinv, double* dhat,
                              double* xold, double* xnew);
void smacofSSDistances(const int* nobj, const int* ndim, const int* ndat,
                       int* iind, int* jind, double* xmat, double* edis);
void matrixPrint(const double* x, const size_t Nrow, const size_t Ncol,
                 const int digits, const int width);

void smacofSSGuttmanTransform(const int* nobj, const int* ndim, const int* ndat,
                              int* iind, int* jind, const int* weighted,
                              double* wght, double* vinv, double* dhat,
                              double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim, Weighted = *weighted;
    double* xtmp = malloc(Nobj * Ndim * sizeof(double));
    double* edis = malloc(Ndat * sizeof(double));
    (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xold, edis);
    for (int k = 0; k < Ndat; k++) {
    	printf(" %6.2f ", edis[k]);
    }
    printf("\n");
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        double ecof = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            int iobj = i + Nobj * s, jobj = j + Nobj * s;
            double add = ecof * (xold[iobj] - xold[jobj]);
            xtmp[iobj] += add;
            xtmp[jobj] -= add;
        }
    }
    for (int k = 0; k < Nobj * Ndim; k++) {
        if (Weighted) {
            xnew[k] = 0.0;
        } else {
            xnew[k] = xtmp[k];
        }
    }
    if (Weighted) {
        for (int k = 0; k < Ndat; k++) {
            int i = iind[k];
            int j = jind[k];
            double ecof = -vinv[k];
            for (int s = 0; s < Ndim; s++) {
                int iobj = i + Nobj * s, jobj = j + Nobj * s;
                double add = ecof * (xold[iobj] - xold[jobj]);
                xnew[iobj] += add;
                xnew[jobj] -= add;
            }
        }
    } else {
        for (int i = 0; i < Nobj; i++) {
            int is = i;
            for (int s = 0; s < Ndim; s++) {
                xnew[is] /= (double)Nobj;
                is += Nobj;
            }
        }
    }
    free(xtmp);
    free(edis);
    return;
}

void smacofSSDistances(const int* nobj, const int* ndim, const int* ndat,
                       int* iind, int* jind, double* xmat, double* edis) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k], js = jind[k];
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(xmat[is] - xmat[js]);
            is += Nobj;
            js += Nobj;
        }
        edis[k] = sqrt(sum);
    }
}

void matrixPrint(const double* x, const size_t Nrow, const size_t Ncol,
                 const int digits, const int width) {
    size_t k = 0;
    for (size_t i = 0; i < Nrow; i++) {
        k = i;
        for (size_t s = 0; s < Ncol; s++) {
            printf(" %+*.*f ", width, digits, x[k]);
            k += Nrow;
        }
        printf("\n");
    }
    printf("\n\n");
}

int Nobj = 4, Ndim = 2, Ndat = 6, Weighted = 0;
int iind[6] = {2, 3, 4, 3, 4, 4};
int jind[6] = {1, 1, 1, 2, 2, 3};
double wght[6] = {1.0};
double dhat[6] = {1.0};
double vinv[6] = {1.0};
double xold[8] = {-1.0, 2.0, -3.0, 4.0, -1.0, -1.0, 1.0, 1.0};
double xnew[8] = {0.0};

int main(void) {
    (void)matrixPrint(xold, Nobj, Ndim, 6, 4);
    (void)smacofSSGuttmanTransform(&Nobj, &Ndim, &Ndat, iind, jind, &Weighted,
                                   wght, vinv, dhat, xold, xnew);
    (void)matrixPrint(xnew, Nobj, Ndim, 6, 4);
    return (EXIT_SUCCESS);
}