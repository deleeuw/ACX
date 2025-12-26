#include "smacofSS.h"

void smacofMPInverseV(const int* nobj, const int* ndat, int* iind, int* jind,
                      double* wght, double* vinv) {
    int Nobj = *nobj, Ndat = *ndat;
    double** v = xcalloc(Nobj, sizeof(double*));
    for (int i = 0; i < Nobj; i++) {
        v[i] = xcalloc(Nobj, sizeof(double));
    }
    for (int k = 0; k < Ndat; k++) {
        int ii = iind[k];
        int jj = jind[k];
        v[ii][jj] = -wght[k];
        v[jj][ii] = -wght[k];
    }
    for (int i = 0; i < Nobj; i++) {
        double sum = 0.0;
        for (int j = 0; j < Nobj; j++) {
            sum -= v[i][j];
        }
        v[i][i] = sum;
    }
    for (int i = 0; i < Nobj; i++) {
        for (int j = 0; j < Nobj; j++) {
            v[i][j] += 1 / (double)Nobj;
        }
    }
    for (int k = 0; k < Nobj; k++) {
        double piv = v[k][k];
        for (int i = 0; i < Nobj; i++) {
            if (i == k) {
                continue;
            }
            for (int j = 0; j < Nobj; j++) {
                if (j == k) {
                    continue;
                }
                v[i][j] = v[i][j] - v[i][k] * v[k][j] / piv;
            }
        }
        for (int i = 0; i < Nobj; i++) {
            if (i == k) {
                continue;
            }
            v[i][k] = v[i][k] / piv;
        }
        for (int j = 0; j < Nobj; j++) {
            if (j == k) {
                continue;
            }
            v[k][j] = v[k][j] / piv;
        }
        v[k][k] = -1.0 / piv;
    }
    int k = 0;
    for (int j = 0; j < Nobj - 1; j++) {
        for (int i = j + 1; i < Nobj; i++) {
            vinv[k] = -v[i][j] - (1.0 / (double)Nobj);
            k++;
        }
    }
    for (int i = 0; i < Nobj; i++) {
        free(v[i]);
    }
    free(v);
    return;
}