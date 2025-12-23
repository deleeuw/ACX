#include "smacofSS.h"

void smacofSSMajorize(int* nobj, int* ndim, int* ndat, int* iind, int* jind, int *iord,
                      int* weighted, double* wght, double* vinv, double* edis,
                      double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim, Iord = *iord;
    double sig = 0.0, sum1 = 0.0, sum2 = 0.0;
    double* xone = xmalloc(Nobj * Ndim * sizeof(double));
    double* xtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* xthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* done = xmalloc(Nobj * Ndim * sizeof(double));
    double* dtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* dthr = xmalloc(Nobj * Ndim * sizeof(double));
    (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted, wght,
                                   vinv, edis, dhat, xold, xone);
    for (int i = 0; i < Nobj * Ndim; i++) {
        done[i] = xone[i] - xold[i];
    }
    (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xone, edis);
    (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted, wght,
                                   vinv, edis, dhat, xone, xtwo);
    for (int i = 0; i < Nobj * Ndim; i++) {
        dtwo[i] = xtwo[i] - 2.0 * xone[i] + xold[i];
    }
    if (Iord == 0) {
        for (int i = 0; i < Nobj * Ndim; i++) {
            xnew[i] = xone[i];
        }
    }
    if (Iord == 1) {
        sum1 = 0.0;
        sum2 = 0.0;
        for (int i = 0; i < Nobj * Ndim; i++) {
            sum1 += xold[i] * done[i];
            sum2 += SQUARE(done[i]);
        }
        sig = sum1 / sum2;
        for (int i = 0; i < Nobj * Ndim; i++) {
            xnew[i] = xold[i] + sig * done[i];
        }
    }
    if (Iord == 2) {
        sum1 = 0.0;
        sum2 = 0.0;
        for (int i = 0; i < Nobj * Ndim; i++) {
            sum1 += done[i] * dtwo[i];
            sum2 += SQUARE(dtwo[i]);
        }
        sig = sum1 / sum2;
        for (int i = 0; i < Nobj * Ndim; i++) {
            xnew[i] = xold[i] + 2 * sig * done[i] + SQUARE(sig) * dtwo[i];
        }
    }
    if (Iord == 3) {
        (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted,
                                       wght, vinv, edis, dhat, xtwo, xthr);
        sum1 = 0.0;
        sum2 = 0.0;
        for (int i = 0; i < Nobj * Ndim; i++) {
            dthr[i] = xthr[i] - 3.0 * xtwo[i] + 3.0 * xone[i] - xold[i];
            sum1 += dtwo[i] * dthr[i];
            sum2 += SQUARE(dthr[i]);
        }
        sig = sum1 / sum2;
        for (int i = 0; i < Nobj * Ndim; i++) {
            xnew[i] = xold[i] + 3.0 * sig * done[i] +
                      3.0 * SQUARE(sig) * dtwo[i] + CUBE(sig) * dthr[i];
        }
    }
    xfree(xone);
    xfree(xtwo);
    xfree(xthr);
    xfree(done);
    xfree(dtwo);
    xfree(dthr);
    (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xnew, edis);
    return;
}

void smacofSSGuttmanTransform(int* nobj, int* ndim, int* ndat, int* iind,
                              int* jind, int* weighted, double* wght,
                              double* vinv, double* edis, double* dhat,
                              double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        if (edis[k] == 0.0) {
            continue;
        }
        int is = iind[k], js = jind[k];
        double elem = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
            is += Nobj;
            js += Nobj;
        }
    }
    for (int k = 0; k < Nobj * Ndim; k++) {
        if (weighted) {
            xnew[k] = 0.0;
        } else {
            xnew[k] = xtmp[k];
        }
    }
    if (weighted) {
        int k = 0;
        for (int j = 0; j < Nobj - 1; j++) {
            for (int i = j + 1; i < Nobj; i++) {
                double elem = vinv[k];
                int is = i, js = j;
                for (int s = 0; s < Ndim; s++) {
                    double add = elem * (xtmp[is] - xtmp[js]);
                    xnew[is] += add;
                    xnew[js] -= add;
                    is += Nobj;
                    js += Nobj;
                }
                k++;
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
    xfree(xtmp);
    return;
}

void smacofSSDistances(int* nobj, int* ndim, int* ndat, int* iind, int* jind,
                       double* xmat, double* edis) {
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