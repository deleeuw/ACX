#include "smacofSS.h"

void smacofSSMajorize(const int* nobj, const int* ndim, const int* ndat,
                      const int* itel, const int* nord, int* iind, int* jind,
                      const int* iord, const int* safe, const int* weighted,
                      double* wght, double* vinv, double* dhat, double* xold,
                      double* xnew) {
    int Nobj = *nobj, Ndim = *ndim, Itel = *itel, Nord = *nord;
    int kord = iord[(Itel - 1) % Nord], ncor = Nobj * Ndim;
    double* xone = xmalloc(Nobj * Ndim * sizeof(double));
    double* xtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* xthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* done = xmalloc(Nobj * Ndim * sizeof(double));
    double* dtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* dthr = xmalloc(Nobj * Ndim * sizeof(double));
    (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted, wght,
                                   vinv, dhat, xold, xone);
    printf("xold\n");
    (void)matrixPrint(xold, Nobj, Ndim, 6, 4);
    printf("xone\n");(void)matrixPrint(xone, Nobj, Ndim, 6, 4);
    if (kord == 0) {
        for (int i = 0; i < ncor; i++) {
            xnew[i] = xone[i];
        }
    } else {
        for (int i = 0; i < ncor; i++) {
            done[i] = xone[i] - xold[i];
        }
        if (kord == 1) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int i = 0; i < ncor; i++) {
                sum1 += xold[i] * done[i];
                sum2 += SQUARE(done[i]);
            }
            double sig1 = fabs(sum1 / sum2);
            for (int i = 0; i < ncor; i++) {
                xnew[i] = xold[i] + sig1 * done[i];
            }
        } else {
            (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind,
                                           weighted, wght, vinv, dhat, xone,
                                           xtwo);
            printf("xtwo\n");
            (void)matrixPrint(xtwo, Nobj, Ndim, 6, 4);
            for (int i = 0; i < ncor; i++) {
                dtwo[i] = xtwo[i] - 2.0 * xone[i] + xold[i];
            }
            printf("dtwo\n");
            (void)matrixPrint(dtwo, Nobj, Ndim, 6, 4);
            if (kord == 2) {
                double sum1 = 0.0, sum2 = 0.0;
                for (int i = 0; i < ncor; i++) {
                    sum1 += done[i] * dtwo[i];
                    sum2 += SQUARE(dtwo[i]);
                }
                double sig2 = fabs(sum1 / sum2);
                for (int i = 0; i < ncor; i++) {
                    xnew[i] =
                        xold[i] + 2 * sig2 * done[i] + SQUARE(sig2) * dtwo[i];
                }
                printf("xnew\n");
                (void)matrixPrint(xnew, Nobj, Ndim, 6, 4);
            } else {
                (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind,
                                               weighted, wght, vinv, dhat, xtwo,
                                               xthr);
                for (int i = 0; i < ncor; i++) {
                    dthr[i] = xthr[i] - 3.0 * xtwo[i] + 3.0 * xone[i] - xold[i];
                }
                if (kord == 3) {
                    double sum1 = 0.0, sum2 = 0.0;
                    for (int i = 0; i < ncor; i++) {
                        sum1 += dtwo[i] * dthr[i];
                        sum2 += SQUARE(dthr[i]);
                    }
                    double sig3 = fabs(sum1 / sum2);
                    for (int i = 0; i < ncor; i++) {
                        xnew[i] = xold[i] + 3.0 * sig3 * done[i] +
                                  3.0 * SQUARE(sig3) * dtwo[i] +
                                  CUBE(sig3) * dthr[i];
                    }
                }
            }
        }
    }
    xfree(xone);
    xfree(xtwo);
    xfree(xthr);
    xfree(done);
    xfree(dtwo);
    xfree(dthr);
    return;
}

void smacofSSGuttmanTransform(const int* nobj, const int* ndim, const int* ndat,
                              int* iind, int* jind, const int* weighted,
                              double* wght, double* vinv, double* dhat,
                              double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    double* edis = xmalloc(Ndat * sizeof(double));
    (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xold, edis);
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
    xfree(edis);
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