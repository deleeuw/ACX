#include "smacofSS.h"
#define DEBUG 0

void smacofSSMajorize(const int* nobj, const int* ndim, const int* ndat,
                      const int* itel, const int* nord, int* iind, int* jind,
                      const int* iord, const int* safe, const int* weighted,
                      double* wght, double* vinv, double* dhat, double* xold,
                      double* xnew) {
    int Nobj = *nobj, Ndim = *ndim, Itel = *itel, Nord = *nord, Ndat = *ndat;
    int kord = iord[(Itel - 1) % Nord], ncor = Nobj * Ndim, Safe = *safe;
    double sone = 0.0, snew = 0.0;
    double* xone = xmalloc(Nobj * Ndim * sizeof(double));
    double* xtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* xthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* done = xmalloc(Nobj * Ndim * sizeof(double));
    double* dtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* dthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* edis = xmalloc(Ndat * sizeof(double));
    (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted, wght,
                                   vinv, dhat, xold, xone);
    if (Safe) {
        (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xnew, edis);
        sone = smacofSSLoss(ndat, edis, dhat, wght);
    }
    if (DEBUG) {
        printf("xold\n");
        (void)matrixPrint(xold, Nobj, Ndim, 6, 4);
        printf("xone\n");
        (void)matrixPrint(xone, Nobj, Ndim, 6, 4);
    }
    if (kord == 0) {
        for (int i = 0; i < ncor; i++) {
            xnew[i] = xone[i];
        }
    } else {
        for (int i = 0; i < ncor; i++) {
            done[i] = xone[i] - xold[i];
        }
        if (DEBUG) {
            printf("done\n");
            (void)matrixPrint(done, Nobj, Ndim, 6, 4);
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
            if (DEBUG) {
                printf("xtwo\n");
                (void)matrixPrint(xtwo, Nobj, Ndim, 6, 4);
            }
            for (int i = 0; i < ncor; i++) {
                dtwo[i] = xtwo[i] - 2.0 * xone[i] + xold[i];
            }
            if (DEBUG) {
                printf("dtwo\n");
                (void)matrixPrint(dtwo, Nobj, Ndim, 6, 4);
            }
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
                if (DEBUG) {
                    printf("xnew\n");
                    (void)matrixPrint(xnew, Nobj, Ndim, 6, 4);
                }
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
    if (Safe) {
        (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xnew, edis);
        snew = smacofSSLoss(ndat, edis, dhat, wght);
        if (sone < snew) {
            for (int i = 0; i < ncor; i++) {
                xnew[i] = xone[i];
            }
        }
    }
    xfree(xone);
    xfree(xtwo);
    xfree(xthr);
    xfree(done);
    xfree(dtwo);
    xfree(dthr);
    xfree(edis);
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
    double* bmat = xmalloc(Nobj * Nobj * sizeof(double));
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        bmat[i + Nobj * j] = bmat[j + Nobj * i] = wght[k] * dhat[k] / edis[k];
    }
    for (int s = 0; s < Ndim; s++) {
        int sobj = s * Nobj;
        for (int i = 0; i < Nobj; i++) {
            for (int j = 0; j < Nobj; j++) {
                xtmp[i + sobj] +=
                    bmat[i + Nobj * j] * (xold[i + sobj] - xold[j + sobj]);
            }
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
        for (int k = 0; k < Ndat; k++) {
            int i = iind[k];
            int j = jind[k];
            bmat[i + Nobj * j] = bmat[j + Nobj * i] = vinv[k];
        }
        for (int s = 0; s < Ndim; s++) {
            int sobj = s * Nobj;
            for (int i = 0; i < Nobj; i++) {
                for (int j = 0; j < Nobj; j++) {
                    xnew[i + sobj] +=
                        bmat[i + Nobj * j] * (xtmp[i + sobj] - xtmp[j + sobj]);
                }
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
    xfree(bmat);
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