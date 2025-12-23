#include "smacofSS.h"

void smacofSSEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                    int* itmax, int* digits, int* width, int* verbose,
                    int* ordinal, int* weighted, double* sold, double* snew,
                    double* eps, int* iind, int* jind, int* blks, double* wght,
                    double* edis, double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    (void)smacofMPInverseV(nobj, ndat, iind, jind, wght, vinv);
    while (true) {
        (void)smacofSSMajorize(nobj, ndim, ndat, iind, jind, weighted, wght,
                               vinv, edis, dhat, xold, xnew);
        double smid = smacofSSLoss(ndat, edis, dhat, wght);
        if (*ordinal) {
            for (int k = 0; k < Ndat; k++) {
                dhat[k] = edis[k];
            }
            (void)smacofSSMonotone(ndat, ties, iind, jind, blks, edis, dhat,
                                   wght);
            (void)smacofSSNormDhat(ndat, dhat, wght);
            *snew = smacofSSLoss(ndat, edis, dhat, wght);
        } else {
            *snew = smid;
        }
        if (*verbose) {
            if (*ordinal) {
                printf("itel %4d sold %*.*f smid %*.*f snew %*.*f\n", *itel,
                       *width, *digits, *sold, *width, *digits, smid, *width,
                       *digits, *snew);
            } else {
                printf("itel %4d sold %*.*f snew %*.*f\n", *itel, *width,
                       *digits, *sold, *width, *digits, *snew);
            }
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        for (int k = 0; k < Nobj * Ndim; k++) {
            xold[k] = xnew[k];
        }
        *sold = *snew;
        *itel += 1;
    }
    xfree(vinv);
    return;
}

double smacofSSLoss(int* ndat, double* edis, double* dhat, double* wght) {
    int Ndat = *ndat;
    double loss = 0.0;
    for (int k = 0; k < Ndat; k++) {
        loss += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    return loss;
}

void smacofSSNormDhat(int* ndat, double* dhat, double* wght) {
    int Ndat = *ndat;
    double norm = 0.0;
    for (int k = 0; k < Ndat; k++) {
        norm += wght[k] * SQUARE(dhat[k]);
    }
    norm = sqrt(norm);
    for (int k = 0; k < Ndat; k++) {
        dhat[k] /= norm;
    }
}