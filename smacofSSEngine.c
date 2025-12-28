#include "smacofSS.h"

void smacofSSEngine(const int* nobj, const int* ndim, const int* ndat,
                    const int* nord, const int* safe, int* itel, int *kord,
                    const int* ties, const int* itmax, const int* digits,
                    const int* width, const int* verbose, const int* ordinal,
                    const int* weighted, double* sold, double* snew,
                    const double* eps, int* iind, int* jind, int* iord,
                    int* blks, double* wght, double* edis, double* dhat,
                    double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    double smid = 0.0, difx = 0.0;
    (void)smacofMPInverseV(nobj, ndat, iind, jind, wght, vinv);
    while (true) {
        (void)smacofSSMajorize(nobj, ndim, ndat, itel, kord, nord, iind, jind, iord,
                               safe, weighted, wght, vinv, dhat, xold, xnew);
        (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xnew, edis);
        smid = smacofSSLoss(ndat, edis, dhat, wght);
        difx = 0.0;
        for (int k = 0; k < Nobj * Ndim; k++) {
            difx = fmax(difx, fabs(xold[k] - xnew[k]));
        }
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
                printf("itel %4d kord %4d difx %*.*f sold %*.*f smid %*.*f snew %*.*f\n",
                       *itel, *kord, *width, *digits, difx, *width, *digits, *sold,
                       *width, *digits, smid, *width, *digits, *snew);
            } else {
                printf("itel %4d kord %4d difx %*.*f sold %*.*f snew %*.*f\n", *itel, *kord,
                       *width, *digits, difx, *width, *digits, *sold, *width,
                       *digits, *snew);
            }
        }
        if ((*itel == *itmax) || (difx < *eps)) {
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

double smacofSSLoss(const int* ndat, double* edis, double* dhat, double* wght) {
    int Ndat = *ndat;
    double loss = 0.0;
    for (int k = 0; k < Ndat; k++) {
        loss += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    return loss;
}

void smacofSSNormDhat(const int* ndat, double* dhat, double* wght) {
    int Ndat = *ndat;
    double norm = 0.0;
    for (int k = 0; k < Ndat; k++) {
        norm += wght[k] * SQUARE(dhat[k]);
    }
    norm = sqrt(norm);
    for (int k = 0; k < Ndat; k++) {
        dhat[k] /= norm;
    }
    return;
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
