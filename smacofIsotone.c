#include "smacofSS.h"

void primaryApproach(const int* ndat, int* blks, double* x, double* w,
                     double* d, int* iind, int* jind) {
    int Ndat = *ndat;
    for (int i = 0; i < Ndat; i++) {
        int blksize = blks[i];
        if (blksize > 0) {
            double* extracx = xmalloc(blksize * sizeof(double));
            double* extracw = xmalloc(blksize * sizeof(double));
            double* extracd = xmalloc(blksize * sizeof(double));
            int* extraci = xmalloc(blksize * sizeof(int));
            int* extracj = xmalloc(blksize * sizeof(int));
            for (int j = 0; j < blksize; j++) {
                extracx[j] = x[i + j];
                extracw[j] = w[i + j];
                extracd[j] = d[i + j];
                extraci[j] = iind[i + j];
                extracj[j] = jind[i + j];
            }
            (void)mySort(extracx, extracw, extracd, extraci, extracj, &blksize);
            for (int j = 0; j < blksize; j++) {
                x[i + j] = extracx[j];
                w[i + j] = extracw[j];
                d[i + j] = extracd[j];
                iind[i + j] = extraci[j];
                jind[i + j] = extracj[j];
            }
            xfree(extracx);
            xfree(extracw);
            xfree(extracd);
            xfree(extraci);
            xfree(extracj);
        }
    }
    double* ww = xmalloc(Ndat * sizeof(double));
    for (int i = 0; i < Ndat; i++) {
        ww[i] = w[i];
    }
    (void)monotone(ndat, x, ww);
    xfree(ww);
    return;
}

void secondaryApproach(const int* ndat, int* blks, double* x, double* w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double* xsum = xmalloc(nblk * sizeof(double));
    double* wsum = xmalloc(nblk * sizeof(double));
    double* xave = xmalloc(nblk * sizeof(double));
    int* csum = xmalloc(nblk * sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k];
        }
        l += csum[k];
    }
    xfree(xsum);
    xfree(wsum);
    xfree(csum);
    xfree(xave);
    return;
}

void tertiaryApproach(const int* ndat, int* blks, double* x, double* w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double* xsum = xmalloc(nblk * sizeof(double));
    double* wsum = xmalloc(nblk * sizeof(double));
    double* xave = xmalloc(nblk * sizeof(double));
    double* yave = xmalloc(nblk * sizeof(double));
    int* csum = xmalloc(nblk * sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    for (int k = 0; k < nblk; k++) {
        yave[k] = xave[k];
    }
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k] + (x[i - 1] - yave[k]);
        }
        l += csum[k];
    }
    xfree(xsum);
    xfree(wsum);
    xfree(xave);
    xfree(yave);
    xfree(csum);
    return;
}

void tieBlockAverages(const int* ndat, const int* nblk, int* blks, double* x,
                      double* w, double* xsum, double* wsum, int* csum,
                      double* xave) {
    int iblk = 0, Ndat = *ndat, Nblk = *nblk;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int l = k; l < k + blks[k]; l++) {
                sum1 += w[l] * x[l];
                sum2 += w[l];
            }
            xsum[iblk] = sum1;
            wsum[iblk] = sum2;
            csum[iblk] = blks[k];
            iblk++;
        }
    }
    for (int i = 0; i < Nblk; i++) {
        xave[i] = xsum[i] / wsum[i];
    }
}

// Function monotone(),
// performs simple linear ordered monotone regression
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv
// dot nl) This function is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation. This program is distributed in the hope that it
// will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with this function. If not, see
// <https://www.gnu.org/licenses/>.

void monotone(const int* n, double* x, double* w) {
    double* rx = &x[-1];
    double* rw = &w[-1];
    int* idx = xcalloc(*n + 1, sizeof(int));
    idx[0] = 0;
    idx[1] = 1;
    int b = 1;
    double xbm1 = rx[b];
    double wbm1 = rw[b];
    for (int i = 2; i <= *n; i++) {
        b++;
        double xb = rx[i];
        double wb = rw[i];
        if (xbm1 > xb) {
            b--;
            double sb = wbm1 * xbm1 + wb * xb;
            wb += wbm1;
            xb = sb / wb;
            while (i < *n && xb >= rx[i + 1]) {
                i++;
                sb += rw[i] * rx[i];
                wb += rw[i];
                xb = sb / wb;
            }
            while (b > 1 && rx[b - 1] > xb) {
                b--;
                sb += rw[b] * rx[b];
                wb += rw[b];
                xb = sb / wb;
            }
        }
        rx[b] = xbm1 = xb;
        rw[b] = wbm1 = wb;
        idx[b] = i;
    }
    int from = *n;
    for (int k = b; k > 0; k--) {
        const int to = idx[k - 1] + 1;
        const double xk = rx[k];
        for (int i = from; i >= to; i--) rx[i] = xk;
        from = to - 1;
    }
    xfree(idx);
}
