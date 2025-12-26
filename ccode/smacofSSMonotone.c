#include "smacofSS.h"

void smacofSSMonotone(const int* ndat, const int* ties, int* iind, int* jind,
                      int* blks, double* edis, double* dhat, double* wght) {
    if (*ties == 1) {
        (void)primaryApproach(ndat, blks, dhat, wght, edis, iind, jind);
    }
    if (*ties == 2) {
        (void)secondaryApproach(ndat, blks, dhat, wght);
    }
    if (*ties == 3) {
        (void)tertiaryApproach(ndat, blks, dhat, wght);
    }
}