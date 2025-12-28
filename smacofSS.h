#ifndef SMACOF_SS_H
#define SMACOF_SS_H

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define true 1
#define false 0

// in smacofIsotone.c

void primaryApproach(const int* ndat, int* blks, double* x, double* w,
                     double* d, int* iind, int* jind);
void secondaryApproach(const int* ndat, int* blks, double* x, double* w);
void tertiaryApproach(const int* ndat, int* blks, double* x, double* w);
void tieBlockAverages(const int* ndat, const int* nblk, int* blks,
                      double* x, double* w, double* xsum,
                      double* wsum, int* csum, double* xave);
void monotone(const int* n, double* x, double* w);
void matrixPrint(const double* x, const size_t Nrow, const size_t Ncol,
                 const int digits, const int width);

// in smacofSort.c

int myComp(const void *, const void *);
void mySort(double* x, double* w, double* d, int* iind, int* jind,
            const int* n);

// in smacofMPInverseV.c

void smacofMPInverseV(const int* nobj, const int* ndat, int* iind, int* jind, double* wght,
                      double* vinv);

// in smacofSSEngine.c

void smacofSSEngine(const int* nobj, const int* ndim, const int* ndat,
                    const int* nord, const int* safe, int* itel, int *kord,
                    const int* ties, const int* itmax, const int* digits,
                    const int* width, const int* verbose, const int* ordinal,
                    const int* weighted, double* sold, double* snew,
                    const double* eps, int* iind, int* jind, int* iord,
                    int* blks, double* wght, double* edis, double* dhat,
                    double* xold, double* xnew);

double smacofSSLoss(const int* ndat, double* edis, double* dhat, double* wght); 

void smacofSSNormDhat(const int* ndat, double* dhat, double* wght);

// in smacofSSMajorize.c

void smacofSSMajorize(const int* nobj, const int* ndim, const int* ndat,
                      const int* itel, int *kord, const int* nord, int* iind, int* jind,
                      const int* iord, const int* safe, const int* weighted,
                      double* wght, double* vinv, double* dhat, double* xold,
                      double* xnew);

void smacofSSGuttmanTransform(const int* nobj, const int* ndim, const int* ndat, 
                              int* iind, int* jind, const int* weighted, 
                              double* wght, double* vinv, double* dhat,
                              double* xold, double* xnew);

void smacofSSDistances(const int* nobj, const int* ndim, const int* ndat, int* iind, int* jind,
                       double* xmat, double* edis);

// in smacofSSMonotone.c

void smacofSSMonotone(const int* ndat, const int* ties, int* iind, int* jind, int* blks,
                      double* edis, double* dhat, double* wght);



static inline void *xmalloc(const size_t size) {
  void *p = malloc(size);
  if (!p && size) {
    fprintf(stderr, "FATAL: malloc(%zu) failed\n", size);
    abort();
  }
  return p;
}

static inline void *xcalloc(const size_t nmemb, const size_t size) {
  if (size && nmemb > SIZE_MAX / size) {
    fprintf(stderr, "FATAL: calloc overflow (%zu,%zu)\n", nmemb, size);
    abort();
  }
  void *p = calloc(nmemb, size);
  if (!p && nmemb && size) {
    fprintf(stderr, "FATAL: calloc(%zu,%zu) failed\n", nmemb, size);
    abort();
  }
  return p;
}

static inline void *xrealloc(void *ptr, const size_t size) {
  void *p = realloc(ptr, size);
  if (!p && size != 0) {
    fprintf(stderr, "FATAL: realloc(%p,%zu) failed\n", ptr, size);
    abort();
  }
  return p;
}

#define xfree(p)                                                               \
  {                                                                            \
    if ((p) != NULL) {                                                         \
      free(p);                                                                 \
      p = NULL;                                                                \
    }                                                                          \
  }

#endif /* SMACOF_SS_H */