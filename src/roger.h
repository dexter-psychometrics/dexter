#ifndef roger_H
#define roger_H

void sampleNRM(double *theta, double *b, int *a, int *i, int *first, int *last, int *m, int *response);
void sampleNRM2(double *theta, double *b, int *a, int *first, int *last, int *nI, int *m, int *score);
//void PV1(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, double *theta);
void PV(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPV, double *theta);
void PVMix(double *b, int *a, int *first, int *last, double *alpha, double *mu, double *sigma, int *score, int *nP, int *nI, int *nPV, double *theta);
void Escore(double *theta, double *score, double *b, int *a, int *first, int *last, int *n);
void PVrecycle(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *scoretb, int *cscoretb, int *nP_, int *nI_, int *maxA_, double *theta);
void recyclePV(double *b, int *a, int *first, int *last, int *nI, int *n, double *mu, double *sigma, double *R);
void recyclePV2(double *b, int *a, int *first, int *last, int *nI, int *n, double *prior, int *nprior, double *R);
void recyclePV_aA(double *b, int *a, int *A, int *first, int *last, int *nI, int *n, double *mu, double *sigma, double *R);
  
#endif
