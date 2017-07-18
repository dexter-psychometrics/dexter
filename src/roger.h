#ifndef roger_H
#define roger_H

void sampleNRM(double *theta, double *b, int *a, int *i, int *first, int *last, int *m, int *response);
//void PV1(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, double *theta);
void PV(double *b, int *a, int *first, int *last, double *mu, double *sigma, int *score, int *pop, int *nP,int *nI, int *nPop, int *nPV, double *theta);
void Escore(double *theta, double *score, double *b, int *a, int *first, int *last, int *n);
  
#endif
