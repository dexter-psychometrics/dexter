#ifndef elsymean_H
#define elsymean_H

void ElSym0(double *b, int *a, int *first, int *last, int *item1, int* item2, int *nI, int *mS, double *g);
void ElSym(double *b, int *a, int *first, int *last, int *item1, int* item2, int *nI, int *mS, double *g);
void meanElSym(double *b, int *a, int *first, int *last, int *item1, int* item2, int *nI, int *mS, double *g, int *idx);
void meanElSym0(double *b, int *a, int *first, int *last, int *item1, int* item2, int *nI, int *mS, double *g, int *idx);
void ittot_mat0(double *b, int *a, double *c, int *first, int *last, int *nI, int *ms, double *pi);
void ittot_mat(double *b, int *a, double *c, int *first, int *last, int *nI, int *ms, double *pi);
void E(double *b, int *a, int *first, int *last, int *scoretab, int *nI, int *mS, double *expect);
void H(double *b, int *a, int* nPar, int *first, int *last, int *scoretab, int *nI, int *mS, double *Hess);
void H0(double *b, int *a, int* nPar, int *first, int *last, int *nscore, int *nI, int *mS, double *Hess);

#endif
