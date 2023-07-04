#ifndef DX_SHARED_
#define DX_SHARED_

inline double SQR(const double v){ return v*v; };
inline long double SQR(const long double v){ return v*v; };
inline int SQR(const int v){ return v*v; };
inline double CUB(const double v){ return v*v*v; };
inline long double CUB(const long double v){ return v*v*v; };
inline int CUB(const int v){ return v*v*v; };
inline double QRT(const double v){ return v*v*v*v; };
inline long double QRT(const long double v){ return v*v*v*v; };
inline int QRT(const int v){ return v*v*v*v; };


#endif