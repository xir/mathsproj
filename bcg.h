#ifndef bcg_h
#define bcg_h
#include"theta.h"
#include"psi.h"
#include<complex>
void init();
std::complex<double> phase(theta &, int, int, int);
void update(theta &, psi &, psi &, int, int);
void update_t(theta &, psi &, psi &, int, int);
void Wmatrix(theta &, psi &, psi &);
void Wmatrix_t(theta &, psi &, psi &);
std::complex<double> dot(matrix &, matrix &);
std::complex<double> dot(psi &, psi &);
void bcg(theta &, psi &, psi &, psi &, int);
void bcg(matrix &, matrix &, matrix &, matrix &, int);
#endif
