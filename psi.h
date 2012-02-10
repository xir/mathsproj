#ifndef psi_h
#define psi_h
#include"matrix.h"
#include<complex>
#include<iostream>
class psi
{
  //size of each matrix object the default constructor of matrix
  //2*1
  public:
    int x, y;
    matrix * p;
    psi(int, int);

    ~psi();

    std::complex<double> get(int m, int n, int r, int c);
    void put(int m, int n, matrix &);
    void set(int m, int n, int r, int c, std::complex<double> z);
    // set p[m][n] values to values of matrix
    void set(int m, int n, matrix &);

    matrix operator()(int m, int n)
    {
      return p[m+n*y];
    }

    void scale(std::complex<double>);
    //result stored in last argument
    void copy(psi &);
    void subtract(psi &, psi &);
    void add(psi &, psi &);

    void print();
    void fill();
    void fill(std::complex<double>);
};
double d(psi &, psi &);
#endif
