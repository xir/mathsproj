#include <vector>
#include <cstdlib>
#include <iostream>
#include "matrix.h"
#include <complex>

using std::complex;
using std::vector;
using std::cout;
using std::endl;

matrix::matrix(int rows, int columns)
{
  this->rows = rows; this->columns = columns;
  m = std::vector<std::vector<complex<double> > > (rows, std::vector<complex<double> >(columns));
}

// prints matrix
void matrix::print()
{
  for (int i=0;i<rows;++i)
  {
    for (int j=0;j<columns;++j)
      cout << get(i,j);
    //cout << "\n";
  }
  //cout << "\n";
}

// fills matrix with random numbers
void matrix::fill()
{
  for (int i=0;i<rows;++i)
    for (int j=0;j<columns;++j)
      set(i,j,complex<double>(10*drand48(),10*drand48()));
}

// fills matrix with z
void matrix::fill(complex<double> z)
{
  for (int i=0;i<rows;++i)
    for (int j=0;j<columns;++j)
      set(i,j,z);
}

void matrix::trans(matrix & A)
{
  if (rows == A.getcols() && columns==A.getrows())
  {
    for (int i=0;i<rows;++i)
      for (int j=0;j<rows;++j)
      {
        A.set(j,i,get(i,j));
      }
  }
}

// return i,j element
std::complex<double> matrix::get (int i, int j)
{
  return m[i][j];
}

// sets i,j element to k
void matrix::set(int i, int j, complex<double> k)
{
  m[i][j] = k;
}

// copies current matrix into x
void matrix::copy(matrix & x)
{
  if (rows==x.getrows() && columns==x.getcols())
  {
    for (int i=0;i<rows;++i)
      for (int j=0;j<columns;++j)
        x.set(i,j,get(i,j));
  }
  else
  {
    cout << "Size mismatch in matrix::copy\n";
    exit(1);;
  }
}

// scales matrix by z
void matrix::scale(complex<double> z)
{
  for (int i=0;i<rows;++i)
    for (int j=0;j<columns;++j)
      set(i,j,get(i,j)*z);
}

// subtracts x from matrix and stores it in y
void matrix::subtract(matrix & x, matrix & y)
{
  if (rows==x.rows && columns==x.columns && rows==y.rows && columns==y.columns)
  {
    int yrows=y.getrows();
    int ycols=y.getcols();

    for (int i=0;i<yrows;++i)
      for (int j=0;j<ycols;++j)
      {
        complex<double> t = get(i,j) - x.get(i,j);
        y.set(i,j,t);
      }
  }
  else
  {
    cout << "Size mismatch in matrix::subtract\n";
    exit(1);
  }
}

// adds x to matrix and stores it in y
void matrix::add(matrix &x, matrix & y)
{
  if (rows==x.rows && columns==x.columns && rows==y.rows && columns==y.columns)
  {
    int yrows=y.getrows();
    int ycols=y.getcols();

    for (int i=0;i<yrows;++i)
      for (int j=0;j<ycols;++j)
      {
        complex<double> t = get(i,j) + x.get(i,j);
        y.set(i,j,t);
      }
  }
  else
  {
    cout << "Size mismatch in matrix::add\n";
    exit(1);
  }
}

// multiplies matrix by x and stores it in y
void matrix::mult(matrix & x, matrix & y)
{
  if (rows==y.rows && columns==x.rows && x.columns==y.columns)
  {
    int yrows=y.getrows();
    int ycols=y.getcols();
    matrix temp(yrows,ycols);
    for (int i=0;i<rows;++i)
      for (int j=0;j<x.columns;++j)
      {
        complex<double> sum = 0;
        for (int k=0;k<columns;++k)
          sum += get(i,k) * x.get(k,j);
        temp.set(i,j,sum);
      }
    for (int i=0;i<yrows;++i)
      for (int j=0;j<ycols;++j)
        y.set(i,j,temp.get(i,j));
  }
  else
  {
    cout << "Size mismatch in matrix::mult\n";
    exit(1);
  }
}
