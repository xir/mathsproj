#include"psi.h"
#include<cstdlib>
#include<complex>
psi::psi(int m, int n)
{
  x = m; y = n;
  p = new matrix[x*y];
}

psi::~psi()
{
  delete[] p;
}

// returns the r,c element of the matrix at m,n with periodic boundary conditions
std::complex<double> psi::get(int m, int n, int r, int c)
{
  return p[((m+x)%x) + (((n+y)%y) * x)].get(r,c);
}

// sets A to be m,n element
void psi::put(int m, int n, matrix & A)
{
  int Arows=A.getrows();
  int Acols=A.getcols();
  for (int i=0;i<Arows;++i)
    for (int j=0;j<Acols;++j)
      A.set(i,j,get(m,n,i,j));
}

// set the r,c of the matrix at m,n to z
void psi::set(int m, int n, int r, int c, std::complex<double> z)
{
  p[((m+x)%x) + (((n+y)%y) * x)].set(r,c,z);
}

// sets the m,n element to be A
void psi::set(int m, int n, matrix & A)
{
  int Arows=A.getrows();
  int Acols=A.getcols();
  for (int i=0;i<Arows;++i)
    for (int j=0;j<Acols;++j)
      set(m,n,i,j,A.get(i,j));
}

// copies current into A
void psi::copy(psi & A)
{
  for (int i=0; i<x;++i)
    for (int j=0;j<y;++j)
      p[i+x*j].copy(A.p[i+x*j]);
}

// multiples by z in-place
void psi::scale(std::complex<double> z)
{
  for (int i=0; i<x;++i)
    for (int j=0;j<y;++j)
      p[i+x*j].scale(z);
}

// subtracts A from current and stores it in B
void psi::subtract(psi & A, psi & B)
{
  if (x==A.x && y==A.y && x==B.x && y==B.y)
  {
    for (int i=0; i<x;++i)
      for (int j=0;j<y;++j)
        p[i+x*j].subtract(A.p[i+x*j],B.p[i+x*j]);
  }
  else
  {
    std::cout << "Size mismatch in psi::subtract\n";
    exit(1);
  }
}

// adds A to current and stores it in B
void psi::add(psi & A, psi & B)
{
  if (x==A.x && y==A.y && x==B.x && y==B.y)
  {
    for (int i=0; i<x;++i)
      for (int j=0;j<y;++j)
        p[i+x*j].add(A.p[i+x*j],B.p[i+x*j]);
  }
  else
  {
    std::cout << "Size mismatch in psi::add\n";
    exit(1);
  }
}

// prints
void psi::print()
{
  for (int i=0; i<x;++i)
  {
    for (int j=0;j<y;++j)
      p[i+x*j].print();
    std::cout << "\n";
  }
  std::cout << "\n";
}

// fills with random numbers
void psi::fill()
{
  for (int i=0;i<x;++i)
    for (int j=0;j<y;++j)
      p[i+x*j].fill();
}

// fills with z
void psi::fill(std::complex<double> z)
{
  for (int i=0;i<x;++i)
    for (int j=0;j<y;++j)
      p[i+x*j].fill(z);
}

double d(psi & A, psi & B)
{
  if (A.x==B.x && A.y==B.y)
  {
    double accum=0;
    for (int i=0; i<A.x;++i)
      for (int j=0;j<A.y;++j)
      {
        accum+=abs(A.get(i,j,0,0)-B.get(i,j,0,0));
        accum+=abs(A.get(i,j,1,0)-B.get(i,j,1,0));
      }
    return accum;
  }
  else
  {
    std::cout << "Size mismatch in psi::subtract\n";
    exit(1);
  }
}
