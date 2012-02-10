#include"psi.h"
#include"theta.h"
#include"matrix.h"
#include"bcg.h"
#include<complex>
using std::complex;

double kappa = 0.5;
//int a[2] = {1,0};
//int b[2] = {0,1};
// rows*columns
int a[2] = {0,1};
int b[2] = {1,0};
matrix r_n[2] = {matrix(2,2),matrix(2,2)};
matrix r_p[2] = {matrix(2,2),matrix(2,2)};

// initalises some matrices
void init()
{
  matrix r(2,2);
  r.set(0,0,1);
  r.set(0,1,0);
  r.set(1,0,0);
  r.set(1,1,1);
  matrix gamma[2] = {matrix(2,2),matrix(2,2)};
  gamma[0].set(0,0,0);
  gamma[0].set(0,1,1);
  gamma[0].set(1,0,1);
  gamma[0].set(1,1,0);
  gamma[1].set(0,0,1);
  gamma[1].set(0,1,0);
  gamma[1].set(1,0,0);
  gamma[1].set(1,1,-1);
  for (int i=0; i<2; ++i)
  {
    r.subtract(gamma[i],r_n[i]);
    r.add(gamma[i],r_p[i]);
  }
}

// e^i theta
complex<double> phase(theta & l, int x, int y, int i)
{
  double angle = l(x,y,i);
  return complex<double>(cos(angle),sin(angle));
}

//updates x,y of p_in using field t stores in p_out
void update(theta & t, psi & p_in, psi & p_out, int x, int y)
{
  matrix accum(2,1);
  matrix mt(2,1);
  matrix pt(2,1);

  for (int i=0;i<2;++i)
  {
    //temp = temp + (r_n[i] * p(x+a[i],y+b[i])) * phase(t,x,y,i)
    //+ (r_p[i] * p(x-a[i],y-b[i])) * conj(phase(t,x,y,i))
    p_in.put(x+a[i],y+b[i],pt);
    r_n[i].mult(pt,mt);
    mt.scale(phase(t,x,y,i));
    accum.add(mt,accum);

    p_in.put(x-a[i],y-b[i],pt);
    r_p[i].mult(pt,mt);
    mt.scale(conj(phase(t,x-a[i],y-b[i],i)));
    accum.add(mt,accum);
  }
  accum.scale(kappa);
  p_in.put(x,y,pt);
  pt.subtract(accum,mt);
  p_out.set(x,y,mt);
}

//updates x,y of p_in using field t stores in p_out (transpose)
void update_t(theta & t, psi & p_in, psi & p_out, int x, int y)
{
  matrix accum(2,1);
  matrix mt(2,1);
  matrix pt(2,1);

  for (int i=0;i<2;++i)
  {
    p_in.put(x+a[i],y+b[i],pt);
    r_p[i].mult(pt,mt);
    mt.scale(phase(t,x,y,i));
    accum.add(mt,accum);

    p_in.put(x+a[i],y+b[i],pt);
    r_n[i].mult(pt,mt);
    mt.scale(conj(phase(t,x,y,i)));
    accum.add(mt,accum);
  }
  accum.scale(kappa);
  p_in.put(x,y,pt);
  pt.subtract(accum,mt);
  p_out.set(x,y,mt);
}

//applies wilson fermion matrix to p_in stores in p_out
void Wmatrix(theta & t, psi & p_in, psi & p_out)
{
  int sizex=p_in.x;
  int sizey=p_in.y;

  for (int i=0; i<sizex;++i)
    for (int j=0; j<sizey;++j)
      update(t,p_in,p_out,i,j);
}

//applies wilson fermion matrix transpose to p_in stores in p_out
void Wmatrix_t(theta & t, psi & p_in, psi & p_out)
{
  int sizex=p_in.x;
  int sizey=p_in.y;

  for (int i=0; i<sizex;++i)
    for (int j=0; j<sizey;++j)
      update_t(t,p_in,p_out,i,j);
}

//dot product of 2 matrices
complex<double> dot(matrix & A, matrix & B)
{
  int rows=A.getrows();
  int cols=A.getcols();
  complex<double> accum=0;
  for (int i=0;i<rows;++i)
    for (int j=0;j<cols;++j)
      accum += A.get(i,j) * B.get(i,j);
  return accum;
}

//dot product of 2 lattices
complex<double> dot(psi & A, psi & B)
{
  int rows=A(0,0).getrows();
  int cols=A(0,0).getcols();
  complex<double> accum=0;
  for (int i=0;i<A.x;++i)
    for (int j=0;j<B.y;++j)
      for (int r=0;r<rows;++r)
        for (int c=0;c<cols;++c)
          accum += A.get(i,j,r,c) * B.get(i,j,r,c);
  return accum;
}

//BCG
//Wx=b
//x_in the initial guess
//x_out the output
//iter number of iterations (no proximity testing)
void bcg(theta & th, psi & b, psi & x_in, psi & x_out, int iter)
{
  int x = b.x;
  int y = b.y;

  complex<double> alpha,beta;
  complex<double> s1, s2;
  psi r(x,y);
  psi r_(x,y);
  psi p(x,y);
  psi p_(x,y);
  psi tmp(x,y);
  x_in.copy(x_out);

  Wmatrix(th,x_in,tmp);
  b.subtract(tmp,r);
  r.copy(r_);
  r.copy(p);
  r.copy(p_);
  //r = b-Wx_0
  //r_ = r
  //p  = r
  //p_ = r_
  for (int i=0; i<iter;++i)
  {
    s1 = dot(r_,r);
    Wmatrix(th,p,tmp);
    s2 = dot(p_,tmp);
    alpha = s1/s2;
    //alpha=(r_*r)/(p_*A*p)

    tmp.scale(alpha);
    r.subtract(tmp,r);
    //r = r-alpha*A*p

    Wmatrix_t(th,p_,tmp);
    tmp.scale(alpha);
    r_.subtract(tmp,r_);
    //r_ = r_ - alpha*A_t*p_

    p.copy(tmp);
    tmp.scale(alpha);
    x_out.add(tmp,x_out);

    //std::cout << d(x_in,x_out) << "\n";
    x_out.copy(x_in);
    x_out.print();
    //psu = psu + p*alpha;

    s2 = dot(r_,r);
    beta = s2/s1;
    //beta=(r_'*r')/(r_*r)

    p.copy(tmp);
    tmp.scale(beta);
    r.add(tmp,p);
    //p=r+p*beta;

    p_.copy(tmp);
    tmp.scale(beta);
    r_.add(tmp,p_);
    //p_=r_+p_*beta;
  }
  //x_out.print();
}

void bcg(matrix & A, matrix & x_in, matrix & b, matrix & x_out, int iter)
{
  int x = b.getrows();
  int y = b.getcols();

  complex<double> alpha,beta;
  complex<double> s1, s2;
  matrix At(A.getcols(),A.getrows());
  matrix r(x,y);
  matrix r_(x,y);
  matrix p(x,y);
  matrix p_(x,y);
  matrix tmp(x,y);
  x_in.copy(x_out);

  b.subtract(tmp,r);
  r.copy(r_);
  r.copy(p);
  r.copy(p_);
  //r = b-Wx_0
  //r_ = r
  //p  = r
  //p_ = r_

  A.trans(At);
  for (int i=0; i<iter;++i)
  {
    s1 = dot(r_,r);
    A.mult(p,tmp);
    s2 = dot(p_,tmp);
    alpha = s1/s2;
    //alpha=(r_*r)/(p_*A*p)

    tmp.scale(alpha);
    r.subtract(tmp,r);
    //r = r-alpha*A*p

    At.mult(p_,tmp);
    tmp.scale(alpha);
    r_.subtract(tmp,r_);
    //r_ = r_ - alpha*A_t*p_

    p.copy(tmp);
    tmp.scale(alpha);
    x_out.add(tmp,x_out);
    x_out.print();
    std::cout << "\n";
    //psu = psu + p*alpha;

    s2 = dot(r_,r);
    beta = s2/s1;
    //beta=(r_'*r')/(r_*r)

    p.copy(tmp);
    tmp.scale(beta);
    r.add(tmp,p);
    //p=r+p*beta;

    p_.copy(tmp);
    tmp.scale(beta);
    r_.add(tmp,p_);
    //p_=r_+p_*beta;
  }
  //x_out.print();
}
