#include<complex>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include"matrix.h"
#include"bcg.h"
#include"theta.h"
#include"psi.h"
using std::cout;
using std::endl;
using std::complex;

int size_x=3, size_y;
int iter=size_x;

int main(int argc, char ** argv)
{

  if (argc > 1)
    size_x = strtol(argv[1],NULL,10);
  size_y=size_x;
  if (argc > 2)
    iter = strtol(argv[2],NULL,10);

  srand48(5338);
  init();
  theta th(size_x,size_y);
  psi psi_known(size_x,size_y);
  psi psi_unknown(size_x,size_y);
  psi base(size_x,size_y);
  psi in(size_x,size_y);
  psi out(size_x,size_y);
  cout << std::setprecision(3);

  // check check Wmatrix()
  int x=1;
  for (int i=0;i<size_x;++i)
    for (int j=0;j<size_x;++j)
    {
      base.set(i,j,0,0,x);
      //base.set(i,j,1,0,x+1);
      in.set(i,j,0,0,x+3);
      //in.set(i,j,1,0,x+3);
      ++x;
    }
  base.print();
  
  Wmatrix(th,base,base);
  base.print();
  bcg(th,base,in,out,iter);


  /* basic Wmatrix
  in.set(1,1,0,0,1);
  in.print();
  Wmatrix(th,in,out);
  out.print();

  in.set(1,1,1,0,2);
  in.print();
  Wmatrix(th,in,out);
  out.print();
  */

  /*test - last and first should be the same
  base.fill();
  base.print();
  in.fill();
  Wmatrix(th,base,base);
  base.print();

  bcg(th,base,in,out,iter);
  */

  /* matrix version
  matrix A(size_x,size_y);
  matrix x(size_x,size_y);
  matrix b(size_x,size_y);
  matrix xout(size_x,size_y);
  A.fill();
  x.fill();
  x.print();
  cout << "\n";
  A.mult(x,b);
  x.fill();
  bcg(A,x,b,xout,iter);
  */
  
  return 0;
}
