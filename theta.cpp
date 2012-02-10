#include"theta.h"
theta::theta(int x, int y)
{
  this->x=x;this->y=y;
  t = new double [2*x*y];
}

theta::~theta()
{
  delete[] t;
}
