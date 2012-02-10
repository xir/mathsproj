#ifndef theta_h
#define theta_h
class theta
{
  public:
    int x, y;
    double * t;
    theta(int x, int y);
    ~theta();

    double & operator() (int m, int n, int i)
    { return t[(m+x)%x + (((n+y)%y) * x) + (i+2)%2];}
};
#endif
