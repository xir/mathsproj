#ifndef matrix_h
#define matrix_h
#include<vector>
#include<complex>
class matrix
{
    int rows, columns;
    std::vector<std::vector<std::complex<double> > > m;
  public:
    matrix(int rows=2, int columns=1);
    void print();
    void fill();
    void fill(std::complex<double>);
    void trans(matrix &);

    std::complex<double> get(int, int);
    void set(int, int, std::complex<double>);

    void scale(std::complex<double>);
    //result stored in last argument
    void copy(matrix &);
    void subtract(matrix &, matrix &);
    void add(matrix &, matrix &);
    void mult(matrix &, matrix &);

    int getrows()
    {return rows; }
    int getcols()
    {return columns;}
};
#endif
