#include <iostream>
#include <cmath>
#include <vector>
#include "operations.h"


int main(int argc, char** argv){
  std::vector<std::vector<double> > M;
  std::vector<double> v;
  v.push_back(0.);
  v.push_back(1.);
  M.push_back(v);
  v[0]=1.;
  v[1]=0.;
  M.push_back(v);

  std::vector<double> x;
  x.push_back(1.);
  x.push_back(2.);
  std::vector<double> b;

  b = matrix_times_vector(M,x);

  std::cout << "["<< x[0]<<","<< x[1]<< "]" <<std::endl;
  std::cout << "["<< b[0]<<","<< b[1]<< "]" <<std::endl;


  return 0;
}
