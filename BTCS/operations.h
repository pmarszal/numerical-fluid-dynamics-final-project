#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <iostream>
#include <cmath>
#include <vector>


std::vector<double> matrix_times_vector(std::vector<std::vector<double> > M, std::vector<double> x){
  //M[i][j] i Spalte, j Zeile
  if(M.size()!=x.size()){
    std::cout<<"Error: Matrix size does not fit the vector size!"<<std::endl;
    return x;
  }
  std::vector<double> y;
  for(int j = 0; j<M[0].size();j++){
    double x_j = 0;
    for(int i = 0; i<x.size(); i++){
      x_j += M[i][j]*x[i];
    }
    y.push_back(x_j);
  }
  return y;
}
double scalar_product(std::vector<double> x, std::vector<double> y){
  if(y.size()!=x.size()){
    std::cout<<"Error: Vector size does not fit the vector size!"<<std::endl;
    return 1.;
  }
  double xy = 0;
  for(int i = 0; i<x.size();i++){
    xy+=x[i]*y[i];
  }
  return xy;
}
std::vector<double> scalar_multiplication(std::vector<double> x, double lambda){
  std::vector<double> y;
  for(int i =0; i<x.size();i++){
    y.push_back(lambda*x[i]);
  }
  return y;
}
std::vector<double> subtract_vector(std::vector<double> x, std::vector<double> y){
  if(y.size()!=x.size()){
    std::cout<<"Error: Vector size does not fit the vector size!"<<std::endl;
    return 1.;
  }
  std::vector<double> r;
  for(int i=0;i<x.size(); i++){
    r.push_back(x[i]-y[i]);
  }
  return r;
}
std::vector<std::vector<double> > inverse_matrix(std::vector<std::vector<double> > M){
  for(int j = 0; j<M[0].size();j++){
    for(int i=1;i<=j;i++){
    }
  }

  return M_inv;
}



#endif
