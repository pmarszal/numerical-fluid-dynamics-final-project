#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>

//Funktion zum Berechnen des Produktes M*x
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

//Skalarprodukt
double scalar_product(std::vector<double> x, std::vector<double> y){
  if(y.size()!=x.size()){
    std::cout<<"Error: Vector size does not fit the vector size!"<<std::endl;
    return 1.;
  }
  double xy = 0.;
  for(int i = 0; i<x.size();i++){
    xy+=x[i]*y[i];
  }
  return xy;
}

//Funktion zum Berechnen des Produktes lambda*x
std::vector<double> scalar_multiplication(std::vector<double> x, double lambda){
  std::vector<double> y;
  for(int i =0; i<x.size();i++){
    y.push_back(lambda*x[i]);
  }
  return y;
}
double magnitude(std::vector<double> x){
  return sqrt(scalar_product(x,x));
}

//Funktion zum Berechnen der Differenz zweier Vektoren x-y.
std::vector<double> subtract_vector(std::vector<double> x, std::vector<double> y){
  if(y.size()!=x.size()){
    std::cout<<"Error: Vector size does not fit the vector size!"<<std::endl;
    return x;
  }
  std::vector<double> r;
  for(int i=0;i<x.size(); i++){
    r.push_back(x[i]-y[i]);
  }
  return r;
}

//Funktion zum Berechnen der Differenz zweier Vektoren x-y.
std::vector<double> add_vector(std::vector<double> x, std::vector<double> y){
  if(y.size()!=x.size()){
    std::cout<<"Error: Vector size does not fit the vector size!"<<std::endl;
    return x;
  }
  std::vector<double> r;
  for(int i=0;i<x.size(); i++){
    r.push_back(x[i]+y[i]);
  }
  return r;
}

//Triangularisiere eine Matrix M
std::vector<std::vector<double> > triangularize(std::vector<std::vector<double> > M){
  vector<std::vector<double> > R=M;
  for(int i =0;i<M.size(); i++){
    for (int j = 0; j < i ; j++) {
      R[i][j]=0.;
    }
  }
  return R;
}


//Funktion zum Berechnen des Produktes M^(-1)*x mittels Vorwärtssubstitution
std::vector<double>  inverse_matrix_multiplication(std::vector<std::vector<double> > M, std::vector<double> x){
  /*
  Zähle die Zeile j hoch, dann die Spalte i bis i=j. Summiere alle Werte M_ij*y_i auf und speicher
  [x_j-SUM(M_ij*y_i)]/M_jj in den Ergebnisvektor.
  */
  if(M.size()!=x.size()||M[0].size()!=M.size()){
    std::cout<<"Error: Matrix size does not fit the vector size!"<<std::endl;
    return x;
  }
  std::vector<double> y;
  for(int j = 0; j<M.size();j++){
    double tmp = x[j];
    if(y.size()>0){
      for(int i = 0; i<j; i++){
        tmp -= M[i][j]*y[i];
      }
    }
    tmp=tmp/M[j][j];
    y.push_back(tmp);
  }
  return y;
}

//Berechne die Lösung eines linearen Gleichungssystems mit Successive Over-Relaxation
int SOR(vector<double>  &T, vector<vector<double> > M, vector<vector<double> > LD, double omega, double end_r){
	std::vector<double> x_old=T;
	std::vector<double> x_n = x_old;
	double r = 1000.;
	int n_gs = 0;

  long int t_start;//TIME LOG
  time(&t_start);
  std::cout<<std::endl;

	while(r>end_r){ //Stoppe wenn ein bestimmter Fehler unterschritten wird
		std::vector<double> r_n = subtract_vector(matrix_times_vector(M,x_old), T);
		r = magnitude(r_n);
		r_n = inverse_matrix_multiplication(LD,r_n);
		r_n = scalar_multiplication(r_n,omega);
		x_n = subtract_vector(x_old, r_n);
		n_gs++;
		x_old = x_n;

    long int t_finished;
    time(&t_finished);
    std::cout<<"\r n: "<< n_gs<<"\t\t, w : "<<omega<<"\t, Time: "<< t_finished-t_start <<" s"<< std::flush;
	}
  std::cout<<std::endl;
	T = x_n;
  return n_gs;
}

//Transponiere eine Matrix
vector<vector<double> > transpose(vector<vector<double> > M){
  vector<vector<double> > T = M;
  for(int k=0;k<M.size();k++){
    for(int l=0;l<M.size();l++){
      T[k][l] = M[l][k];
    }
  }
  return T;
}

//M[i][j] i Spalte, j Zeile
#endif
