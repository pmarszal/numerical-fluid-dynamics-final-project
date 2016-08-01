#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"

int main(int argc, char** argv){

  //load_conf(argv[1]);
  int N = 100;
  vector<double> T(N, 0.0);
  vector<double> b(N, 0.0);
  vector<double> r(N, 0.0);
  vector<double> x(N, 0.0);


  vector<vector<double> > M(N,T);
  vector<vector<double> > LD(N,T);

  for(int i = 0;i<T.size(); i++){
    for(int j = 0;j<T.size(); j++){
      M[i][j]= 1./double(N*N)*(0.5*i*j+j);
    }
    M[i][i] = 2.;
  }
  cout << "Dominance\n";
  for(int i = 0;i<T.size(); i++){
    double sum = 0;
    for(int j = 0;j<T.size(); j++){
      sum+=sqrt(M[i][j]*M[i][j]);
    }
    sum -= sqrt(M[i][i]*M[i][i]);
    cout << sum/sqrt(M[i][i]*M[i][i]) << "\t";
  }
  cout << endl;
  LD = triangularize(M);

  for(int i = 0; i<T.size(); i++){
    T[i]= double(std::rand())/RAND_MAX;
  }

  x = matrix_times_vector(LD,T);

  print_vector(T);

  b = matrix_times_vector(M, T);
  //print_vector(b);
  x = b;
  r = subtract_vector(matrix_times_vector(M,x),b);
  double r_betrag = magnitude(r);
  double omega = 1.;

  SOR(x,M,LD,omega);
  std::cout << "New T"<< std::endl;
  print_vector(x);





return 0;

}
/*

cout << "l k :  ";
for(int k = 0; k<M.size(); k++){if(k<10){cout << k<<"  :   ";}
    else{
      cout << k<<"  :  ";
    }}
cout << endl;
for(int l=0; l<M.size();l++){
  cout << l <<" :";
	for(int k=0;k<M.size();k++){
		int j = (l%(u_0.size()-2)+k%(u_0.size()-2))%(u_0.size()-2)+1;
		int i = (l/(u_0.size()-2)+k/(u_0.size()-2))%u_0.size();
		cout <<" : "<<i<<"," <<j<<" ";}
    cout << endl;
  }

*/
