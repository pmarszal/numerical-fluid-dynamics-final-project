/*
#include <iostream>
#include <cmath>
#include <vector>
*/
#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"
//M[i][j] i Spalte, j Zeile
//T[i][j] i x, j y


int main(int argc, char** argv){

  load_conf(argv[1]);

  vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;

  init_u(u_0, v_0);
  init_T(T);

  //Create a Vector from 2D-Field
  std::vector<double> T_Vec = reshape_vector(T);
  //Calculate the implicit Matrix
  std::vector<std::vector<double> > M = BCTS_implicit_Matrix(u_0,v_0);
  std::vector<std::vector<double> > LD = triangularize(M);


/*
  cout << "Dominance\n";
  for(int i = 0;i<T.size(); i++){
    double sum = 0;
    for(int j = 0;j<T.size(); j++){
      sum+=sqrt(M[i][j]*M[i][j]);
    }
    sum -= sqrt(M[i][i]*M[i][i]);
    cout << sum/sqrt(M[i][i]*M[i][i]) << "\t";
  }
  cout<<endl;
*/
//  print_matrix(M);
//MAtrix Indizes
/*
  cout << "l k :  ";
  for(int k = 0; k<M.size(); k++){if(k<10){cout << k<<"  :   ";}
      else{
        cout << k<<"  :  ";
      }}
  cout << endl;
  for(int l=0; l<M.size();l++){
    if(l<10){cout << l <<"   ";}
    else{cout << l <<"  ";}
  	for(int k=0;k<M.size();k++){
  		int j = (l%(u_0.size()-2)+k%(u_0.size()-2))%(u_0.size()-2)+1;
  		int i = (l/(u_0.size()-2)+k/(u_0.size()-2))%u_0.size();
  		cout <<" : "<<i<<"," <<j<<" ";}
      cout << endl;
    }
*/

  double omega=1.0;
  ofstream outfile;
  outfile.open("iterationszahl.txt");

  int i_t = 0;
  for(int n=0; n*dt<t_fin; n++){
    impose_dirichlet(T_Vec,u_0,v_0);
    SOR(T_Vec,M,LD, omega, 0.0001);
    //Snapshots machen
    if( (n+1)*dt >= t_snap[i_t] && (n+1)*dt<t_snap[i_t+1]){
      ostringstream snap_name;
      snap_name <<dirname<< (n+1)*dt << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
      save_data(shape_back(T_Vec, T_unten,T_oben), snap_name.str().c_str());
    }
    if((n+1)*dt>=t_snap[i_t]){
      i_t+=1;
    }
  }
  T=shape_back(T_Vec, T_unten, T_oben);
  //print_array(T);
  save_data(T,"aktuell.txt");
  return 0;
}
//M[i][j] i Spalte, j Zeile

/*
  while(omega <1.){//Schleife über Omega
    int n = 0;
    cout<< "W: "<<omega<<endl;
    vector<double> r_n=T_Vec;
    std::vector<double> T_it = T_Vec;

    double absolute_r= sqrt(scalar_product(r_n, r_n));
    while(absolute_r>0.1){//Schleife über n
      r_n = subtract_vector(matrix_times_vector(M,T_it), T_Vec);

      T_it = subtract_vector(T_it, scalar_multiplication(
        inverse_matrix_multiplication(
        LD, r_n),omega));
      absolute_r = sqrt(scalar_product(r_n, r_n));
      cout<< "\r"<< absolute_r << std::flush;
      n++;
    }
    cout << endl;
    outfile << omega<< " "<< n<< endl;
    omega += 0.2;
  }
  outfile.close();
*/
