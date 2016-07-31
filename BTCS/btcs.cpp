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
  init_u();
  //print_array(u_0);
  init_T();
  //print_matrix(T);

  //Create a Vector from 2D-Field

  std::vector<double> T_Vec = reshape_vector(T);
  std::vector<double> T_it = T_Vec;

  //low = -1 ; up=+1
//  double Sj_low = 4;// = -(dt/dx/dx+dt*Pe/2./dx*v_0[i][j])
//  double Sj_up = 6;// = -(dt/dx/dx-dt*Pe/2./dx*v_0[i][j])
//  double Si_low = 14;// = -(dt/dx/dx+dt*Pe/2./dx*u_0[i][j])
//  double Si_up = 16;// = -(dt/dx/dx-dt*Pe/2./dx*u_0[i][j])
//  double S_diag = 5.;// = (1.+4*dt/dx/dx)
/*
  std::vector<std::vector<double> > M = BCTS_implicit_Matrix();
  std::vector<std::vector<double> > LD = GaussSeidel_Matrix();
  std::cout<< M.size()<<" "<<M[0].size()<<" : "<<T_Vec.size()<< std::endl;
  for(int n=0; n*dt<t_fin; n++){
    for(int n_gs = 0; n_gs<20; n_gs++){
      T_it = add_vector(T_it, scalar_multiplication(inverse_matrix_multiplication(LD,subtract_vector(matrix_times_vector(M,T_it), T_Vec)),1.));
    }
    T_Vec = T_it;
  }
*/
  return 0;
}
//M[i][j] i Spalte, j Zeile
