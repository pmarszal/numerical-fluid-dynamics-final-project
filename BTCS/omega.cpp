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


  double omega=1.596;

  ofstream outfile;
  std::vector<double> x(T_Vec.size(),0.0);
  std::vector<double> y(T_Vec.size(),0.0);
  y = T_Vec;
  impose_dirichlet(y, u_0, v_0);
  while (omega<1.599){
    int n_gs=0;
    for(int k = 0; k<T_Vec.size(); k++){x[k]=y[k];}
    outfile.open(outname.c_str(), std::ios_base::app);
    n_gs=SOR(x, M, LD, omega, r_end);
    outfile<<omega<<" "<<n_gs<<endl;
    outfile.close();
    omega+=0.001;
    std::cout<<std::endl;
  }



  return 0;
}
