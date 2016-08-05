
#include "interface.h"
#include "operations.h"
#include "integration.h"
#include <ctime>
//M[i][j] i Spalte, j Zeile
//T[i][j] i x, j y


int main(int argc, char** argv){

  load_conf(argv[1]);

  //Initialisiere Felder
  vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;
  init_u(u_0, v_0);
  init_T(T);

  //Erstelle den Temperaturvektor
  std::vector<double> T_Vec = reshape_vector(T);
  //Berechne die BTCS-Matrix
  std::vector<std::vector<double> > M = BCTS_implicit_Matrix(u_0,v_0);
  //Berechne die untere Dreiecksmatrix
  std::vector<std::vector<double> > LD = triangularize(M);

  //Integration
  long int t_start;//TIME LOG
  time(&t_start);
  cout << "W = " << omega<< endl;

  int i_t = 0;//Zähler für die Snapshots
  for(int n=0; n*dt<t_fin; n++){
    //Füge Dirichlet-Randbedingungen in den Vektor ein
    impose_dirichlet(T_Vec,u_0,v_0);
    //Löse das Gleichungssystem
    SOR(T_Vec,M,LD, omega, r_end);

    //Snapshots
    if( (n+1)*dt >= t_snap[i_t] && (n+1)*dt<t_snap[i_t+1]){
      ostringstream snap_name;
      snap_name <<dirname<< (n+1)*dt << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
      save_data(shape_back(T_Vec, T_unten,T_oben), snap_name.str().c_str());
    }
    if((n+1)*dt>=t_snap[i_t]){
      i_t+=1;
    }
  }
  cout<<endl;
  T=shape_back(T_Vec, T_unten, T_oben);
  //TIME LOG
  long int t_finished;
  time(&t_finished);
  cout << "\n\n"<< t_finished-t_start<<endl;
  //print_array(T);
  save_data(T,"aktuell.txt");
  return 0;
}
