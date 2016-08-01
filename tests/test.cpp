#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"

int main(int argc, char** argv){

  load_conf(argv[1]);

  vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;

  init_u(u_0, v_0);
  init_T(T);
  print_array(T);
  vector<double> Tv = reshape_vector(T, u_0, v_0);
  cout<< "Tv"<<endl;
  print_vector(Tv);
  cout << "M"<<endl;
  std::vector<std::vector<double> > M = BCTS_implicit_Matrix(u_0,v_0);
  print_matrix(M);

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

  cout<< "M*Tv"<<endl;
  print_vector(matrix_times_vector(M,Tv));

  vector<vector<double> > T_new = shape_back(Tv, T_unten, T_oben);
  cout<<"T_new"<<endl;
  print_array(T_new);


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
