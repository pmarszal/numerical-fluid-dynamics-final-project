#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"

void print_matrix_template(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0);
void print_matrix_template_new(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0);


int main(int argc, char** argv){

  load_conf(argv[1]);
  vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;

  init_u(u_0, v_0);
  init_T(T);

  std::vector<double> T_Vec = reshape_vector(T, u_0, v_0);
  std::vector<std::vector<double> > M = BCTS_implicit_Matrix(u_0,v_0);
  std::vector<std::vector<double> > LD = triangularize(M);
  std::vector<double> T_Vec2 = matrix_times_vector(M,T_Vec);
  print_vector(T_Vec);
  print_matrix_template(M, u_0);
  cout<<endl;
  print_matrix_template_new(M,u_0);
  SOR(T_Vec2, M,LD,1.0,0.0001);

  print_vector(T_Vec2);
  print_vector(T_Vec);
return 0;

}

void print_matrix_template(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0){
  cout << "l k :  ";
  for(int k = 0; k<M.size(); k++){if(k<10){cout << k<<"   :  ";}
      else{
        cout << k<<"  :  ";
      }}
  cout << endl;
  for(int l=0; l<M.size();l++){
    if(l<10)cout << l <<"  :";
    else cout << l <<" :";
  	for(int k=0;k<M.size();k++){
  		int j = (l%(u_0.size()-2)+k%(u_0.size()-2))%(u_0.size()-2)+1;
  		int i = (l/(u_0.size()-2)+k/(u_0.size()-2))%u_0.size();
      if(l==k){cout <<":|"<<i<<"," <<j<<"| ";}
      else{cout <<": "<<i<<"," <<j<<"  ";}
    }
      cout << endl;
    }
}
void print_matrix_template_new(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0){
  cout << "l k :  ";
  for(int k = 0; k<M.size(); k++){if(k<10){cout << k<<"   :  ";}
      else{
        cout << k<<"  :  ";
      }}
  cout << endl;
  for(int l=0; l<M.size();l++){
    if(l<10)cout << l <<"  :";
    else cout << l <<" :";
  	for(int k=0;k<M.size();k++){
  		int j = (k%(u_0.size()-2))%(u_0.size()-2)+1;
  		int i = (k/(u_0.size()-2))%u_0.size();
      if(l==k){cout <<":|"<<i<<"," <<j<<"| ";}
      else{cout <<": "<<i<<"," <<j<<"  ";}
    }
      cout << endl;
    }
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
