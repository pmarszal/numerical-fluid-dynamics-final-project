#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"
#include <string>
#include <sstream>

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

void print_matrix_template(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0);
void print_matrix_template_new(std::vector<std::vector<double> > M,std::vector<std::vector<double> > u_0);
void print_matrix_str(std::vector<std::vector< string> > T);
vector<vector<string> > BCTS_implicit_Matrix_str(vector<vector<double> > u_0, vector<vector<double> > v_0);



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
  cout<<endl;
  print_matrix_str(BCTS_implicit_Matrix_str(u_0,v_0));
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

vector<vector<string> > BCTS_implicit_Matrix_str(vector<vector<double> > u_0, vector<vector<double> > v_0){
	//Create (L+D)-Matrix
	vector<string> xtemp((u_0.size())*(u_0.size()-2), "0");
	vector<vector<string> > MLD((u_0.size())*(u_0.size()-2),xtemp);
	//M[k][l] k Spalte, l Zeile

	//Schreibe die Komponenten der DifferenzenGl. in die Matrix.
	for(int k=0; k<MLD.size();k++){
		for(int l=0;l<MLD.size();l++){
			int j = (k%(u_0.size()-2))%(u_0.size()-2)+1;
			int i = (k/(u_0.size()-2))%u_0.size();
			//cout << k<<" "<<l <<" : "<<i<<" " <<j << endl;
      ostringstream stream;
      stream << patch::to_string(i)<<","<<patch::to_string(j);
			double Diff= dt/dx/dx;
			double Adv = dt*Pe/2./dx;
			if((k-l)==0){ // k==l => Diagonale
				MLD[k][l]=stream.str();//S_diag
			}
			else if(l-k==1){ // l-k==1 => eins links von der Diagonalen
				if(j!=u_0.size()-2){//An den Dirichlet Randbedingungen ist die Matrix 0
					MLD[k][l]=stream.str();//Sj_low
				}//Nur weil es nur eins entfernt ist heißt es nicht dass es gekoppelt ist hierdran
			}
			else if(l-k==-1){ // l-k = -1 => eins rechts von der Diagonalen
				if(j!=1){//An den Dirichlet Randbedingungen ist die Matrix 0
					MLD[k][l]=stream.str();//Sj_up
				}
			}
			else if(l-k==u_0.size()-2){ // l-k==Ny-2 => ein Block links von der Diagonalen
				if(i==0){//Neumann Randbed.
					MLD[k][l]=stream.str();
				}
				else if(i!=u_0.size()-1){
					MLD[k][l]=stream.str();//Si_low
				}
			}
			else if(l-k==-u_0.size()+2){// l-k==-(Ny-2) => ein Block rechts von der Diagonalen
				if(i==u_0.size()-1){//Neumann Randbed.
					MLD[k][l]=stream.str();
				}
				else if(i!= 0){
					MLD[k][l]=stream.str();//Si_up
				}
			}

		}
	}
	return MLD;
}

void print_matrix_str(std::vector<std::vector< string> > T){
	//std::cout<<std::setprecision(2);
	for(int j = 0; j<T.size(); j++){

		for(int i=0;i<T.size();i++){
			std::cout << T[i][j] << "  ";
		}
		std::cout << std::endl;
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
