#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int Nx = 10;
int Ny = 10;
double dt = 1;
double t_fin = 2;
double dx = 1./Nx;
double dy = 1./Ny;


void timestep(vector<vector<double> > T){
	int i = 0;
	int j = 0;
	int Ni = T.size();
	int Nj = T[0].size();
	
	
	
}

void print_array(vector<vector< double> > T){
	for(int j = T[0].size()-1; j>=0; j--){
		
		for(int i=0;i<T.size();i++){
			cout << T[i][j] << "  ";
		}
		cout << endl;
	}
}

int main(int argc, char ** argv){
	


	vector<vector<double> > u_0;
	vector<vector<double> > v_0;
	for(int i = 0; i< Nx+1; i++){
		vector<double> u_F;
		vector<double> v_F;
		for(int j = 0; j<Ny+1;j++){
			u_F.push_back(M_PI*sin(2.*M_PI*i*dx)*cos(M_PI*j*dy));
			v_F.push_back(-2*M_PI*cos(2*M_PI*i*dx)*sin(M_PI*j*dy));
		}
		u_0.push_back(u_F);
		v_0.push_back(v_F);
	}


	/*
	Definiere die Anfangsbedingung fuer T. 
	*/
	vector<vector<double> > T;
	for(int i = 0; i< Nx+1; i++){
		vector<double> F;
		for(int j = 0; j<Ny+1;j++){
			F.push_back(j*dy);
		}
		T.push_back(F);
	}
	print_array(T);

	/* 
	Integration 
	*/
	/* 
	Schleife ueber die Zeit
	*/
	for(int n=0; n < t_fin/dt; n++){
		/*
		Kopiere das Feld zum aktuellen Zeitpunkt, fuer die 
		Berechnung des Feldes zum naechsten Zeitpunktes.
		*/
		vector<vector<double> > T_old;
		for(int i = 0; i<T.size(); i++){
			vector<double> F;
			for(int j= 0; j<T[i].size();j++){
				F.push_back(T[i][j]);
			}  
			T_old.push_back(F);
		}
		//DEBUG
		/*
		print_array(T);
		cout << "Old_array:" << endl;
		print_array(T_old);
		*/
		
		
	}



	cout << "It works!" << endl;
	return 0;
}
