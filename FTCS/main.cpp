#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

int Nx = 30;
int Ny = 30;
double dt = 0.0002;
double t_fin;
double dx = 1./Nx;
double dy = 1./Ny;
double Pe;

/*
Ein Integrationsschritt. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
Zielarray, welches ihr uebergeben wird.
*/
void timestep(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0){
	vector<vector<double> > T_old = T;
	for(int i =1; i< Nx; i++){
		for(int j = 1; j<Ny;j++){
			//Advektionsterm: 2te Ordnung Eulerschritt zentriert
			double Adv = dt*Pe/2.* (u_0[i][j]/dx*(T_old[i+1][j]-T_old[i-1][j]) + v_0[i][j]/dy*(T_old[i][j+1]-T_old[i][j-1]) );
			//Diffusionsterm: 2te Ordnung zentrierte Differenzen
			double Diff = dt/dx/dx*(T_old[i+1][j]-2.*T_old[i][j]+T_old[i-1][j])+dt/dy/dy*(T_old[i][j+1]-2.*T_old[i][j]+T_old[i][j-1]);
			
			T[i][j] = T_old[i][j]-Adv+Diff;
		}
	}
	/* Schleife um den Rand */
	for(int j = 1; j< Ny; j++){
		//Vorwaertsdifferenz zweiter Ordnung
		T[0][j]= 1./3.*(0.+4.* T[1][j]-T[2][j]);
		//Rueckwaertsdifferenz zweiter Ordnung
		T[Nx][j]= 1./3.*(4.*T[Nx-1][j]-T[Nx-2][j]);
	}
	
	
}
/*
Funktion zum Ausgeben des 2D Arrays zum Debuggen.
*/
void print_array(vector<vector< double> > T){
	for(int j = T[0].size()-1; j>=0; j--){
		
		for(int i=0;i<T.size();i++){
			cout << T[i][j] << "  ";
		}
		cout << endl;
	}
}
/*
Funktion zum Speichern des Ergebnisses
*/
void save_data(vector<vector<double> > T, const char* fname){
	ofstream output ;
	output.open(fname);
	for(int j = T[0].size()-1; j>=0; j--){
		
		for(int i=0;i<T.size();i++){
			output << T[i][j] << "  ";
		}
		output << endl;
	}

}


int main(int argc, char ** argv){
	/*
	Parameter von Konsole erfahren
	*/
	sscanf(argv[1], "%lf", &t_fin);
	sscanf(argv[2],"%lf", &Pe);
	/*
	Geschwindigkeitsfelder initialisieren
	*/
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
	/* 
	Integration 
	*/
	/* 
	Schleife ueber die Zeit
	*/
	for(int n=0; n < t_fin/dt; n++){

		timestep(T,u_0,v_0);		
	}

	print_array(T);

	save_data(T,argv[3]);
	return 0;
}
