#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

int Nx = 30;
int Ny = 30;
double dt = 0.0002;
double t_fin;
double dx;
double dy;
double Pe;
double * t_snap;
bool b_Q;


string outname;

/*
Ein Integrationsschritt. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
Zielarray, welches ihr uebergeben wird.
*/
void timestep(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0);
/*
Funktion zum Berechnen des Fehlers
*/
double SSE (vector<vector< double> > T ) ;
/*
Funktion um Qij zu berechnen
*/
double Qij(double x, double y);

/*
Angepasster Zeitschritt mit Q zur Fehlerbestimmung
*/
void timestep_with_Q(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0, vector<vector<double> > Q);
/*
Funktion zum Ausgeben des 2D Arrays zum Debuggen.
*/
void print_array(vector<vector< double> > T);
/*
Funktion zum Speichern des Ergebnisses
*/
void save_data(vector<vector<double> > T, const char* fname);
/*
Nur zur Uebersichtlichkeit
*/
void usage_message(char** params);
/*
Laedt eine configurationsdatei 
*/
void load_conf(const char* cfg_name);
void store_line(string key, string value);


int main(int argc, char ** argv){

	if(argc != 2){
		usage_message(argv);
		return 1;
	}
	/*
	Parameter aus dem .cfg file einlesen
	*/
	load_conf(argv[1]);

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
	Erzeuge Q(x,y)
	*/
	vector<vector<double> > Q;
	if(argv[3]=="y"){
		for(int i = 0; i<Nx+1; i++){
			vector<double> F;
			for(int j=0; j<Ny+1;j++){
				F.push_back(Qij(double(i*dx), double(j*dy)));
			}
			Q.push_back(F);
		}
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
	Schleife ueber die Zeit, je nachdem ob mit Quellterm oder ohne laut .cfg
	*/
	if(!b_Q){
		for(int n=0; n < t_fin/dt; n++){
			timestep(T,u_0,v_0);		
		}
	}
	else{
		for(int n=0; n < t_fin/dt; n++){
			timestep_with_Q(T,u_0,v_0, Q);		
		}
	}


	print_array(T);

	save_data(T,outname.c_str());
	return 0;
}


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
Funktion zum Berechnen des Fehlers
*/
double SSE (vector<vector< double> > T ) {
	return 0.;
}

/*
Funktion um Qij zu berechnen
*/
double Qij(double x, double y){
	double sum1, sum2, sum3, sum4;
	sum1 = -Pe*M_PI*M_PI*sin(2.*M_PI*x)*cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*y) ;
	sum2 = -2.*Pe*M_PI*M_PI*cos(2*M_PI*x)*cos(M_PI*y)*cos(M_PI*x)*sin(M_PI*y);
	sum3 = -2.*Pe*M_PI*cos(2.*M_PI*x)*sin(M_PI*y);
	sum4 = 2.*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y);
	double Qij = sum1+sum2+sum3+sum4;
	return Qij;
}
/*
Angepasster Zeitschritt mit Q zur Fehlerbestimmung
*/
void timestep_with_Q(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0, vector<vector<double> > Q){
	vector<vector<double> > T_old = T;
	for(int i =1; i< Nx; i++){
		for(int j = 1; j<Ny;j++){
			//Advektionsterm: 2te Ordnung Eulerschritt zentriert
			double Adv = dt*Pe/2.* (u_0[i][j]/dx*(T_old[i+1][j]-T_old[i-1][j]) + v_0[i][j]/dy*(T_old[i][j+1]-T_old[i][j-1]) );
			//Diffusionsterm: 2te Ordnung zentrierte Differenzen
			double Diff = dt/dx/dx*(T_old[i+1][j]-2.*T_old[i][j]+T_old[i-1][j])+dt/dy/dy*(T_old[i][j+1]-2.*T_old[i][j]+T_old[i][j-1]);
			
			T[i][j] = T_old[i][j]-Adv+Diff+Q[i][j];
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
	output.close();

}
/*
Nur zur Uebersichtlichkeit
*/
void usage_message(char ** params){
	cout<< "Usage: "<< params[0]<< " [config-file]" << endl;
}
/*
Laedt eine Konfigurationsdatei
*/
void load_conf(const char* cfg_name){
	ifstream is_file;
	is_file.open(cfg_name);

	std::string line;
	while( std::getline(is_file, line) )
	{
		std::istringstream is_line(line);
		std::string key;
		if( std::getline(is_line, key, '=') )
		{
			std::string value;
	    	if( std::getline(is_line, value) ) 
	      		store_line(key, value);
	  		}
		}
}
void store_line(string key, string value){
	if(key==  "Nx"){
		Nx = atoi(value.c_str());
		dx = 1./Nx;
		cout << "Nx=" << Nx<< endl;
	}
	else if(key==  "Ny"){
		Ny = atoi(value.c_str());
		dy=1./Ny;
		cout << "Ny=" << Ny<< endl;
	}
	else if(key==  "dt"){
		dt = atof(value.c_str());
		cout << "dt=" << dt<< endl;
	}
	else if(key==  "t_fin"){
		t_fin = atof(value.c_str());
		cout << "t_fin=" << t_fin<< endl;
	}
	else if(key==  "Pe"){
		Pe = atof(value.c_str());
		cout << "Pe=" << Pe<< endl;
	}
	else if(key==  "b_Q"){
		b_Q = bool(atoi(value.c_str()));
		cout << "b_Q=" << b_Q<< endl;
	}
	else if(key==  "outname"){
		outname = string(value);
		cout << "outname=" << outname<< endl;
	}
	else if(key==  "t_snap"){
		vector<double> elements;
		string element;
		std::istringstream line;
		while(std::getline(line, element, ',')){
			elements.push_back(atof(element.c_str()));
		}
		t_snap = new double[elements.size()];
		for(int i=0;i<elements.size();i++){
			t_snap[i]=elements[i];
		}
	}
}


