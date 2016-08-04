#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <iomanip>
using namespace std;

bool DEBUG=false;

//Allgemeine Integrationsvariablen
int Nx;//Anzahl der Gitterpunkte
int Ny;
double dt;//Zeitschrittlänge
double t_fin;//Endzeitpunkt
double dx;//Gitterkonstanten
double dy;
double Pe;//Pecletzahl
double T_unten = 0.;//Randbedingungen
double T_oben = 1.;

//BTCS Optionen
double r_end;//Abbruchbedingung
double omega;//Relaxationsparameter

//Integration mit Quellterm?
bool b_Q;

//Ausgabe des Feldes
vector<double> t_snap;//Liste der zu speichernden Zeitpunkte
string dirname;//Pfad in dem die Ausgabe gespeichert werden soll
//Ausgabe des Ergebnisses der Abweichung von T zu T*
string outname;



/*
Funktion zum Ausgeben der Mehrdimensionalen Objekte
*/
void print_array(vector<vector< double> > T);
void print_vector(std::vector<double> x);
void print_matrix(std::vector<std::vector< double> > T);

/*
Funktion zum Speichern von Ergebnissen.
*/
void save_data(vector<vector<double> > T, const char* fname);
/*
Nur zur Uebersichtlichkeit.
*/
void usage_message(char** params);
/*
Laedt eine Konfigurationsdatei und weist jeder globalen Variablen den Wert aus der Datei zu.
*/
void load_conf(const char* cfg_name);//Liest die .cfg Datei
void store_line(string key, string value);//Beinhaltet Ordnet jeder globalen Variablen einen Wert aus der .cfg Datei zu


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
Gibt einen Vektor ins Terminal aus (zum Debuggen).
*/
void print_vector(std::vector<double> x){
  for(int i = 0; i<x.size(); i++){
    std::cout<<x[i]<<" ";
  }
  std::cout<<std::endl;
}
/*
Gibt eine Matrix ins Terminal aus (zum Debuggen).
*/
void print_matrix(std::vector<std::vector< double> > T){
	//std::cout<<std::setprecision(2);
	for(int j = 0; j<T.size(); j++){

		for(int i=0;i<T.size();i++){
			std::cout << T[i][j] << "  ";
		}
		std::cout << std::endl;
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
	while( std::getline(is_file, line) ){
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
		if(DEBUG){cout << "Nx=" << Nx<< endl;}
	}
	else if(key==  "Ny"){
		Ny = atoi(value.c_str());
		dy=1./Ny;
		if(DEBUG){cout << "Ny=" << Ny<< endl;}
	}
	else if(key==  "dt"){
		dt = atof(value.c_str());
		if(DEBUG){cout << "dt=" << dt<< endl;}
	}
	else if(key==  "t_fin"){
		t_fin = atof(value.c_str());
		if(DEBUG){cout << "t_fin=" << t_fin<< endl;}
	}
	else if(key==  "Pe"){
		Pe = atof(value.c_str());
		if(DEBUG){cout << "Pe=" << Pe<< endl;}
	}
	else if(key==  "T_unten"){
		T_unten = atof(value.c_str());
		if(DEBUG){cout << "T_unten=" << T_unten<< endl;}
	}
	else if(key==  "T_oben"){
		T_oben = atof(value.c_str());
		if(DEBUG){cout << "T_oben=" << T_oben<< endl;}
	}
	else if(key==  "r_end"){
		r_end = atof(value.c_str());
		if(DEBUG){cout << "r_end=" << r_end<< endl;}
	}
	else if(key==  "omega"){
		omega = atof(value.c_str());
		if(DEBUG){cout << "omega=" << omega<< endl;}
	}
	else if(key==  "b_Q"){
		b_Q = bool(atoi(value.c_str()));
	}
	else if(key==  "outname"){
		outname = string(value);
		if(DEBUG){cout << "outname=" << outname<< endl;}
	}
	else if(key==  "dirname"){
		dirname = string(value);
		if(DEBUG){cout << "dirname=" << dirname<< endl;}
	}
	else if(key==  "t_snap"){
		vector<double> elements;
		string element;
		std::istringstream line(value);
		while(std::getline(line, element, ',')){
			double snap = atof(element.c_str());
			if(DEBUG){cout << snap << endl;}
			elements.push_back(snap);
		}
		for(int i=0;i<elements.size();i++){
			t_snap.push_back(elements[i]);
		}
	}
}




#endif
