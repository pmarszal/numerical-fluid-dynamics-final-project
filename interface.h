#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
using namespace std;

bool DEBUG=false;
int Nx;
int Ny;
double dt;
double t_fin;
double dx;
double dy;
double Pe;
vector<double> t_snap;
bool b_Q;
string dirname;
string outname;

vector<vector<double> > u_0;
vector<vector<double> > v_0;
vector<vector<double> > T;
vector<vector<double> > Q;
vector<vector<double> > T_star;
double T_unten = 0.;
double T_oben = 1.;

/*
Funktion zum Ausgeben des 2D Arrays zum Debuggen.
*/
void print_array(vector<vector< double> > T);
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
void load_conf(const char* cfg_name);
void store_line(string key, string value);


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
	for(int j = 0; j<T.size(); j++){

		for(int i=0;i<T.size();i++){
			std::cout << T[i][j] << "  ";
		}
		std::cout << std::endl;
	}
}



#endif
