#include "integration.h"
/*
Ein Integrationsschritt. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
Zielarray, welches ihr uebergeben wird.
*/
void FTCS(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0){
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
		T[0][j]= 1./3.*(4.* T[1][j]-T[2][j]);
		//Rueckwaertsdifferenz zweiter Ordnung
		T[Nx][j]= 1./3.*(4.*T[Nx-1][j]-T[Nx-2][j]);
	}
}
/*
Funktion zum Berechnen des Fehlers
*/
double SSE (vector<vector< double> > T , vector<vector<double> > T_star) {
	double Sum_SE = 0.;
	for(int i = 0; i<T.size();i++){
		for(int j = 0; j<T[i].size();j++){
			Sum_SE=Sum_SE+(T[i][j]-T_star[i][j])*(T[i][j]-T_star[i][j]);
		}
	}
	return sqrt(Sum_SE);
}

/*
Funktion um Qij zu berechnen
*/
double Qij(double x, double y){
	double sum1, sum2, sum3, sum4;
	sum1 = Pe*M_PI*M_PI*sin(2.*M_PI*x)*cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*y) ;
	sum2 = 2.*Pe*M_PI*M_PI*cos(2*M_PI*x)*cos(M_PI*y)*cos(M_PI*x)*sin(M_PI*y);
	sum3 = 2.*Pe*M_PI*cos(2.*M_PI*x)*sin(M_PI*y);
	sum4 = 2.*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y);
	double Qij = sum4-sum2-sum3-sum1;
	return Qij;
}
/*
Funktion um T_star zu berechnen
*/
double T_star_ij(double x, double y){
	return cos(M_PI*x)*sin(M_PI*y)+y;
}
/*
Angepasster Zeitschritt mit Q zur Fehlerbestimmung
*/
void FTCS_with_Q(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0, vector<vector<double> > Q){
	vector<vector<double> > T_old;
	T_old=T;
	for(int i =1; i< Nx; i++){
		for(int j = 1; j<Ny;j++){
			//Advektionsterm: 2te Ordnung Eulerschritt zentriert
			long double Adv = dt*Pe/2.* (u_0[i][j]/dx*(T_old[i+1][j]-T_old[i-1][j]) + v_0[i][j]/dy*(T_old[i][j+1]-T_old[i][j-1]) );
			//Diffusionsterm: 2te Ordnung zentrierte Differenzen
			long double Diff = dt/dx/dx*(T_old[i+1][j]-2.*T_old[i][j]+T_old[i-1][j])+dt/dy/dy*(T_old[i][j+1]-2.*T_old[i][j]+T_old[i][j-1]);
			T[i][j] = T_old[i][j]-Adv+Diff+dt*Q[i][j];
		}

	}
	/* Schleife um den Rand */
	for(int j = 1; j< Ny; j++){
		//Vorwaertsdifferenz zweiter Ordnung
		T[0][j]= 1./3.*(4.* T[1][j]-T[2][j]);
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
