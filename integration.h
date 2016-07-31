//integration.h
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <iostream>
#include <cmath>
#include <vector>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <cstdlib>
#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
using namespace std;

/*
Ein Integrationsschritt mit dem FTCS-Verfahren. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
Zielarray, welches ihr uebergeben wird.
*/
void FTCS(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0);
/*
Funktion berechnet SSE(Sum of squared errors) eines Vektors T von einem Vektor T_star.
*/
double SSE (vector<vector< double> > T, vector<vector< double> > T_star);
/*
Funktion um Qij zu berechnen, nach der bestimmten Formel
*/
double Qij(double x, double y);
/*
Funktion um T_star aus Aufgabe 3 zu berechnen
*/
double T_star_ij(double x, double y);
/*
Angepasster Zeitschritt mit Q zur Fehlerbestimmung. Unterschied zu FTCS() ist die Addition von dt*Q.
*/
void FTCS_with_Q(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0, vector<vector<double> > Q);

/*
Definitionen der Funktionen
*/

/*
Variablen Initialisierungen
*/
/*
Geschwindigkeitsfelder initialisieren
*/
void init_u(){
for(int i = 0; i< Nx+1; i++){
		vector<double> u_F;
		vector<double> v_F;
		for(int j = 0; j<Ny+1;j++){
			u_F.push_back(M_PI*sin(2.*M_PI*i*dx)*cos(M_PI*j*dy));
			v_F.push_back(-2.*M_PI*cos(2.*M_PI*i*dx)*sin(M_PI*j*dy));
		}
		u_0.push_back(u_F);
		v_0.push_back(v_F);
	}
}

void init_T(){
	for(int j=0;j<Nx+1;j++){
		std::vector<double> v;
		for(int i = 0; i<Ny+1;i++){
			v.push_back(double(i*1./10.));
		}
		T.push_back(v);
	}
}

void init_Q_T_star(){
	if(b_Q){
		for(int i = 0; i<Nx+1; i++){
			vector<double> F;
			for(int j=0; j<Ny+1;j++){
				F.push_back(Qij(double(i*dx), double(j*dy)));
			}
			Q.push_back(F);
		}
		for(int i = 0; i<Nx+1; i++){
			vector<double> F;
			for(int j=0; j<Ny+1;j++){
				F.push_back(T_star_ij(double(i*dx), double(j*dy)));
			}
			T_star.push_back(F);
		}
	}
}


/*
Ein Integrationsschritt mit dem FTCS-Verfahren. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
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
Angepasster FTCS-Zeitschritt mit Q zur Fehlerbestimmung
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
Ein Integrationsschritt mit dem BTCS-Verfahren. Beinhaltet die Schleife ueber x und y. Funktion veraendert das
Zielarray, welches ihr uebergeben wird.
*/
void BTCS(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0){

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
Funktion die aus einem 2D-Feld einen 1D Vektor macht. Berücksichtigt auch
die Dirichlet Randbedingungen.
*/
std::vector<double> reshape_vector(std::vector<std::vector<double> > T){
	std::cout<<T.size()*T.size()<<std::endl;
	std::vector<double> Vec_M;
	for(int i=0; i<T.size();i++){ //Schleife x
    for(int j=0; j<T.size();j++){ //Schleife y
      //Randbed:
      if(j==0){
				double Sj_low = -(dt/dx/dx+dt*Pe/2./dx*v_0[i][j]);
        Vec_M.push_back(T[i][j+1]-Sj_low*T_unten);//j+1
        j++;
      }
      else if(j==T.size()-2){
				double Sj_up = -(dt/dx/dx-dt*Pe/2./dx*v_0[i][j]);
        Vec_M.push_back(T[i][j]-Sj_up*T_oben);
        j++;
      }
      else{
        Vec_M.push_back(T[i][j]);
      }
    }
  }
	std::cout<<Vec_M.size()<<std::endl;
	return Vec_M;
}

vector<vector<double> > BCTS_implicit_Matrix(){
	//Create (L+D)-Matrix
	vector<double> xtemp((T.size()-2)*(T.size()-2), 0.0);
	vector<vector<double> > MLD((T.size()-2)*(T.size()-2),xtemp);
	//M[k][l] k Spalte, l Zeile

	//Schreibe die Komponenten der DifferenzenGl. in die Matrix.
	for(int k=0; k<MLD.size();k++){
		for(int l=0;l<MLD.size();l++){
			int j = l%(T.size()-2)+k%(T.size()-2)+1;
			int i = l/(T.size()-2)+k/(T.size()-2);
			std::cout<<k<<" "<<l<<": "<<i<<" "<<j<<endl;
			if(i<v_0.size()){
				if(-(k-l)==0){
					MLD[k][l]=(1.+4*dt/dx/dx);//S_diag
				}
				else if(-(k-l)==1){
					MLD[k][l]=-(dt/dx/dx+dt*Pe/2./dx *v_0[i][j]);//Sj_low
				}
				else if(-(k-l)==-1){
					MLD[k][l]=-(dt/dx/dx-dt*Pe/2./dx*v_0[i][j]);//Sj_up
				}
				else if(-(k-l)==T.size()-2){
					MLD[k][l]=-(dt/dx/dx+dt*Pe/2./dx*u_0[i][j]);//Si_low
				}
				else if(-(k-l)==-T.size()+2){
					MLD[k][l]=-(dt/dx/dx-dt*Pe/2./dx*u_0[i][j]);//Si_up
				}
			}
		}
	}
	return MLD;
}
vector<vector<double> > GaussSeidel_Matrix(){
	//Create (L+D)-Matrix
	vector<double> xtemp((T.size()-2)*(T.size()-2), 0.0);
	vector<vector<double> > MLD((T.size()-2)*(T.size()-2),xtemp);
	//M[k][l] k Spalte, l Zeile

	//Schreibe die Komponenten der DifferenzenGl. in die Matrix.
	for(int l=0;l<MLD.size();l++){
		for(int k=0; k<=l;k++){
			int j = l%(T.size()-2)+k%(T.size()-2)+1;
			int i = l/(T.size()-2)+k/(T.size()-2);
			std::cout<<k<<" "<<l<<": "<<i<<" "<<j<<endl;
			if(i<v_0.size()){
				if(-(k-l)==0){
					MLD[k][l]=(1.+4*dt/dx/dx);//S_diag
				}
				else if(-(k-l)==1){
					MLD[k][l]=-(dt/dx/dx+dt*Pe/2./dx *v_0[i][j]);//Sj_low
				}
				else if(-(k-l)==-1){
					MLD[k][l]=-(dt/dx/dx-dt*Pe/2./dx*v_0[i][j]);//Sj_up
				}
				else if(-(k-l)==T.size()-2){
					MLD[k][l]=-(dt/dx/dx+dt*Pe/2./dx*u_0[i][j]);//Si_low
				}
				else if(-(k-l)==-T.size()+2){
					MLD[k][l]=-(dt/dx/dx-dt*Pe/2./dx*u_0[i][j]);//Si_up
				}
			}
		}
	}
	return MLD;
}

#endif
