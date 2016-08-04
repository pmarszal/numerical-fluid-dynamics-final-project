//integration.h
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
using namespace std;


//Erzeugt das Geschwindigkeitsfeld
void init_u(vector<vector<double> > &u_0, vector<vector<double> > &v_0);

//Erzeugt das Temperaturfeld zum Zeitpunkt 0
void init_T(vector<vector<double> > &T);

//Erzeugt das Quellfeld und die stationäre Lösung
void init_Q_T_star(vector<vector<double> > &Q, vector<vector<double> > &T_star);

//Ein Integrationsschritt mit dem FTCS-Verfahren.
void FTCS(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0);

//Funktion berechnet SSE(Sum of squared errors) eines Vektors T von einem Vektor T_star.
double SSE (vector<vector< double> > T, vector<vector< double> > T_star);

//Berechnet das Quellfeld an einem Ort
double Qij(double x, double y);

//Berechnet die Stationäre Lösung
double T_star_ij(double x, double y);

//Ein Integrationsschritt mit der modifizierten Gleichung
void FTCS_with_Q(vector<vector<double> > &T, vector<vector<double> > u_0, vector<vector<double> > v_0, vector<vector<double> > Q);

//Wandelt das 2D Temperaturfeld in einen Vektor um
vector<double> reshape_vector(vector<vector<double> > T);

//Fügt die Dirichlet-Randbedingungen in den Temperaturvektor ein
void impose_dirichlet(vector<double>  &T, vector<vector<double> > u_0, vector<vector<double> > v_0);

//Wandelt einen Temperaturvektor zurück in ein 2D Array
vector<vector<double> > shape_back(vector<double> T, double T_low, double T_high);

//Erzeugt die BTCS-Matrix
vector<vector<double> > BCTS_implicit_Matrix(vector<vector<double> > u_0, vector<vector<double> > v_0);




/*
Definitionen der Funktionen
*/


void init_u(vector<vector<double> > &u_0, vector<vector<double> > &v_0){
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

void init_T(vector<vector<double> > &T){
	for(int j=0;j<Nx+1;j++){
		vector<double> v;
		for(int i = 0; i<Ny+1;i++){
			v.push_back(double(i*1./Nx));
		}
		T.push_back(v);
	}
}

void init_Q_T_star(vector<vector<double> > &Q, vector<vector<double> > &T_star){
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

/*
Funktion die aus einem 2D-Feld einen 1D Vektor macht. Berücksichtigt auch
die Dirichlet Randbedingungen.
*/
vector<double> reshape_vector(vector<vector<double> > T){
	vector<double> Vec_M;
	for(int i=0; i<T.size();i++){ //Schleife x
    for(int j=0; j<T.size();j++){ //Schleife y
      //Randbed:
      if(j==0){//Überspringe wenn am Dirichlet-Rand
        Vec_M.push_back(T[i][j+1]);
        j++;
      }
      else if(j==T.size()-2){
        Vec_M.push_back(T[i][j]);
        j++;
      }
      else{
        Vec_M.push_back(T[i][j]);
      }
    }
  }
	return Vec_M;
}
void impose_dirichlet(vector<double>  &T, vector<vector<double> > u_0, vector<vector<double> > v_0){
	double Diff = dt/dx/dx;
	double Adv = dt*Pe/2./dx;
	for(int k=0; k<T.size();k++){	//Iteriere über das vektorisierte T-feld

		int j = k%(u_0.size()-2)+1;	//Berechne die Orts-koordinaten
		int i = k/(u_0.size()-2)%u_0.size();

		if(j==1){	//Wenn j==1 dann entspricht T(k) dem Gitterpunkt, der an den oberen
							//Rand gekoppelt ist.
							//
			T[k]=T[k]+(Diff+Adv*v_0[i][v_0.size()-2])*T_unten;
		}
		else if(j==u_0.size()-2){//Wenn j==Ny-2 dann entspricht T(k) dem Gitterpunkt,
														//der an den unteren Rand gekoppelt ist.
														//
			T[k]=T[k]+(+Diff-Adv*v_0[i][1])*T_oben;
		}
	}
}

vector<vector<double> > shape_back(vector<double> T, double T_low, double T_high){
	int T_2d_size = 1+sqrt(1+T.size());
	vector<double> t_tmp(T_2d_size, 0.0);
	vector<vector<double> > T_2d(T_2d_size, t_tmp);
	//Randterme setzen
	for(int i=0; i<T_2d_size; i++){
		T_2d[i][0]=T_low;
		T_2d[i][T_2d_size-1]=T_high;
	}
	//Rest setzen
	for(int l=0; l<T.size(); l++){
		int j = l%(T_2d_size-2)+1;
		int i = l/(T_2d_size-2)%T_2d_size;
		T_2d[i][j]=T[l];
	}
	return T_2d;

}

vector<vector<double> > BCTS_implicit_Matrix(vector<vector<double> > u_0, vector<vector<double> > v_0){
	vector<double> xtemp((u_0.size())*(u_0.size()-2), 0.0);
	vector<vector<double> > MLD((u_0.size())*(u_0.size()-2),xtemp);
	//M[k][l] k Spalte, l Zeile

	//Schreibe die Komponenten der DifferenzenGl. in die Matrix.
	for(int k=0; k<MLD.size();k++){
		for(int l=0;l<MLD.size();l++){
			int j = (k%(u_0.size()-2))+1;
			int i = (k/(u_0.size()-2))%u_0.size();
			double Diff= dt/dx/dx;
			double Adv = dt*Pe/2./dx;
			if((k-l)==0){ // k==l => Diagonale
				MLD[k][l]=(1.+4.*Diff);//S_diag
			}
			else if(l-k==1){ // l-k==1 => eins links von der Diagonalen
				if(j!=u_0.size()-2){	//Bedeutet, dass das Diagonalelement hierzu örtlich am anderen Ende liegt,
															//also nicht gekoppelt ist an diesen Eintrag. => 0.
															//Abfrage betrachtet nur an die Diagonale gekoppelte Gitterpunkte.
					MLD[k][l]=(-Diff-Adv*v_0[i][j+1]);//Geschwindigkeit am Ort des rechten Eintrags
				}
			}
			else if(l-k==-1){ // l-k == -1 => eins rechts von der Diagonalen
				if(j!=1){	//Genauso wie eins drüber.
					MLD[k][l]=(-Diff+Adv*v_0[i][j-1]);//Geschwindigkeit am Ort des linken Eintrags
				}
			}
			else if(l-k==u_0.size()-2){ // l-k==Ny-2 => ein Block links von der Diagonalen
				if(i==u_0.size()-2){	//Das bedeutet, dass das Diagonalelement T(N,j) ist.
															//Die Neumann Randbedingungen liefern in der Gleichung für
															//T(N,j) den Term -2*Diff*T(N-1,j).
					MLD[k][l]=(-2.*Diff);
				}
				else if(i!=u_0.size()-1){	//Element muss an die Diagonale gekoppelt sein
					MLD[k][l]=(-Diff-Adv*u_0[i+1][j]);
				}
				else{MLD[k][l]=0;}
			}
			else if(l-k==-u_0.size()+2){	// l-k==-(Ny-2) => ein Block rechts von der Diagonalen
				if(i==1){	//Das bedeutet, dass das Diagonalelement T(0,j) ist.
									//Die Neumann Randbedingungen liefern in der Gleichung für
									//T(0,j) den Term -2*Diff*T(1,j).
					MLD[k][l]=(-2.*Diff);
				}
				else if(i!= 0){	//Element muss an die Diagonale gekoppelt sein
					MLD[k][l]=(-Diff+Adv*u_0[i-1][j]);
				}
				else{MLD[k][l]=0;}
			}
			else{MLD[k][l]=0;}

		}
	}
	return MLD;
}

#endif
