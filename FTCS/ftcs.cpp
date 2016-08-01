/*#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
*/
#include "/home/marszal/Projects/Num_Str_final/integration.h"


using namespace std;

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
	init_u();
	/*
	Erzeuge Q(x,y) und T*(x,y)
	*/
	init_Q_T_star();

	/*
	Definiere die Anfangsbedingung fuer T.
	*/
	init_T();

	/*
	Integration
	*/
	/*
	Schleife ueber die Zeit, je nachdem ob mit Quellterm oder ohne laut .cfg
	*/
	int i_t = 0; //Snapshotcounter
	int i_t_max = t_snap.size(); //Snapshotgrenze
	if(!b_Q){
		for(int n=0; n < t_fin/dt; n++){
			FTCS(T,u_0,v_0);
			if(i_t<i_t_max){
				//Snapshots machen
				if( (n+1)*dt >= t_snap[i_t] && (n+1)*dt<t_snap[i_t+1]){
					ostringstream snap_name;
					snap_name <<dirname<< (n+1)*dt << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
					save_data(T, snap_name.str().c_str());
				}
				if((n+1)*dt>=t_snap[i_t]){
					i_t+=1;
				}
			}
		}
	}
	else{
		ofstream outfile;
		outfile.open(outname.c_str());
		for(int n=0; n < t_fin/dt; n++){
			outfile << n*dt << " " << SSE(T,T_star)<<endl;
			FTCS_with_Q(T,u_0,v_0, Q);
			if(i_t<i_t_max-1){
				//Snapshots machen
				if( (n+1)*dt >= t_snap[i_t] && (n+1)*dt<t_snap[i_t+1]){
					ostringstream snap_name;
					snap_name <<dirname<< (n+1)*dt << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
					save_data(T, snap_name.str().c_str());
				}
				if((n+1)*dt>=t_snap[i_t]){
					i_t+=1;
				}
			}

		}
		outfile.close();
	}
	ostringstream snap_name;
	snap_name <<dirname<< t_fin << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
	save_data(T, snap_name.str().c_str());
	return 0;
}
