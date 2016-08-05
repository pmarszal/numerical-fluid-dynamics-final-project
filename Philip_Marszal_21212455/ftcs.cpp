
#include "integration.h"
#include <ctime>


using namespace std;

int main(int argc, char ** argv){
	if(argc != 2){
		usage_message(argv);
		return 1;
	}
	//Lies die .cfg Datei ein
	load_conf(argv[1]);

	//Initialisiere Felder
	vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;
	init_u(u_0, v_0);
	init_T(T);

	//Initialisiere Q und T_star nur wenn gewollt
	vector<vector<double> > Q;
	vector<vector<double> > T_star;
	if(b_Q){
		init_Q_T_star(Q,T_star);
	}


	long int t_start;//TIME LOG
  time(&t_start);

	int i_t = 0; //Snapshotcounter
	int i_t_max = t_snap.size(); //Snapshotgrenze

	//Integration falls ohne Quellterm
	if(!b_Q){
		for(int n=0; n < t_fin/dt; n++){
			//Führe einen Zeitschritt durch
			FTCS(T,u_0,v_0);

			//Snapshots machen
			if(i_t<i_t_max){
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
	//Integration mit Quellterm
	else{
		ofstream outfile;
		outfile.open(outname.c_str());
		for(int n=0; n < t_fin/dt; n++){
			//Speichere die SSE
			outfile << n*dt << " " << SSE(T,T_star)<<endl;
			//Führe einen Zeitschritt mit Quellterm aus
			FTCS_with_Q(T,u_0,v_0, Q);
		}
		outfile.close();
	}

	//TIME LOG
	long int t_finished;
	time(&t_finished);
	cout << "\n\n"<< t_finished-t_start<<endl;

	return 0;
}
