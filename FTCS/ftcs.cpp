
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

	vector<vector<double> > u_0;
  vector<vector<double> > v_0;
  vector<vector<double> > T;

	/*
	Geschwindigkeitsfelder initialisieren
	*/
	init_u(u_0, v_0);

	/*
	Erzeuge Q(x,y) und T*(x,y)
	*/
	vector<vector<double> > Q;
	vector<vector<double> > T_star;
	if(b_Q);
	init_Q_T_star(Q,T_star);

	/*
	Definiere die Anfangsbedingung fuer T.
	*/
	init_T(T);

	/*
	Integration
	*/
	/*
	Schleife ueber die Zeit, je nachdem ob mit Quellterm oder ohne laut .cfg
	*/

	long int t_start;//TIME LOG
  time(&t_start);

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
			cout << "\r t/t_fin :  "<<n*dt<<"/"<<t_fin<<std::flush;
		}
		cout<<endl;
	}
	else{
		ofstream outfile;
		outfile.open(outname.c_str());
		for(int n=0; n < t_fin/dt; n++){
			outfile << n*dt << " " << SSE(T,T_star)<<endl;
			FTCS_with_Q(T,u_0,v_0, Q);
			// if(i_t<i_t_max-1){
			// 	//Snapshots machen
			// 	if( (n+1)*dt >= t_snap[i_t] && (n+1)*dt<t_snap[i_t+1]){
			// 		ostringstream snap_name;
			// 		snap_name <<dirname<< (n+1)*dt << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
			// 		save_data(T, snap_name.str().c_str());
			// 	}
			// 	if((n+1)*dt>=t_snap[i_t]){
			// 		i_t+=1;
			// 	}
			// }

		}
		outfile.close();
	}
	//TIME LOG
	long int t_finished;
	time(&t_finished);
	cout << "\n\n"<< t_finished-t_start<<endl;

	ostringstream snap_name;
	snap_name <<dirname<< t_fin << "_" << Pe << "_"<< Nx<<"_"<<Ny<<"_"<<dt<<"_"<<b_Q<<".txt";
	save_data(T, snap_name.str().c_str());
	return 0;
}
