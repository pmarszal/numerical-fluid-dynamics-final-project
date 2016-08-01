#include "/home/marszal/Projects/Num_Str_final/interface.h"
#include "/home/marszal/Projects/Num_Str_final/operations.h"
#include "/home/marszal/Projects/Num_Str_final/integration.h"

int main(int argc, char** argv){

  load_conf(argv[1]);
  string test;
  vector<double> T(10,0);
	double Diff = dt/dx/dx;
	double Adv = dt*Pe/2./dx;
  for(int i=0; i<T.size();i++){ //Schleife x
    for(int j=0; j<T.size();j++){ //Schleife y
      //Randbed:
      if(j==0){
        cout<<i<<","<<j+1<<endl;
        j++;
      }
      else if(j==T.size()-2){
        cout<<i<<","<<j<<endl;
        j++;
      }
      else{
        cout<<i<<","<<j<<endl;
      }
    }
  }


return 0;

}
/*

cout << "l k :  ";
for(int k = 0; k<M.size(); k++){if(k<10){cout << k<<"  :   ";}
    else{
      cout << k<<"  :  ";
    }}
cout << endl;
for(int l=0; l<M.size();l++){
  cout << l <<" :";
	for(int k=0;k<M.size();k++){
		int j = (l%(u_0.size()-2)+k%(u_0.size()-2))%(u_0.size()-2)+1;
		int i = (l/(u_0.size()-2)+k/(u_0.size()-2))%u_0.size();
		cout <<" : "<<i<<"," <<j<<" ";}
    cout << endl;
  }

*/
