# numerical-fluid-dynamics-final-project

##Über die Arbeit

Im Zuge dieses Projektes wird das Verhalten der BTCS und FTCS am Beispiel der Integration der Konvektions-Diffusions-Gleichung behandelt. Es wurden zwei C++ Anwendungen zu Integration der Gleichung erstellt, die benutzt wurden um Fehler, Stabilität und Laufzeit der beiden Verfahren zu untersuchen. Der Fehler des FTCS Verfahren nimmt mit feinerer Aufteilung des Gitters ab, führt aber bei gleichbleibender Zeitdiskretisierung bei hohen Peclet-Zahlen zu Instabilität. Die numerische Suche nach der Instabilitätsgrenze stößt bei begrenzten Integrationszeiträumen an ihre Grenze.

Für das in der BTCS-Methode verwendete SOR-Verfahren wurden optimale Werte des Relaxationsparameters bestimmt. Für einen Zeitschritt von 100 und einer Peclet-Zahl von 10, wurde dieser bei einem Wert von etwa 1.7 gefunden.

Der Vergleich zwischen BTCS und FTCS ergibt wenig Abweichung in Genauigkeit. Das BTCS-Verfahren ist jedoch mit einem größeren Programmieraufwand verbunden und braucht zur Integration des gleichen Zeitraums, bei kleinen Zeiten, ca. 60 mal so viel Zeit.

##Die Programmstruktur
Für die beiden Integrationsverfahren gibt es zwei separate Programme. Beide Programme benötigen die Header Dateien "interface.h" und "integration.h", zusätzlich baut btcs.cpp auf dem "operations.h" Header auf.

Zur Ausführung der kompilierten Programme müssen diese mit dem Pfad zur Konfigurationsdatei aufgerufen werden. Diese muss folgendermaßen aufgebaut sein:

Nx=30
Ny=30
dt=100.
t_fin=2.
Pe=10.
T_unten=0.
T_oben=1.
r_end=0.0001
omega=1.0
b_Q=0
outname=/home/marszal/Projects/Num_Str_final/BTCS/data/SSE/SSE_over_time.txt
dirname=/home/marszal/Projects/Num_Str_final/BTCS/data/snapshots/
t_snap=0.000002,0.00555733333333,0.0111126666667,0.016668,0.0222233333333,0.0277786666667,0.033334,0.0388893333333,0.0444446666667,0.05

Hinweis: Es dürfen keine Leerzeichen verwendet werden. Nx und Ny sind Integer Zahlen. dt, Pe, T_unten, T_oben, r_end und omega müssen als float eingegeben werden. b_Q kann 1 oder 0 sein. t_snap ist die Liste der Zeitpunkte zu denen das Temperaturfeld gespeichert werden soll. Sofern die Schrittweite es zulässt wird in jedem Intervall zwischen den Werten ein Snapshot gemacht.
dirname gibt den Pfad zu dem Ordner an, in dem die Snapshots gespeichert werden sollen. Dieser muss mit einem / enden.
outname gibt den Pfad zu der Datei an in die der SSE output gespeichert werden soll.
