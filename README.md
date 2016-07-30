# numerical-fluid-dynamics-final-project

Aufgabe 1
In einem Rechteck wird (etwa durch einen Ventilator) ein Geschwindigkeitsfeld v
angefacht. Zwei (gegen uberliegende) Seiten des Rechtecks sind w ̈armeisolierend, die
anderen beiden Seiten haben festgelegte Temperaturen. In dieser Aufgabe geht es
darum, die zweidimensionale Temperaturverteilung in dem Rechteck zu bestimmen.
Nach einer geeigneten Entdimensionalisierung lautet das Problem fur das Temperaturfeld T (t, x, y) mit der Peclet-Zahl P e:

Ist v 0 divergenzfrei? Skizzieren Sie v 0 .

Aufgabe 2
Schreiben Sie ein Programm, das diese dimensionslose Gleichung mit dem FTCS-
Schema integriert. Achten Sie darauf, daß die r ̈aumliche Diskretisierung insgesamt
zweiter Ordnung bleibt. Sinnvolle Parameter f ̈
ur Tests sind P e = 2, ∆t = 0.001,
N x = N y = 10. Dabei ist ∆t die Gr ̈oße des Zeitschrittes und das zu berechnende
Gebiet wird auf einem (N x + 1) × (N y + 1) Gitter diskretisiert, wobei die Punkte
auf dem Rand die Indizes 0 und N x bzw. N y haben.

Aufgabe 3
Es stellt sich nun die Frage, ob das Programm auch fehlerfrei ist. Da wir das
urspr ̈
ungliche Problem nicht analytisch l ̈osen k ̈onnen, w ̈ahlen wir ein beliebiges
T ∗ (x, y) das die Randbedingungen erf ̈
ullt, und setzen es in die Gleichung ein. T ∗ ist
nur dann eine L ̈osung, wenn wir noch einen Quellterm Q(x, y) einf ̈
uhren:
Wie muß Q fur gewahlt werden? Fuhren Sie das Q in
Ihr Programm ein und testen Sie, ob dann die Zeitintegration auch wirklich in die
stationare Losung lauft. Den Fehler in der numerischen L ̈osung berechnen wi
als uber alle Gitterpunkte lauft. Wie verhalt, wobei die Summe u
sich der Fehler als Funktion der Aufl ̈osung N x falls N x = N y ?

Aufgabe 4
Ab jetzt ist wieder Q = 0 und es geht um das System (1)-(4). Visualisieren Sie das
Temperaturfeld f ̈
ur P e = 2, N x = N y = 30, ∆t = 2 × 10 −4 zu den Zeiten t = 0.005,
0.05, und 0.5. Wiederholen Sie dasselbe f ̈
ur P e = 10.

Aufgabe 5
Testen Sie mit Ihrem Programm, welches die maximal zul ̈assigen Zeitschritte bei
P e = 0.1, 1 und 10 sind. Was beschrankt die Große des Zeitschrittes?

Aufgabe 6
Nun soll die Gleichung mit einem voll impliziten Euler-Schritt, d.h. mit dem BTCS-
Schema simuliert werden. Dabei f ̈allt die L ̈osung eines linearen Gleichungssystems
 ̈
M x = b an. Dieses System wird mit dem Gauß-Seidel-Verfahren und Uberrelaxation
gel ̈ost. Um das L ̈osungsverfahren zu testen erzeugen Sie erst ein Hilfsprogramm, in
dem zu einem beliebig gew ̈ahlten x das b ausgerechnet wird, und das dann ausgehend
von diesem b das x wieder findet.
Sei M = L + D + U , wobei L und U untere und obere Dreiecksmatritzen und D
eine Diagonalmatrix ist. Die Iterierten innerhalb des Verfahrens seien x n mit den
zugeh ̈origen Residuen r n = M x n − b. Mit dem Relaxationsparameter ω wird das
Verfahren zu:
Bestimmen Sie f ̈
ur P e = 10, N x = N y = 30, ∆t = 100 das optimale ω, indem Sie
|r n | als Funktion von n berechnen und die Anzahl Iterationen bestimmen, die n ̈otig
sind, damit |r n | < 10 −4 . Tragen Sie diese Iterationszahl als Funktion von ω auf. Was
ist das optimale ω?
Verwenden Sie dann dieses L ̈osungsverfahren im BTCS-Schema, um die Ergebnisse
von Aufgabe 4 f ̈
ur P e = 10 zu reproduzieren. Wie verhalten sich die Laufzeiten des
expliziten und impliziten Verfahrens zueinander?
