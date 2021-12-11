## à faire

* Celine: il faut programmer également les variables de controle (antithétiques) et l'intervalle de confiance
* vérifier la simulation de S pas à pas (comparer la variance simulée avec la calculée)
* ecdf au liru de histogramme?
* traduire en VBA
  + écrire un test pour VBA, Excel et son graphe ne supporte plus que 250 series
* creer un dashboard
* sauvegarder des capture d'écran du dashboard
* écrire la documentation 20 pages
  + méthode de reduction de variance
  + décrire les résultats et les graphiques (Nd=8 ~ 40j, t=[0,1], nt = 1000, S0 = 40, r=.05, S1=42, K=41)
  + décrire la fonctionalité des programmes (integrale methode des trapezes, simulation pas a pas, fonctionalité de X_N)
  + montrer que les estimateurs sont non-biaisés
  + détailler les différences entre X_inf et X_N
  
## réunion M. Deutsch 09 Nov

Je crois que j'ai rempli tous les points de la réunion, le 24 Nov [Matthias)

* ne utiliser pas S, au lieu de cela des vecteurs par trajectoire
* calculer aussi la variance de l'estimateur
* n bas, mais nt haut
* afficher tout pour le debogage
* calcul de l'integral: ne pas utiliser une double boucle, 0.5\*S_0 + S_1 +...+ S_n + 0.5\*S_{n+1}
* sigma est calculé par S_t -> ajuster sigma (desormais bas)


## QUESTIONS : 
est ce que la méthode choisie est valide ? 
quelle sens pour la matrice ? 
est ce que les calculs sont bons ? 
quel type de dashboard attendez vu sur excel ? 
mieux utiliser X au lieu de C pour intervalle et variable de controle?
Y ne dèpend pas de lambda?
variable de contrôle pour C ou pour S?
Calculation en double boucle dure 14/0.008 sec = 1700 fois (pour n=100, nt=1000) plus longtemps qu'en matrice <-> 3600/7.7 KiB RAM 

n= 1.000;
nt = 100.000;
en matrice => 2,4sec; 702.600 KiB (Mais arrays en plus les limites RAM poseront de problèmes)
en boucle => durée estimée de la calculation: 44mn

Puisque nous n'avons pas trouvé un autre solution de la simulation que la simulation pas a pas, il faut décider entre RAM ou temps = > RAM est moins cher pour une exactitude suffisante.
Le code en boucle est eǵalement sauvegarder, voir appendix ref{X}. - S_simule_boucle.m
Discussion similaire pour la calculation de dWt, autre décision:

% max/min tictoc dWt simulé pas a pas (3Mo)
% Elapsed time is 0.081486 seconds.
% S_simule_matrice
% Elapsed time is 0.037845 seconds.

% max/min tictoc dWt simulé comme matrice (7ko)
% Elapsed time is 0.074900 seconds.
% S_simule_matrice
% Elapsed time is 0.037370 seconds. (indexer plus vite que simuler chaque fois, mais factor ~1.1)

Auch erwähnen das Plot S nur für 15 geplottet wird

Anhand des Boxplots sieht man dass mit dem gewählten K Median Gewinn möglich ist (Mean sogar noch höher) -> K zu hoch als erwartungswert S zu T

IC_gauss wurde gemacht da TCL (siehe auch Bild 2), bootstrap bei C weil keine normal verteilung des Schätzer mehr angenommen werden kann.


% l'erreur relatif à S_t dépend que de dt, alors avec dt desormais
% petit pas de souci (d'ailleurs: interpolation lineaire, simuler,
% choisir N n/a a partie de G, prendre le prix dernier)


Vorteil variable antithetique: 2 Werte pro boucle pas a pas


figure 6 et 6.5: 
% Avantage: IC tres etroite
        % Probleme: K est loin hors de ic, mais EY est dedans ?
        % Peut-etre pcq outliers a cause de la variance(sqrt(S)) ?
        input('\n');
        % explication: pour X grand Y est plus petit
E!=E <> utiliser EX au lieu de EY, variance ist in eine Richtung größer!
Sie auch Grafik 6 wo Z am Ende einen Schwanz hat: Wenn sich X positiv von der Mitte entfernt, dann mit größerer Varianz als X_a nach unten!
Nach der gleichen Rechnung ist K kein fairer Preis
wenn man dS_anti l. 74 ändert, sodass X_a gleiche Varianz wie X hat, hat Z eine normalere Form. Allerdings ist K immer noch verzerrt.
