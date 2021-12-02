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
