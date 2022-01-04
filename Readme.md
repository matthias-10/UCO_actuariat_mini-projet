## à faire

* mettre à jour le VBA
    + copier le VBA_Script.txt dans le .xlsm
    + VC a une variance fausse: remplacer par sigma2/n^2
    + `C_N` est faux, la méthode de trapezes n'applique pas; `C_inf != C_N`
    + VC avec C au lieu de X
* mettre à jour la partie plots de matlab
* le rapport
  + parmi autres:
  + décrire les résultats et les graphiques (Nd=8 ~ 40j, t=[0,1], nt = 1000, S0 = 40, r=.05, S1=42, K=41)
  + décrire la fonctionalité des programmes (integrale methode des trapezes, simulation pas a pas, fonctionalité de X_N)
  + montrer que les estimateurs sont non-biaisés
  + détailler les différences entre X_inf et X_N
* sauvegarder des capture de l'écran du dashboard, des graphiques génerés, de VBA, et du code output
* ajouter curseur pour les paramètres du dashboard
* il faut que le tic toc inclus le calcul du lambda
* VC: p-lambda(VC1 - E[VC1])-eta(VC2 - E[VC2]) -- lambda/eta comme coeff. de régresssion
* `exp(lambda*M_t-lambda^2/2<M,M>_t)` = martingale => cela sert à calculer l'espérance


*** 

Quéstion:

* A-t-on dit que la variable de contrôle doit être indépendant de X ?
