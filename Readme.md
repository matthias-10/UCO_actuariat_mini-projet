## à faire

* mettre à jour le VBA
    + il y a des problèmes par rapport aux VC: les ICs sont trop grandes quand on joue avec K
    + les IC sont en général trop grandes
    + lambda de VC_1 pour K bas haut
    + pour un K trop élevé c'est clair, pcq Cov(C, VC) est preque 0. Mais pour un K très petit il me semble que la Covariance croisse plus vite.
    + pour sigma ~0.01 tout est bien, pour sigma ~0.03 C_inf -> inf, pour simga >0.05 runtime overflow
* mettre à jour la partie plots de matlab
* le rapport
  + les bases mathématiques (voir le cours)
  + la méthode de trapezes
  + décrire les résultats et les graphiques (Nd=8 ~ 40j, t=[0,1], nt = 1000, S0 = 40, r=.05, S1=42, K=41)
  + décrire la fonctionalité des programmes (integrale methode des trapezes, simulation pas a pas, fonctionalité de X_N)
  + décrire les VC
  + montrer que les estimateurs sont non-biaisés
  + détailler les différences entre X_inf et X_N
  + ajouter des graphiques mis á jour
* sauvegarder des capture de l'écran du dashboard, des graphiques génerés, de VBA, et du code output
* ajouter curseur pour les paramètres du dashboard
* il faut que le tic toc inclus le calcul du lambda
* VC: p-lambda(VC1 - E[VC1])-eta(VC2 - E[VC2]) -- lambda/eta comme coeff. de régresssion
* Data Validation
* Protect, Protect Shee, Allow Edit Range, Formatting
* hide formula bar

* `exp(lambda*M_t-lambda^2/2<M,M>_t)` = martingale => cela sert à calculer l'espérance
* anti-feature VBA: calculer Var/Cov dure très longtemps pour les VC

VC_2: Variance de la VC = T
*** 

Quéstion:

* A-t-on dit que la variable de contrôle doit être indépendant de X ?
