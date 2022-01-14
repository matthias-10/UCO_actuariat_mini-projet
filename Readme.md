## à faire
* mettre à jour la partie plots de matlab
* nettoyer l'affichage de matlab
* en matlab: VC 4 nécessite son espérance !
* le rapport
  + les bases mathématiques (voir le cours)
  + la méthode de trapezes
  + décrire les résultats et les graphiques (Nd=8 ~ 40j, t=[0,1], nt = 1000, S0 = 40, r=.05, S1=42, K=41)
  + décrire la fonctionalité des programmes (integrale methode des trapezes, simulation pas a pas, fonctionalité de X_N)
  + décrire les VC
  + montrer que les estimateurs des VCs sont non-biaisés
  + détailler les différences entre X_inf et X_N
  + ajouter des graphiques mis á jour
  + ajouter l'affichage des scriptes
* sauvegarder des capture de l'écran du dashboard, des graphiques génerés, de VBA, et du code output
* anti-feature VBA: calculer Var/Cov dure très longtemps pour les VC
* feature Excel Dashboard: Sheet and Workbook protected, conditional formating, data validation, curseur

VC_2: Variance de la VC = T

## messages à Valentin 11-12.01.2022

***

Variable de contrôle 5


$\mathrm{e}^{\frac{W_T}{S_0}}$, avec $W_0 = 0$.

\begin{align*}

\mathbb{E}[ \mathrm{e}^{\frac{W_T}{S_0}} ] = 
\mathrm{e}^{\mathbb{E}[\frac{W_T}{S_0}}] + 0.5 Var(\frac{W_T}{S_0}})} = 
\mathrm{e}^{0 +  \frac{1}{2S_0^2}Var(W_T)} =
\mathrm{e}^{\frac{T}{2S_0^2}} =

\end{align*}

***

VC 6 on peut l'effacer.
VC 3 est à toi à prouver l'ésperance, mais c'est facile.

* pour X_prime il y avait une erreur, le boucle ligne 110: for i = 2:(n+1)
* quelque part c'est écrit "simules" pourtant il fallait "simulees"
* les explications pourquoi X_N est plus grand que X_inf sont faites avec la courbe continu et les bornes pour [1*n/N, 2*n/N, ...], qui dépassent toujours la courbe si bien elle monte.
* est-ce que tu peux vérifier et implémenter les résulats suivants? Merci

*

Variable de contrôle 4


On se met dans le modèle Black-Scholes, et on a choisi comme Variable de Contrôle 4 le payoff d'une option d'achat avec la diffusion $dS_t = S_t(rdt + \sigma sqrt{S_0} dW_t$, qui est très pareil à la diffusion donnée dans le sujet.
On peut utiliser le modèle Black-Scholes avec $\beta_{BS} = r$ et $\sigma_{BS} = \sigma sqrt{S_0}$ (les indices indiquent les variables du modèle comme vu dans le cours).
Il n'est même pas nécessaire d'utiliser le théorème de Girsanov, car on sait déjà que $W_t$ est un mouvement brownien.
Alors $C_0 = \mathbb{E}[ \mathrm{e}^{-rT} C_T ]$.

% Ici on peut ajouter bcp de trucs de III Formule de valorisation en temps continu et suivant

Donc $S_T = S_0 \mathrm{e}^{(r-\frac{\sigma^2 S_0}{2})T + \sigma sqrt{S_0} W_T}$

% Dans le code il fallait sauvegarder dans le boucle aussi W_T comme en VC 2 (ligne 293).
% en code matlab lignes 484 et 529: 
% E_C_VC = S0 * exp((r-0.5*sigma^2*S0)T + sigma*sqrt(S0)*W_T)

Alors E(VC_4) = E(C_0) = E[exp(-rT)*max(S_T - K, 0)] ~ exp(-rT)*max(E[S_T] - K, 0), si S_T assez plus grand que K. Plus éxact: E(C_0) = \int_0^inf exp(-rT)(S_T - K) f(S_T) dS_T; avec S_T qui suit une loi log-normale, mais qui peut être approximé par N( S_T, sqrt(Var(S_T) ).
Donc la distribution est connu, on peut faire en sorte, par exemple: 
E(C_0) = exp(-rT)( max(E[S_T] - K, 0) - \int_{-inf}^K f(S_T) * (S_T - K) )

Make a warning that shows or stops VC 4 if (S_T - K)/Var(S_T) < 2
***

R
```
E_S_T = 50
K = 46.5

layout(c(1,2))

# la loi de distribution de S_T
X = rnorm(100000)*sqrt(5)+E_S_T
hist(X-K,breaks= 100, freq= F, main = "histogramm de S_T - K", xlim = c(-3,10), ylim = c(0,0.6))
abline(v=mean(X-K), col="red")
text( -1, 0.13, col= "red", 
      label= paste( "mean(S_T - K) = ", round(mean(X-K),1)))

# la loi de distribution de C_T
C = ifelse(X > K, X-K, 0)
hist(C, breaks=100, freq= F, main = "histogramm de (S_T - K)^+", xlim = c(-3,10), ylim = c(0,0.6))
abline(v=mean(C), col="red")
text(mean(C)+2, 0.2,col="red", 
     labels = paste("mean(C) = ", round(mean(C),1)) )

```
![image](https://user-images.githubusercontent.com/66843529/149406038-9953c24f-68cb-4ba8-b9c3-536ac8ae5ace.png)

