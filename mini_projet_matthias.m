%% Parametres
tic
r = 0.05;
sigma = 0.04;
T = 3;
n=2^10;
nt = 100; % nombre de trajectoires
N = 5; % pour le cas discret, verifier que N << n
S0 = 40;
t0=0;

% calculer K

%K = S0*exp(r*T);

syms func(x)
obligation(x) = S0*(1+r)^(x-t0);
%K = int(obligation,t0,T)/(T-t0);

%% Simulation
% simulation en forme matrice:

dt = ((T-t0)/n);
t = t0:dt:T;

%simulation de nt trajectoires
S = zeros(length(t),nt);
for i = 1:nt
    S(:,i) = brownmo(S0, r, sigma ,t0, T, n); %brownmo est definie en bas
end


%% calcul de X_t

vecX_t = zeros(1,nt);
for i = 1:nt
    for j = (1:n)
        vecX_t(i) = vecX_t(i) + (S(j,i)+S(j+1,i))/2;
    end
end

vecX_t = vecX_t/n;
X_t = mean(vecX_t(:,1));
vecC_inf = vecX_t-K;

vecC_inf = vecC_inf .* ( vecC_inf >= 0 );

C_inf = mean(vecC_inf);
%C_inf * exp(-rT) est une martingale donc E[exp(-rT)*C_inf]= C_inf(S_0)
C_inf_0 = exp(-r*T)*C_inf;

%% calcul de X_t_prim

vecX_t_prim = S(1,:);
for i = 2:(n+1)
    vecX_t_prim = vecX_t_prim + S(i,:);
end
vecX_t_prim = vecX_t_prim/(n+1);
X_t_prim = mean(vecX_t_prim);
vecC_N = vecX_t_prim-K;

vecC_N = vecC_N .* ( vecC_N >= 0 );


C_N = mean(vecC_N);
%C_N * exp(-rT) est une martingale donc E[exp(-rT)*C_N]= C_N(S_0)
C_N_0 = exp(-r*T)*C_N;

%% Histogramm
% E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum {C(T)}


%% Plot

duree= toc; 
% Matthias HPPavilion: ~ 0.22, le script R: ~ 0.09 secs; nt=30 n=2^8
fprintf('%d trajectoires plotees\n', nt);
fprintf('Avec S0 de %d, K = %0.5g \n', S0, K);
fprintf('L''integrale par (t_0 - T) de X_t, le prix estime C(T) = %0.5g \n', C_inf_0);
fprintf('La moyenne des X_t: C(T) = %0.5g \n', C_N_0);
fprintf('Fini en %0.5g\n', duree);


tiledlayout(2,1)
nexttile
hold on 

x_ax = t; % x-axe
axis([0 T 0.8*min(min(S)) 1.5*max(max(S))]) %x-axe limits

plot(x_ax, S)

% pour comparison, si j'epargne pour le taux r:
%plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
fplot(obligation, [t0 T], "-k"); 

legend("les prix S_t des actions", "sans risque");

hold off

nexttile

histogram( vecC_inf );
title("Histogramm des C(T) pour X_{infinie}");

%legend("","l'estime pour X_inf","l'estime pour X_N'")
%plot([C_chapeau C_chapeau], [-0.2*nt/sqrt(sigma*nt) nt/sqrt(sigma*nt)], "-k");
%plot([t0 T], [0 0], ":k"); % y=zero

%% fonctions

function S = brownmo(X0, mu, sigma, t0, T, n) %x0 
  delta = (T-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((T-t0)/n):T;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp( (mu-(sigma^2)/2) * (tseq-t0) + sigma*W );
end

