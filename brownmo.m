%% Parametres

r = 0.05;
sigma = 0.02;
T = 100;
n=2^4;
nt = 5; % nombre de trajectoires
N = 5; % pour le cas discret, verifier que N << n
S0 = 40;
t0=0;

%% Simulation

tic

% simulation en forme matrice:

dt = ((T-t0)/n);
t = t0:dt:T;

S_sim = zeros(length(t),nt);
for i = 1:nt
    S_sim(:,i) = brownmo(S0, r, sigma ,t0, T, n); %brownmo est definie en bas
end

duree= toc; 
% Matthias HPPavilion: ~ 0.22, le script R: ~ 0.09 secs; nt=30 n=2^8
fprintf('%d trajectoires plotees\n', nt);
fprintf('Fini en %0.5g\n', duree);


% plot
tiledlayout(3,1)

x_ax = t; % x-axe

axis([0 T 0.8*min(min(S_sim)) 1.5*max(max(S_sim))]) %x-axe limits
%axis([0 T S0-sigma*n S0*(1+r)^(T-t0)+sigma*n]); 

nexttile
hold on

plot(x_ax,S_sim);

% pour comparison, si j'epargne pour le taux r:
%plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
syms func(x)
obligation(x) = S0*(1+r)^(x-t0);
fplot(obligation, [t0 T], "-k"); 

legend("valeur d'une obligation","les prix S_t des actions")
hold off
nexttile
hold on

% calculer K

%K = S0*exp(r*T);
K = int(obligation,t0,T)/(T-t0);
fplot(obligation, [t0 T], "-k"); 

fprintf('Avec S0 de %d, K = int(obligation,t0,T)/(T-t0) = %0.5g \n', S0, K);

% calculer X_t

%fplot(obligation, [t0 T], "-k"); % obligation
plot([t0 T], [K K], "--k"); % y=K

plot(x_ax(1:(n+1)), MoyMob(S_sim, n), ':');
%plot(x_ax,S_sim);

% calculer X_pr


% calculer le prix de l'option C_t = max(0, X_t-S_t)

Xn = MoyMob(S_sim, n);

legend("K","moyenne a temps t")
hold off
nexttile
hold on

plot(x_ax, Xn-S_sim);
plot([t0 T], [0 0], ":k"); % y=zero

legend("","difference X_t - S_t", 'Location','southwest')
hold off
nexttile
hold on

%C = ( (Xn - S_sim) >= 0 ) .* (Xn - S_sim);
K_m = ones(n+1,nt)*K;
%C = (Xn - K_m); %.*( (Xn - K_m) >= 0 );
C = K_m;

% actualiser
for i = 1:(n+1)
    C(i,:) = (Xn(i,:) - K_m(i,:) * exp( -r*((T-t0)/n)*i )); %* ( 1/(1+r) ^ (((T-t0)/n)*i) );
end
% C_0 = C .* ( C >= 0 );

plot(x_ax, C);
plot([t0 T], [0 0], ":k"); % y=zero

legend("C (valeur a temps t0)")
hold off

% Il y a trop de valeurs >0 ? r?
% pour T>10 il y a trop de valeurs = 0?


%% Monte-Carlo

% supprime la valeur initiale
Xn    = Xn(2:end,:);
S_sim = S_sim(2:end,:);
C     = C(2:end,:);

% Calculer C(S_0)
mean(C, 2);


%%%%%%%%%%%%%%%
%% fonctions %%
%%%%%%%%%%%%%%%

function S = brownmo(X0, mu, sigma, t0, T, n) %x0 

  delta = (T-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((T-t0)/n):T;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp( (mu-(sigma^2)/2) * (tseq-t0) + sigma*W );

end


function X_t = MoyMob(M, t_m)
    X_t = cumsum(M,1);
    t_m2= t_m+1;
    X_t = X_t(1:t_m2,:);
    for i = 1:t_m2
        X_t(i,:) = X_t(i,:)/i; %vectoriel?
    end
end








