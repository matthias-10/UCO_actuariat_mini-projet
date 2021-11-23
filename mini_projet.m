%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% UTF-8                           %
% 30.11.2021                      %
% Valentin DE CRESPIN DE BILLY    %
% Matthias LANG                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% ~~~~~ Mathematiques financieres: Mini-projet 1 ~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


%% ~~~~~~~~~~~~~~~~~~~~ Parametres ~~~~~~~~~~~~~~~~~~~~~ %%


S0 = 40;                % Prix initial du sous jacent
N = 5;                  % Nombre des sous-intervalles 
 % verifier que N << n => a faire: ecrire un test
K = 46;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.04/sqrt(S0);  % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 2^4;               % Nombre de intervalles
T = 3;                  % Fin de la periode
nt = 5;               % Nombre de trajectoires


starttime = datetime('now');
fprintf('La programme a demarre a %s \n', starttime);
fprintf('%d -> Nombre de trajectoires \n', nt);
fprintf('%d -> Prix initial du sous jacent \n', S0)
syms func(x)
obligation(x) = S0*(1+r)^(x-t0);
%K = int(obligation,t0,T)/(T-t0);
bonds_T = obligation(T);
fprintf('%0.5g -> Prix d''une obligation a T\n', bonds_T)
fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf(' . . . ')
tic


%% ~~~~~~~~~~~~~~~~~~~~ Simulation ~~~~~~~~~~~~~~~~~~~~~ %%


dt = (T-t0)/n;
t = t0:dt:T;

S = zeros(n+1,nt);
S(1,:) = S0;

for i = 2:(n+1)
    dW_t = normrnd(zeros(1,nt),sqrt(dt));
    dSi = S(i-1,:).*( r*dt + sigma*sqrt(S(i-1,:)).*dW_t );
    S(i,:) = S(i-1,:) + dSi;
end


%% ~~~~~~~~~~~~~~~~~~ calcul de X_t ~~~~~~~~~~~~~~~~~~~~ %%


X_T = 0.5*S0 + sum(S(2:n,:),1) + 0.5*S(n+1,:);
X_T = X_T/n;

C_inf = X_T-K;
C_inf = C_inf .* ( X_T - K >= 0 );

C_inf_0 = exp(-r*T)*C_inf;

% ~ Estimateur ~
%C_inf * exp(-rT) est une martingale donc 
% E[exp(-rT)*C_inf]= C_inf(S_0)

C_inf_est = mean(C_inf_0);
C_inf_est_var = var(C_inf_0);

fprintf('L''estimateur du C a t0 = %0.5g\n', ...
 C_inf_est);
fprintf('Son écart type = %0.5g\n', sqrt(C_inf_est_var));

%% ~~~~~~~~~~~~~~~~ calcul de X_t_prim ~~~~~~~~~~~~~~~~~ %%

vecX_t_prim = S(1,:);
for i = 2:(n+1)
    vecX_t_prim = vecX_t_prim + S(i,:);
end
vecX_t_prim = vecX_t_prim/(n+1);
X_t_prim = mean(vecX_t_prim);
vecC_N = vecX_t_prim-K;

vecC_N = vecC_N .* ( vecC_N >= 0 );


C_N = mean(vecC_N);
%C_N * exp(-rT) est une martingale donc 
% E[exp(-rT)*C_N]= C_N(S_0)
C_N_0 = exp(-r*T)*C_N;

% Histogramm
% E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum {C(T)}


%% ~~~~~~~~~~~~~~~~~~~~~~~ Plot ~~~~~~~~~~~~~~~~~~~~~~~~ %%


duree= toc;
fprintf('%d trajectoires plotees\n', nt);
fprintf('Avec S0 = %d, K = %0.5g \n', S0, K);
fprintf('L''integrale par (t_0 - T) de X_t, ');
    fprintf('le prix estime C(T) = %0.5g \n', C_inf_0);
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
%plot([C_chapeau C_chapeau], ...
%     [-0.2*nt/sqrt(sigma*nt) nt/sqrt(sigma*nt)], "-k");
%plot([t0 T], [0 0], ":k"); % y=zero


%% ~~~~~~~~~~~~~~~~~~~~ fonctions ~~~~~~~~~~~~~~~~~~~~~~ %%


%% a effacer, aine 

% for i = 1:nt
%     S(:,i) = brownmo(S0, r, sigma ,t0, T, n);
% end

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

% Xn = MoyMob(S, n);

% calculer X_T'
% X_pr = zeros(N,nt);
% for i = 1:N %vectoriel?
%     X_pr(i,:) = S(floor(i*n/N),:);
% end
% X_pr = (1/N)*sum(X_pr, 1);


