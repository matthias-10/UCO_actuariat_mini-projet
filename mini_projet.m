% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Valentin DE CRESPIN DE BILLY                      UTF-8 %
% Matthias LANG                                30.11.2021 %
% exige:                                                  %
% - Statistics and Machine Learning Toolbox               %
% - Symbolic Math Toolbox                                 %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ~~~~~~ Mathematiques financieres: Mini-projet 1 ~~~~~~~ %

%% ~~~~~~~~~~~~~~~~~~~~ Parametres ~~~~~~~~~~~~~~~~~~~~~ %%

S0 = 40;                % Prix initial du sous jacent
%K = 41;                 % Prix d'exercice de l'option
% Le prix sera calculer automatiquement plus tard

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode
Nd = 8;                 % Nombre des sous-intervalles 

nt = 1000;              % Nombre de trajectoires

alpha = 0.05;           % niveau au risque


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

if Nd > n/2-1
warning("Le nombre de sous-intervalles est tres petit")
fprintf('Il fallait Nd << n')
end

starttime = datetime('now');
fprintf('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n');
fprintf('La programme a demarre a %s \n', starttime);
fprintf('%d -> Prix initial du sous jacent \n', S0)

% K basé sur le prix moyen d'une obligation sans risque
syms func(x) 
obligation(x) = S0*(1+r)^(x-t0);
K = double( int(obligation,t0,T)/(T-t0) );
bonds_T = obligation(T);
fprintf('%0.5g -> Prix univers risque neutre a T\n',bonds_T)

fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf('calculation en cours . . .\n')
tic

dt = (T-t0)/n;
t = t0:dt:T;


%% ~~~~~~~~~~~~~~~ simuler pas a pas ~~~~~~~~~~~~~~~~~~~ %%

S = zeros(nt, n+1);

% Methode de Euler

S(:, 1) = S0;
S_anti = S; 
for i = 2:(n+1)
    dWt = normrnd(zeros(nt,1),sqrt(dt));
    dSi = S(:,i-1).* ...
          ( r*dt + sigma*sqrt(S(:,i-1)).*dWt );
    S(:,i)      =    S(:,i-1) + dSi;
    
    % variables antithetiques
    dWt_a = -1*dWt;
    dS_anti = S_anti(:,i-1).* ...
             ( r*dt + sigma*sqrt(S_anti(:,i-1)).*dWt_a );
    S_anti(:,i) = S_anti(:,i-1) + dS_anti;
    %S_anti(:,i) = S_anti(:,i-1) - dSi;
end

%% ~~~~~~~~~~~~~~ C_inf: calcul avec X_T ~~~~~~~~~~~~~~~ %%

% integral: l'aire de t0 a T sous S
X = (0.5*S0 + sum(S(:,2:n),2) + 0.5*S(:,n+1))/n;

C = (X - K) .* logical( X - K >= 0 );
C_0 = exp(-r*T)*C;

% avec variables antithetiques
X_a = (0.5*S0 + sum(S_anti(:,2:n),2) + 0.5*S_anti(:,n+1))/n;

%% ~~~~~~~~~~~~ C_N: calcul avec X_T_prim ~~~~~~~~~~~~~~ %%

%1/N * sum_1^N S_{kT/N}
% => kT n'est pas un numero entier, il faut arrondir

index = fliplr(1:n);
warn_id = 'MATLAB:colon:nonIntegerIndex';
warning('off', warn_id);
% ^supprime Warning a cause de arrondir:
index = index(1:(n/Nd):end); 
X_prim = sum(S(:,index),2)/Nd;

C_prim = (X_prim - K) .* logical( X_prim - K >= 0 );
C_0_prim = exp(-r*T)*C_prim;


%% ~~~~~~~~~~~~~~~~~~~ estimateurs ~~~~~~~~~~~~~~~~~~~~~ %%
    
% ~ Estimateur ~
% C_inf * exp(-rT) est une martingale donc 
% E[exp(-rT)*C_inf]= C_inf(S_0)

% C
C_0_est = mean(C_0);
C_est_var = var(C_0)/nt; %/nt ?

X_mu = mean(X);
C_mu = mean(C);

% C_N
C_0_prim_est = mean(C_0_prim);
C_prim_est_var = var(C_0_prim)/nt; %/nt?

X_prim_mu = mean(X_prim);


%% ~~~~~~~~~~~~~~ intervalle de confiance ~~~~~~~~~~~~~~ %%
%                  (seulement pour N=inf)                 %

v = nt/(nt-1)*var(X); % variance d'echantillonnage

%%% variable supossée normale
X_IC_gauss = [X_mu + sqrt(v/nt)*norminv(alpha/2) ...
              X_mu + sqrt(v/nt)*norminv(1-alpha/2) ];

% variable antithetique
X_a_mu = mean([X;X_a]);
na = 2*nt;
va = na/(na-1)*var([X;X_a]);
X_a_IC_gauss = [X_a_mu + sqrt(va/na)*norminv(alpha/2) ...
                X_a_mu + sqrt(va/na)*norminv(1-alpha/2)];

%%% bootstrap pour C
sims = 10^3;
y = zeros(1, sims);
for i = 1:sims
    y(i) = mean(randsample(C,nt,true)) - C_mu;
end

C_IC_boot =  [C_mu + quantile(y,alpha/2) ...
              C_mu + quantile(y,1-alpha/2) ];


%% ~~~~~~~~~~~~~~ variable de controle ~~~~~~~~~~~~~~~~~ %%
%                (seulement pour X_inf)                   %


% on pourrai utiliser au lieu de la var. antithetique la
% variable de controle suivante

% E(Y) = E(X) [=~ X_a_mu]
Y = X_a;
EY = X_a_mu;

p = corr(X, Y); 
%bien entendu, les deux sont au-peu-pres -1 correles

% optimum: lambda =~ corr(X,Y)*(Var(X)/Var(Y))^.5
lambda = p*(var(X)/var(Y))^.5;
Z = X - lambda * (Y - EY);

Z_a_mu = mean(Z); %([X;Z]);
na = nt;%2*nt;
va = na/(na-1)*var(Z); %([X;Z]);
Z_a_IC_gauss = [Z_a_mu + sqrt(va/na)*norminv(alpha/2) ...
                Z_a_mu + sqrt(va/na)*norminv(1-alpha/2)];

plot(sort(Z))
hold on 
plot(sort(X))
plot([1 na],[K K], '--k', 'LineWidth',2)
hold off
title("X vs variable de controle Z")
legend("Z","X","Z")

% Avantage: IC tres etroite
% Probleme: K est loin hors de ic, mais EY est dedans ?
% Peut-etre pcq outliers a cause de la variance(sqrt(S)) ?

%% ~~~~~~~~~~~~ affichage des estimateurs ~~~~~~~~~~~~~~ %%

duree= toc;
fprintf('\n')
fprintf('%d trajectoires simules\n', nt);
fprintf('Fini en %0.5g secondes\n', duree);
fprintf('\n')

fprintf(' ~ Les estimateurs Monte-Carlo: ~ \n')

fprintf('L''estimateur du C_inf a t0 = %0.5g\n', ...
 C_0_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_est_var));


fprintf(['L''estimateur du C_N a t0, avec ' ...
    '%d sous-intervalles = \n%0.5g\n'], ...
    Nd, C_0_prim_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_prim_est_var));

fprintf('\n ~ Des intervalles de Confiance ~ \n');
fprintf('L''intervalle de confiance de X (normal):\n');
X_IC_gauss
fprintf('La meme intervalle avec var. antithetiques:\n');
X_a_IC_gauss
fprintf('L''intervalle de confiance de Z (normal):\n');
Z_a_IC_gauss
fprintf('L''intervalle de confiance de C (bootstrap):\n');
C_IC_boot

%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%

nt_a = 5; % graphes de S affiches
% 1:   graphes de S; 
% 2-3: ecdf de C_inf et C_N; 
% 4-5: boxplot des estimateurs

G = "g";
P = input(['\n' ...
    'Pour afficher n''importe quel graphique, tapez ' ...
    'son numero <1-5> ou [Enter]. \n' ...
    'Pour quitter tapez plusieures fois [Enter]:\n'] );

if isstring(P) || isempty(P)
    P = 1;
else 
    if ~ismember(P,1:6)
        P = 1;
    end
end

while G~="q"
    disp("[Enter] pour continuer")
    switch P
    case 1
        fprintf('< 1: quelques premiers graphes de S >\n')
        figure(1)
        plot([t0 T],[K K], ':k', 'LineWidth',2)
        hold on
        plot(t, S(1:nt_a,:),'b')
        plot(t, S_anti(1:nt_a,:), 'r:')
        % probleme si nt < nt_a
        plot([t0 T],[K K], '--k', 'LineWidth',2)
        hold off
        % pour comparison, si j'epargne pour le taux r:
        %plot([t0 T], [S0 S0*(1+r)^(T-t0)],"--k"); %obl.
        %1% fplot(obligation, [t0 T], "-k"); 
        legend("K, le prix d''exercice", ...
               "les prix S_t des actions",...
               "les variables antithetiques",...
               "Location","northwest");
        if n*nt > 5000*5000; G="q"; end
        P=P+1; input('\n');
    
    case 2
        if n*nt > 5000*5000; G="q"; end
        fprintf(['< 2: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_infinie >\n'])
        figure(1)
        % E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum{C(T)}
        %histogram( C_inf );
        ecdf( X );
        hold on 
        plot([K K],[0 1], 'k')
        plot([min(X) max(X)], [.5 .5],':b')
        x_ax = [min(X):.1:max(X)]; 
        % probleme si max-min < .1
        nor = normcdf(x_ax,X_mu,sqrt(v));
        plot(x_ax,nor,':r')
        hold off
        legend("ecdf", "K", "P=50%", "cdf normal")
        title("ecdf X(T) pour X_{infinie}");
        P=P+1; input('\n');

    case 3
        if n*nt > 5000*5000; G="q"; end
        fprintf(['< 3: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_N >\n'])
        figure(1)
        ecdf( X_prim );
        hold on 
        plot([K K],[0 1], 'k')
        plot([min(X_prim) max(X_prim)], [.5 .5],':b')
        hold off
        legend("ecdf", "K", "P=50%")
        title("ecdf X(T) pour X_{N}");
        P=P+1; input('\n');

    case 4
        if n*nt > 5000*5000; G="q"; end
        fprintf(['< 4: boxplot de l''estimateur ' ...
                 'C_{infinie} >\n'])
        figure(1)
        boxplot( C_0 );
        title('boxplot de C_{infinie} a T')
        ylabel('C_T, valeurs actualisees')
        P=P+1; input('\n');

    case 5
        if n*nt > 5000*5000; G="q"; end
        fprintf('< 5: boxplot de l''estimateur C_{N} >\n\n')
        figure(1)
        boxplot ( C_0_prim );
        title('boxplot de C_{N} a T')
        P=P+1;

    case 6
        if n*nt > 5000*5000; G="q"; end
        P=input(['\n ' ...
            'Pour afficher n''importe quel graphique, ' ...
            'tapez son numero <1-5> \n']);
        if ismember(P, 1:5)
            fprintf("Vous avez choisi: ")
        else
            G="q";
        end
    otherwise
        G="q";
    end
end

if n*nt > 5000*5000
    warning("Donnees trop grandes pour affichage"); 
end

fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n")
fprintf(" ~   MERCI POUR VOTRE ATTENTION    ~")
fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n")