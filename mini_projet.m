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
r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

n = 2^6;                % Nombre de intervalles
T = 5;                  % Fin de la periode/exercice = tau
Nd = 8;                 % Nombre des sous-intervalles 

nt = 1000;            % Nombre de trajectoires

alpha = 0.05;           % niveau au risque


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

nlambda = 500; % sim.s pour determiner le lambda des VC

if Nd > n/2-1
    warning("Le nombre de sous-intervalles est tres petit")
    fprintf('Il fallait Nd << n')
end

tic
starttime = datetime('now');
fprintf('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n');
fprintf('La programme a demarre a %s \n', starttime);

% K base sur le prix moyen d'une obligation sans risque
syms func(x) 
obligation(x) = S0*(1+r)^x;
K = double( int(obligation,0,T)/T);
bonds_T = obligation(T);

fprintf('%d -> Prix initial du sous jacent \n', S0)
fprintf('%0.5g -> Prix univers risque neutre a T\n',bonds_T)
fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf('calculation en cours . . .\n')

dt = T/n;
t = 0:dt:T;


%% ~~~~~~~~~~~~~~~~~ Monte-Carlo pur ~~~~~~~~~~~~~~~~~~~ %%

tic

S = S0 * ones(nt, 1);
X = S/2;

for i = 2:(n+1)
    dW_t = normrnd(0, sqrt(dt), nt, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
end

X = (X - S/2)/n;

C = exp(-r*T) * max(X-K,0);

% ~ Estimateur ~
% C_inf * exp(-rT) est une martingale donc 
% E[exp(-rT)*C_inf]= C_inf(S_0)

% C
C_est = mean(C);
C_est_var = var(C)/nt; %/nt ?

C_IC_inf = C_est + sqrt(C_est_var)*norminv(alpha/2);
C_IC_sup = C_est + sqrt(C_est_var)*norminv(1-alpha/2);
L = C_IC_sup-C_IC_inf;
tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('%d trajectoires simules\n', nt);

fprintf('\n')
fprintf('estimateur Monte-Carlo: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',C_est),...
{' IC = ['},sprintf('%05.3f',C_IC_inf),...
{' , '},sprintf('%05.3f',C_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));


%% ~~~~~~~~~~~~ C_N: calcul avec X_T_prim ~~~~~~~~~~~~~~ %%

tic

S = S0 * ones(nt, 1);
X = zeros(nt, 1);
l=1;

for i = 2:(n+1)
    dW_t = normrnd(0, sqrt(dt), nt, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    if (i/n) > (l/Nd)
        X = X+S;
        l = l+1;
    end
end

X = (X+S)/Nd;
C = exp(-r*T) * max(X-K,0);

% C
C_est = mean(C);
C_est_var = var(C)/nt; %/nt ?

C_IC_inf = C_est + sqrt(C_est_var)*norminv(alpha/2);
C_IC_sup = C_est + sqrt(C_est_var)*norminv(1-alpha/2);
L = C_IC_sup-C_IC_inf;
tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('pour X_prime: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',C_est),...
{' IC = ['},sprintf('%05.3f',C_IC_inf),...
{' , '},sprintf('%05.3f',C_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));


%% ~~~~~~~~~~~~~~~ variable antithetique ~~~~~~~~~~~~~~~ %%

tic

S = S0 * ones(nb_trajectoires,1);
Sa = S;
X = S/2;
Xa = X;

% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nt, 1);
    S = S .* (1 + r*dt + sigma*sqrt(abs(S)).*dW_t );
    Sa = Sa .* (1 + r*dt - sigma*sqrt(abs(Sa)).*dW_t );
    X = X + S;
    Xa = Xa + Sa;
end

X = (X - S/2)/n;
Xa = (Xa - Sa/2)/n;

Ca = exp(-r*T)*(max(Xa - K,0) + max(X-K,0))/2;
Ca_est = mean(Ca);
Ca_est_var = var(Ca)/(2*nt);

Ca_IC_inf = Ca_est + sqrt(Ca_est_var)*norminv(alpha/2);
Ca_IC_sup = Ca_est + sqrt(Ca_est_var)*norminv(1-alpha/2);

L = Ca_IC_sup-Ca_IC_inf;
tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('avec une variable antithetique: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',Ca_est),...
{' IC = ['},sprintf('%05.3f',Ca_IC_inf),...
{' , '},sprintf('%05.3f',Ca_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));


% efficace ?
co = cov([X Xa]);
fprintf(['La covariance entre X et la variable '...
        'antithetique est: %0.5g'], co(2,2))


%% ~~~~~~~~~~~~~ variable de controle 1 ~~~~~~~~~~~~~~~~ %%

% Variable de controle - aire sous un mouvement brownien



%% ~~~~~~~~~~~~~ variable de controle 2 ~~~~~~~~~~~~~~~~ %%

% Variable de controle - somme de dW_t: 
% correspond a un mouvement brownien
% ~~~~~~~~~~~~~~~~ Calculer lambda ~~~~~~~~~~~~~~~~~~~~~~ %

S = S0 * ones(nlambda,1);
VC = S0 * ones(nlambda, 1);
X = S/2;
for i = 2:(n+1)
    dW_t = normrnd(0, sqrt(dt), nlambda, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t; % succesion des petits accroissements du brownien => VC et S sont correles
end
X = (X - S/2)/n;

A = cov(X, VC);
lambda = A(1,2)/A(2,2);

% Deutsch a utilise le payoff au lieu de X


% ~~~~~~~~~~~~~~~~~~~~~ simuler ~~~~~~~~~~~~~~~~~~~~~~~~~ %

tic;

S = S0 * ones(nt,1);
VC = S0 * ones(nt, 1);
X = S/2;
for i = 2:(n+1)
    dW_t = normrnd(0, sqrt(dt), nt, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t;
end
X = (X - S/2)/n;

Z = X - lambda * (VC - S0);
%%%%%
C = exp(-r*T) * max(Z-K,0);

C_est = mean(C);
C_est_var = var(C)/nt;

C_IC_inf = C_est + sqrt(C_est_var)*norminv(alpha/2);
C_IC_sup = C_est + sqrt(C_est_var)*norminv(1-alpha/2);
L = C_IC_sup - C_IC_inf;

tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('Variable de controle - somme de dW_t: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',C_est),...
{' IC = ['},sprintf('%05.3f',C_IC_inf),...
{' , '},sprintf('%05.3f',C_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));


%% ~~~~~~~~~~~~~ variable de controle 3 ~~~~~~~~~~~~~~~~ %%
    dWt_vc = normrnd(zeros(nt,1),sqrt(dt));
    dWt_vc = dWt_vc + dWt_vc .*sign(dWt).*(sign(dWt_vc) - sign(dWt));
    dSi = VC(:,i-1).* ...
          ( r*dt + sigma*sqrt(S(:,i-1)).*dWt_vc );
    VC(:,i)      =    VC(:,i-1) + dSi;
    VC(:,i) = VC(:,i) .* (VC(:,i) >= 0);


%% ~~~~~~~~~~~~~~ C_inf: calcul avec X_T ~~~~~~~~~~~~~~~ %%

% integral: l'aire de 0 a T sous S
X = (0.5*S0 + sum(S(:,2:n),2) + 0.5*S(:,n+1))/n;

C = (X - K) .* logical( X - K >= 0 );
C_0 = exp(-r*T)*C;

% avec la variable antithetique
X_a = (0.5*S0 + sum(S_anti(:,2:n),2) + 0.5*S_anti(:,n+1))/n;
% avec la variable de controle
X_vc = (0.5*S0 + sum(VC(:,2:n),2) + 0.5*VC(:,n+1))/n;






%% ~~~~~~~~~~~~~~~~~~~ estimateurs ~~~~~~~~~~~~~~~~~~~~~ %%
    
% ~ Estimateur ~
% C_inf * exp(-rT) est une martingale donc 
% E[exp(-rT)*C_inf]= C_inf(S_0)

% C
C_0_est = mean(C_0);
C_est_var = var(C_0)/nt; %/nt ?

X_mu = mean(X);

% C_N
C_0_prim_est = mean(C_0_prim);
C_prim_est_var = var(C_0_prim)/nt; %/nt?

X_prim_mu = mean(X_prim);


%% ~~~~~~~~~~~~~~ intervalle de confiance ~~~~~~~~~~~~~~ %%
%                 (seulement pour N=inf)                  %

v = nt/(nt-1)*var(X); % variance d'echantillonnage

%%% variable supossee normale
X_IC_gauss = [X_mu + sqrt(v/nt)*norminv(alpha/2) ...
              X_mu + sqrt(v/nt)*norminv(1-alpha/2) ];

% variable antithetique
X_ab_mu = mean([X;X_a]);
na = 2*nt;
va = na/(na-1)*var([X;X_a]);
X_a_IC_gauss = [X_ab_mu + sqrt(va/na)*norminv(alpha/2) ...
                X_ab_mu + sqrt(va/na)*norminv(1-alpha/2)];

% efficace ?
co = cov([X X_a]);
fprintf(['\nLa covariance entre X et la variable '...
        'antithetique est: %0.5g\n'], co(2,2))

%%% bootstrap pour C
sims = 10^3;
y = zeros(1, sims);
for i = 1:sims
    y(i) = mean(randsample(C_0_est,nt,true)) - C_0_est;
end

C_IC_boot =  [C_0_est + quantile(y,alpha/2) ...
              C_0_est + quantile(y,1-alpha/2) ];


%% ~~~~~~~~~~~~~~ variable de controle ~~~~~~~~~~~~~~~~~ %%
%                (seulement pour X_inf)                   %


% on pourrai utiliser au lieu de la var. antithetique la
% variable de controle suivante et vice-versa

% E(Y) ~= E(X) ~= E(Z) =~ mean(X_a)

EY_vc = mean(X_vc);
X_vc;

p = corr(X, X_vc); 
% optimum: lambda =~ corr(X,Y)*(Var(X)/Var(Y))^.5
lambda = p*(var(X)/var(X_vc))^.5;
Z_vc = X - lambda * (X_vc - EY_vc);

na = nt;
va = na/(na-1)*var(Z_vc);

Z_IC_gauss = [EY_vc + sqrt(va/na)*norminv(alpha/2) ...
              EY_vc + sqrt(va/na)*norminv(1-alpha/2)];

% efficace ?
fprintf(['\nLa correlation entre X et la variable '...
        'de controle X_vc est: %0.5g\n'], p)


% utilisant X_a

EY_a = mean(X_a);
Y = 2*EY_a - X_a; % X_mu + EY_a - X_a ???

p = corr(X, Y); 
%bien entendu, les deux sont au-peu-pres 1 correles

% optimum: lambda =~ corr(X,Y)*(Var(X)/Var(Y))^.5
lambda = p*(var(X)/var(Y))^.5;
Z_a = X - lambda * (Y - EY_a);

na = nt;
va = na/(na-1)*var(Z_a);

Z_a_IC_gauss = [EY_a + sqrt(va/na)*norminv(alpha/2) ...
                EY_a + sqrt(va/na)*norminv(1-alpha/2)];

% efficace ?
fprintf(['\nLa correlation entre X et la variable '...
        'de controle X_a est: %0.5g\n'], corr(X,X_a))

%% ~~~~~~~~~~~~ affichage des estimateurs ~~~~~~~~~~~~~~ %%

duree= toc;

fprintf('Fini en %0.5g secondes\n', duree);
fprintf('\n')

fprintf(' ~ Les estimateurs Monte-Carlo: ~ \n')

fprintf('L''estimateur du C_inf ajd = %0.5g\n', ...
 C_0_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_est_var));


fprintf(['L''estimateur du C_N ajd, avec ' ...
    '%d sous-intervalles = %0.5g\n'], ...
    Nd, C_0_prim_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_prim_est_var));

% fprintf('\n ~ Des intervalles de Confiance ~ \n');
% fprintf('L''intervalle de confiance de X (normal):\n');
% X_IC_gauss
% fprintf('La meme intervalle avec var. antithetiques:\n');
% X_a_IC_gauss
% fprintf(['L''intervalle de confiance '...
%         'de Z avec VC (normal):\n']);
% Z_IC_gauss
% fprintf(['L''intervalle de confiance '...
%         'de Z avec X_a (normal):\n']);
% Z_a_IC_gauss
fprintf('L''intervalle de confiance de C_N (bootstrap):\n');
C_IC_boot


%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%


% 1:   graphes de S; 
% 2-3: ecdf de C_inf et C_N; 
% 4-5: boxplot des estimateurs
% 6:   deux graphiques qui demontrent une problematique
% 7:   intervalles de confiance

G = "g";
P = input(['\n' ...
    'Pour afficher n''importe quel graphique, tapez ' ...
    'son numero <1-8> ou [Enter]. \n' ...
    'Pour quitter tapez plusieures fois [Enter]:\n'] );

if isstring(P) || isempty(P)
    P = 1;
else 
    if ~ismember(P,1:8)
        P = 1;
    end
end

while G~="q"
    disp("[Enter] pour continuer")
    switch P
    case 1
        fprintf('< 1: quelques premiers graphes de S >')
        figure(1)
        plot([0 T],[K K], ':k', 'LineWidth',2)
        hold on
        plot(t, obligation(t))
        plot(t, S_aff,'b')
        plot(t, S_anti_aff, 'r:')
        plot(t, VC_aff(1:nt_a,:),'g--')
        % probleme si nt < nt_a
        plot([0 T],[K K], '--k', 'LineWidth',2)
        hold off
        % pour comparison, si j'epargne pour le taux r:
        %plot([0 T], [S0 S0*(1+r)^T],"--k"); %obl.
        %1% fplot(obligation, [0 T], "-k"); 
        xlabel("t")
        legend("K, le prix d''exercice", ...
               "obligation (sans risque)", ...
               "les prix S_t des actions",...
               "les variables antithetiques",...
               "les variables de controle",...
               "Location","northwest");
        
        if n*nt > 100*1000000; P=7; end
        P=P+1; input('\n\n');
    
    case 2
        if n*nt > 100*1000000; G="q"; end
        
        fprintf(['< 2: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_infinie >'])
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
        
        P=P+1; input('\n\n');

    case 3
        if n*nt > 100*1000000; G="q"; end
        
        fprintf(['< 3: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_N >'])
        figure(1)
        ecdf( X_prim );
        hold on 
        plot([K K],[0 1], 'k')
        plot([min(X_prim) max(X_prim)], [.5 .5],':b')
        hold off
        legend("ecdf", "K", "P=50%")
        title("ecdf X(T) pour X_{N}");
        
        P=P+1; input('\n\n');

    case 4
        if n*nt > 100*1000000; G="q"; end
        
        fprintf(['< 4: boxplot de l''estimateur ' ...
                 'C_{infinie} >'])
        figure(1)
        boxplot( C_0 );
        xticks({})
        title('boxplot de C_{infinie} a T')
        ylabel('C_T, valeurs actualisees')
        
        P=P+1; input('\n\n');

    case 5
        if n*nt > 100*1000000; G="q"; end
        
        fprintf('< 5: boxplot de l''estimateur C_{N} >')
        figure(1)
        boxplot ( C_0_prim );
        xticks({})
        title('boxplot de C_{N} a T')
        
        P=P+1; input('\n\n');
    
    case 6
        if n*nt > 100*1000000; G="q"; end
        
        fprintf('< 6: L''IC de la variable de controle')
        fprintf('\n suivant pour Z a aide de VC')
        
        plot(sort(Z_vc))
        hold on 
        plot(sort(X))
        plot([1 na],[K K], '--k', 'LineWidth',1)
        hold off
        title("X vs variable de controle Z")
        xlabel("nt")
        legend("Z","X","K")
        
        input('\n... 6.5 < scatter >');
        
        scatter(X,X_vc);
        hold on; 
        plot([min(X) max(X)],[min(X) max(X)],'-k');
        plot(X_mu,EY_vc,'*r','LineWidth',2);
        legend("X-X_{vc} en pair",...
               "X=X_{vc}",...
               "les moyennes"); 
        hold off
        xlabel("X")
        ylabel("X_{vc}")
        
        P=P+1; input('\n');

    case 7
        if n*nt > 100*1000000; G="q"; end
        
        fprintf('< 7: L''IC de la variable de controle ')
        fprintf('\n suivant pour Z a aide de X_a')
        
        plot(sort(Z_a))
        hold on 
        plot(sort(X))
        plot([1 na],[K K], '--k', 'LineWidth',1)
        hold off
        title("X vs variable de controle Z")
        xlabel("nt")
        legend("Z","X","K")
        
        input('\n... 7.5 < scatter >');
        
        scatter(X,Y);
        hold on; 
        plot([min(X) max(X)],[min(X) max(X)],'-k');
        plot(X_mu,EY_a,'*r','LineWidth',2);
        legend("X-Y en pair","X=Y","les moyennes"); 
        hold off
        xlabel("X")
        ylabel("Y avec laquelle la v.c. est construite")
        
        P=P+1; input('\n');
    case 8
        
        fprintf('< 8: ICs (normales) >')
        plot([X_mu X_ab_mu EY_vc EY_a],[1 2 3 4], 'x')
        line([K K],[0 5],'Color','green','LineStyle','--')
       
        line(X_IC_gauss,[1 1])
        line(X_a_IC_gauss,[2 2])
        line(Z_IC_gauss,[3 3])
        line(Z_a_IC_gauss,[4 4])

        legend("estimateurs","K","Z_a","Z_{vc}","X_a","X")
        L1 = X_IC_gauss(2)-X_IC_gauss(1);
        L2 = X_mu - K;
        limf = X_IC_gauss + max(L1,L2)*[-1 1];
        xlim(limf)
        ylim([0 5])
        yticks(1:4)

        yticklabels({'X','X_a', 'Z (vc)', 'Z (a)'})
        title('Intervalles de confiance (sauf C)')

        P=P+1;
    case 9
        if n*nt > 100*1000000; G="q"; end
        
        P=input(['\n ' ...
            'Pour afficher n''importe quel graphique, ' ...
            'tapez son numero <1-8> \n']);
        if ismember(P, 1:8)
            fprintf("Vous avez choisi: ")
        else
            G="q";
        end
    otherwise
        G="q";
    end
end

if n*nt > 100*1000000
    warning("Donnees trop grandes pour affichage"); 
end

fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
fprintf(  " ~   MERCI POUR VOTRE ATTENTION    ~")
fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")

% shows that X follow maybe more of a lognormal distr
histogram(X)
line([X_mu X_mu], [0 4500])


mu = 0;
sigma=.25;
M = normrnd(ones(1000),va);
histogram(X,'Normalization', 'pdf');
line([0 0], [0 0.1])
hold on
histogram(exp(M),'DisplayStyle','stairs','Normalization', 'pdf');
hold off
