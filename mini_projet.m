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

nt = 10000;            % Nombre de trajectoires

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
C_N_est = mean(C);
C_N_est_var = var(C)/nt; %/nt ?

C_N_IC_inf = C_N_est + sqrt(C_N_est_var)*norminv(alpha/2);
C_N_IC_sup = C_N_est + sqrt(C_N_est_var)*norminv(1-alpha/2);
L = C_N_IC_sup-C_N_IC_inf;
tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('pour X_prime: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',C_N_est),...
{' IC = ['},sprintf('%05.3f',C_N_IC_inf),...
{' , '},sprintf('%05.3f',C_N_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));


%% ~~~~~~~~~~~~~~~ variable antithetique ~~~~~~~~~~~~~~~ %%

tic

S = S0 * ones(nt, 1);
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
%co = cov([X Xa]);
%fprintf(['La covariance entre X et la variable '...
%        'antithetique est: %0.5g \n'], co(2,2))


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

C = exp(-r*T) * max(X-K,0);
A = cov(C, VC);
lambda = A(1,2)/A(2,2); 

% en utilisant les est. empiriques
% comparer avec les moments analytiques? -> voir pdf/VBA


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
C = exp(-r*T) * max(X-K,0);

Z = C - lambda * (VC - S0);

C_2_est = mean(Z);
C_2_est_var = var(Z)/nt;

C_2_IC_inf = C_2_est + sqrt(C_2_est_var)*norminv(alpha/2);
C_2_IC_sup = C_2_est + sqrt(C_2_est_var)*norminv(1-alpha/2);
L = C_2_IC_sup - C_2_IC_inf;

tps = toc;

% fonction d'affichage 

fprintf('\n')
fprintf('Variable de controle - somme de dW_t: \n');

disp(strcat(...
{' C = '},sprintf('%05.3f',C_2_est),...
{' IC = ['},sprintf('%05.3f',C_2_IC_inf),...
{' , '},sprintf('%05.3f',C_2_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',L),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',L * sqrt(tps))));

%p = corr(X, VC); 
% optimum: lambda =~ corr(X,Y)*(Var(X)/Var(Y))^.5
% efficace ?
%fprintf(['La correlation entre X et la variable '...
%        'de controle est: %0.5g\n'], p)


%% ~~~~~~~~~~~~~ variable de controle 3 ~~~~~~~~~~~~~~~~ %%

% mouvement brownien de la meme signe

%dWt_vc = normrnd(zeros(nt,1),sqrt(dt));
%dWt_vc = dWt_vc + dWt_vc .*sign(dWt).*(sign(dWt_vc) - sign(dWt));
%dSi = VC(:,i-1).* ...
%      ( r*dt + sigma*sqrt(S(:,i-1)).*dWt_vc );
%VC(:,i)      =    VC(:,i-1) + dSi;
%VC(:,i) = VC(:,i) .* (VC(:,i) >= 0);


%% ~~~~~~~~~~~~~ variable de controle 4 ~~~~~~~~~~~~~~~~ %%

% calculer des option europeenes, mais avec des actions:
% S = S .*(1 +r*dt + sigma*dW_t );
% alors qu'on puisse calculer leur payoff avec Balckes-Sch.


%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%

% 1:   graphes de S; 
% 2-3: ecdf de C_inf et C_N; 
% 4-5: boxplot des estimateurs
% 6:   deux graphiques qui demontrent une problematique
% 7:   intervalles de confiance

G = "g";
%G = "q" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P = input(['\n' ...
%    'Pour afficher n''importe quel graphique, tapez ' ...
%    'son numero <1-8> ou [Enter]. \n' ...
%    'Pour quitter tapez plusieures fois [Enter]:\n'] );
P=8;

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
        
        w = 3; % nombre des IC affichees
        %fprintf('< 8: ICs (normales) >')
        plot([C_est Ca_est C_2_est], 1:w, 'x')
        %line([K K],[0 5],'Color','green','LineStyle','--')
       
        line([C_IC_inf C_IC_sup],[1 1])
        line([Ca_IC_inf Ca_IC_sup],[2 2])
        line([C_2_IC_inf C_2_IC_sup],[3 3])
        
        legend("estimateurs",'C','C_a','C_{VC2}')
        L1 = C_IC_sup - C_IC_inf;
        L2 = C_est - K;
        limf = [C_IC_inf C_IC_sup] + max(L1,L2)*[-1 1];
        %xlim(limf)
        ylim([0 w+1])
        yticks(1:w)

        yticklabels({'C','C_a','C_{VC2}'})
        %title('Intervalles de confiance (sauf C)')
        
        %P=P+1;
        G="q";
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

%fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
%fprintf(  " ~   MERCI POUR VOTRE ATTENTION    ~")
fprintf("\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n")
