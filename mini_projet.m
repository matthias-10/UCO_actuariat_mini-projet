% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Valentin DE CRESPIN DE BILLY                      UTF-8 %
% Matthias LANG                                30.11.2021 %
% requires:                                               %
% - Statistics and Machine Learning Toolbox               %
% - Symbolic Math Toolbox                                 %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ~~~~~~ Mathematiques financieres: Mini-projet 1 ~~~~~~~ %

%% ~~~~~~~~~~~~~~~~~~~~ Parametres ~~~~~~~~~~~~~~~~~~~~~ %%

S0 = 40;                % Prix initial du sous jacent
K = 41;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode
Nd = 8;                 % Nombre des sous-intervalles 

nt = 1000;              % Nombre de trajectoires


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

if Nd > n/2-1
    warning("Le nombre des sous-intervalles est trop petit")
    fprintf('Il fallait Nd << n')
end

starttime = datetime('now');
fprintf('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n');
fprintf('La programme a demarre a %s \n', starttime);
fprintf('%d -> Prix initial du sous jacent \n', S0)

%1% syms func(x) %1% requires Symbolic Math Toolbox.
%1% obligation(x) = S0*(1+r)^(x-t0);
%K = int(obligation,t0,T)/(T-t0);
%1% bonds_T = obligation(T);
%1% fprintf('%0.5g -> Prix d''une obligation a T\n',bonds_T)

fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf(' . . .\n')
tic


%% ~~~~~~~~~~~~~~~~~~~~ Simulation ~~~~~~~~~~~~~~~~~~~~~ %%

dt = (T-t0)/n;
t = t0:dt:T;

% les premieres nt_a valeurs pour l'affichage plus tard
nt_a = 15;
S = zeros(nt_a, n+1);


%% ~~~~~~~~~~~~~~~~ prix de l'option C ~~~~~~~~~~~~~~~~~ %%

C_inf = zeros(nt,1);
C_N = zeros(nt,1);
X_inf_v = zeros(nt,1);
X_N_v = zeros(nt,1);
for j = 1:nt
    S_vec = zeros(1, n);
    S_vec(:, 1) = S0;
    

    %% ~~~~~~~~~~~~~ simuler pas a pas ~~~~~~~~~~~~~~~~~ %%
   
    for i = 2:(n+1)
        dW_t = normrnd(0,sqrt(dt));
        dSi = S_vec(i-1)* ...
              ( r*dt + sigma*sqrt(S_vec(i-1))*dW_t );
        S_vec(i) = S_vec(i-1) + dSi;
    end

    % sauvegarder les premieres nt_a actions
    if j <= nt_a
        S(j,:) = S_vec;
    end


    %% ~~~~~~~~~~~~ C_inf: calcul avec X_T ~~~~~~~~~~~~~ %%

    % integral: l'aire de t0 a T sous S
    X_T = 0.5*S0 + sum(S_vec(2:n)) + 0.5*S_vec(n+1);
    X_T = X_T/n; %ou (n+1)?

    C_inf_j = (X_T - K) * ( X_T - K >= 0 );
    C_inf_0 = exp(-r*T)*C_inf_j;

    % ~ Estimateur ~
    % C_inf * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_inf]= C_inf(S_0)
    
    X_inf_v(j)=X_T;
    C_inf(j)=C_inf_0;


    %% ~~~~~~~~~~ C_N: calcul avec X_T_prim ~~~~~~~~~~~~ %%

    %1/N * sum_1^N S_{kT/N}
    % => kT n'est pas un numero entier, il faut arrondir

    index = fliplr(1:n);
    warn_id = 'MATLAB:colon:nonIntegerIndex';
    warning('off', warn_id);
    % ^supprime Warning a cause de arrondir:
    index = index(1:(n/Nd):end); 
    X_T_prim = sum(S_vec(index))/Nd;

    C_N_j = (X_T_prim - K) * ( X_T_prim - K >= 0 );

    % C_N * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_N]= C_N(S_0)
    C_N_0 = exp(-r*T)*C_N_j;

    X_N_v(j)=X_T_prim;
    C_N(j)=C_N_0;

end


%% ~~~~~~~~~~~~~~~~~~~ estimateurs ~~~~~~~~~~~~~~~~~~~~~ %%

% C_inf
C_inf_est = mean(C_inf);
C_inf_est_var = var(C_inf)/nt; %/nt ?

% X_inf
X_inf_mu = mean(X_inf_v);

% C_N
C_N_est = mean(C_N);
C_N_est_var = var(C_N)/nt; %/nt?

% X_N
X_N_mu = mean(X_N_v);


%% ~~~~~~~ estimateurs et intervalle de confiance ~~~~~~ %%
%                 (seulement pour C_inf)                  %

alpha = 0.05; % niveau au risque

v = nt/(nt-1)*var(X_inf_v); %variance d'echantillonnage

%%% variable supossée normale
IC_gauss = [X_inf_mu - sqrt(v/nt)*norminv(1-alpha/2) ...
            X_inf_mu - sqrt(v/nt)*norminv(alpha/2) ];

%%% bootstrap
sims = 10^3;
y = zeros(1, sims);
for i = 1:sims
    y(i) = mean(randsample(X_inf_v,nt,true)) - X_inf_mu;
end

IC_boot =  [X_inf_mu - quantile(y,1-alpha/2) ...
            X_inf_mu - quantile(y,alpha/2) ];


%% ~~~~~~~~~~~~~~ variable de controle ~~~~~~~~~~~~~~~~~ %%

% a faire


%% ~~~~~~~~~~~~ affichage des estimateurs ~~~~~~~~~~~~~~ %%

duree= toc;
fprintf('\n')
fprintf('%d trajectoires simules\n', nt);
fprintf('Fini en %0.5g\n', duree);
fprintf('\n')

fprintf('Les estimateurs Monte-Carlo:\n')

fprintf('L''estimateur du C_inf a t0 = \n%0.5g\n', ...
 C_inf_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_inf_est_var));

fprintf(['L''estimateur du C_N a t0, avec ' ...
    '%d sous-intervalles = \n%0.5g\n'], ...
    Nd, C_N_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_N_est_var));


%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%

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
        plot(t, S)
        plot([t0 T],[K K], ':k', 'LineWidth',2)
        hold off
        % pour comparison, si j'epargne pour le taux r:
        %plot([t0 T], [S0 S0*(1+r)^(T-t0)],"--k"); %obl.
        %1% fplot(obligation, [t0 T], "-k"); 
        legend("K, le prix d''exercice", ...
               "les prix S_t des actions",...
               "Location","northwest");
        P=P+1; input('\n');

    case 2
        fprintf(['< 2: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_infinie >\n'])
        figure(1)
        % E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum{C(T)}
        %histogram( C_inf );
        ecdf( X_inf_v );
        hold on 
        plot([K K],[0 1], 'k')
        plot([min(X_inf_v) max(X_inf_v)], [.5 .5],':b')
        x = [min(X_inf_v):.1:max(X_inf_v)];
        nor = normcdf(x,X_inf_mu,v);
        plot(x,nor,':r')
        hold off
        legend("ecdf", "K", "P=50%", "cdf normal")
        title("ecdf X(T) pour X_{infinie}");
        P=P+1; input('\n');

    case 3
        fprintf(['< 3: fonction de distribution ' ...
            'cumulative estime'  ...
            '\n C(T) pour X_{infinie} de C_N >\n'])
        figure(1)
        ecdf( X_N_v );
        hold on 
        plot([K K],[0 1], 'k')
        plot([min(X_N_v) max(X_N_v)], [.5 .5],':b')
        hold off
        legend("ecdf", "K", "P=50%")
        title("ecdf X(T) pour X_{N}");
        P=P+1; input('\n');

    case 4
        fprintf(['< 4: boxplot de l''estimateur ' ...
                 'C_{infinie} >\n'])
        figure(1)
        boxplot( C_inf );
        title('boxplot de C_{infinie} a T')
        ylabel('C_T, valeurs actualisees')
        P=P+1; input('\n');

    case 5
        fprintf('< 5: boxplot de l''estimateur C_{N} >\n\n')
        figure(1)
        boxplot ( C_N );
        title('boxplot de C_{N} a T')
        P=P+1;

    case 6
        P=input([' ~ ~ ~ ~ ~ ~\n ' ...
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
