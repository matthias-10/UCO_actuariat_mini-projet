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
K = 42;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01/sqrt(S0);  % Variance partie fixe

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
%1% fprintf('%0.5g -> Prix d''une obligation a T\n', bonds_T)
fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf(' . . .\n')
tic


%% ~~~~~~~~~~~~~~~~~~~~ Simulation ~~~~~~~~~~~~~~~~~~~~~ %%

dt = (T-t0)/n;
t = t0:dt:T;

S = zeros(n+1,nt);
S(1,:) = S0;

% Simulation pas a pas
for i = 2:(n+1)
    dW_t = normrnd(zeros(1,nt),sqrt(dt));
    dSi = S(i-1,:).*( r*dt + sigma*sqrt(S(i-1,:)).*dW_t );
    S(i,:) = S(i-1,:) + dSi;
end


%% ~~~~~~~~~~~~~~~ prix de l'option C ~~~~~~~~~~~~~~~~~~ %%

C_inf = zeros(1,nt);
C_N = zeros(1,nt);
for j = 1:nt
    S_vec = S(:,j);
    %% ~~~~~~~~~~~~~~~ calcul avec X_t ~~~~~~~~~~~~~~~~~ %%

    % integral: l'aire de t0 a T sous S
    X_T = 0.5*S0 + sum(S_vec(2:n,:),1) + 0.5*S_vec(n+1,:);
    X_T = X_T/n; %ou (n+1)?

    C_inf_j = X_T - K .* ( X_T - K >= 0 );
    C_inf_0 = exp(-r*T)*C_inf_j;

    % ~ Estimateur ~
    % C_inf * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_inf]= C_inf(S_0)
    
    C_inf(j)=C_inf_0;

    %% ~~~~~~~~~~~~~ calcul avec X_t_prim ~~~~~~~~~~~~~~ %%

    %1/N * sum_1^N S_{kT/N}
    % => kT n'est pas un numero entier, il faut arrondir

    index = fliplr(1:n);
    warn_id = 'MATLAB:colon:nonIntegerIndex';
    warning('off', warn_id);
    % ^supprime Warning a cause de arrondir:
    index = index(1:(n/Nd):end); 
    X_T_prim = sum(S_vec(index,:),1)/Nd;

    C_N_j = X_T_prim - K .* ( X_T_prim - K >= 0 );

    % C_N * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_N]= C_N(S_0)
    C_N_0 = exp(-r*T)*C_N_j;

    C_N(j)=C_N_0;

end

%% ~~~~~~~~~~~~ affichage des estimateurs ~~~~~~~~~~~~~~ %%

duree= toc;
fprintf('\n')
fprintf('%d trajectoires simules\n', nt);
fprintf('L''integrale par (t_0 - T) de X_t, ');
    fprintf('le prix estime C(T) = %0.5g \n', C_inf_0);
fprintf('La moyenne des X_t: C(T) = %0.5g \n', C_N_0);
fprintf('Fini en %0.5g\n', duree);
fprintf('\n')

fprintf('Les estimateurs Monte-Carlo:\n')

% C_inf
C_inf_est = mean(C_inf);
C_inf_est_var = var(C_inf);

fprintf('L''estimateur du C_inf a t0 = \n%0.5g\n', ...
 C_inf_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_inf_est_var));

% C_N
C_N_est = mean(C_N);
C_N_est_var = var(C_N);

fprintf(['L''estimateur du C_N a t0, avec ' ...
    '%d sous-intervalles = \n%0.5g\n'], ...
    Nd, C_N_est_mat);
fprintf('Son ecart type = %0.5g\n', sqrt(C_N_est_var));
%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%

% 1:   graphe de S; 
% 2-3: histogrammes de C_inf et C_N; 
% 4-5: boxplot des estimateurs

P = input( ...
    ['Pour afficher n''importe quelle graphique, tapez' ...
    'son numero <1-5>. \n' ...
    'Pour quitter tapez [q] ou plusieures fois [Enter]'] );
while G~="q"
P = input(' ');
switch P
    case 1
        disp('Flor')
        P=P+1;
    case 2
        disp('Circulo')
        P=P+1;
    case 3
        disp('Onda senoidal')
        P=P+1;
    case 4
        disp('Elipse')
        P=P+1;
    case 5
        disp('Nube')
        P=P+1;
    case 6
        P=input('')
        if ismember(P, 1:5)
            P = P
        else
            G="q"
        end
    otherwise
        G="q")
    end
end
fprintf('\n< 1: graphe de S >\n')
input('Tapez [Enter] pour afficher le graphe\n')

plot(t, S)

% pour comparison, si j'epargne pour le taux r:
%plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
%1% fplot(obligation, [t0 T], "-k"); 
legend("les prix S_t des actions", "sans risque",...
       "Location","northwest");


fprintf('\n< 2: histogramme de C_inf >\n')
input('Tapez [Enter] pour afficher le graphe\n')

% E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum {C(T)}
histogram( C_inf );
title("Histogramm des C(T) pour X_{infinie}");


fprintf('\n< 3: histogramme de C_N >\n')
input('Tapez [Enter] pour afficher le graphe\n')

histogram( C_N );
title("Histogramm des C(T) pour X_{N}");


fprintf('\nboxplot des estimateurs')
fprintf('\n< 4: boxplot de l''estimateur C_{infinie} >\n')
input('Tapez [Enter] pour afficher le graphe\n')

boxplot( C_inf );
title('boxplot de C_{infinie} a T')
ylabel('C_T, valeurs actualisees')


fprintf('\n< 5: boxplot de l''estimateur C_{N} >\n')
input('Tapez [Enter] pour afficher le graphe\n')

boxplot ( C_N );
title('boxplot de C_{N} a T')

