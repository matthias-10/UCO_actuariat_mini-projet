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
K = 0;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01/sqrt(S0);  % Variance partie fixe

Nd = 5;                  % Nombre des sous-intervalles 
 % verifier que Nd << n => a faire: ecrire un test
t0 = 0;                 % Debut de la periode
n = 2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode %% dès T=70 S a elements imaginaires?
nt = 1000;              % Nombre de trajectoires


starttime = datetime('now');
fprintf('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n');
fprintf('La programme a demarre a %s \n', starttime);
fprintf('%d -> Nombre de trajectoires \n', nt);
fprintf('%d -> Prix initial du sous jacent \n', S0)

%1% syms func(x) %1% requires Symbolic Math Toolbox.
%1% obligation(x) = S0*(1+r)^(x-t0);

%K = int(obligation,t0,T)/(T-t0);
%1% bonds_T = obligation(T);
%1% fprintf('%0.5g -> Prix d''une obligation a T\n', bonds_T)
fprintf('%0.5g -> Prix d''exercice de l''option \n', K);
fprintf(' . . . \n\n')
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
%size(S)
% >ans =      513        1000
C_inf = zeros(1,nt);
C_val = zeros(1,nt);
C_mat = zeros(1,nt);
for j = 1:nt
    S_vec = S(:,j);
    %% ~~~~~~~~~~~~~~~ calcul avec X_t ~~~~~~~~~~~~~~~~~ %%
    
    %%%%
    % debogage: C_inf = 0, sauf que la premiere et derniere valeur de C_inf

    % integral: l'aire de t0 a T sous S
    X_T = 0.5*S0 + sum(S_vec(2:n,:),1) + 0.5*S_vec(n+1,:);
    X_T = X_T/n; %ou (n+1)?

    C_inf_j = X_T - K .* ( X_T - K >= 0 );
    C_inf_0 = exp(-r*T)*C_inf_j;

    % ~ Estimateur ~
    % C_inf * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_inf]= C_inf(S_0)
    
    C_inf(j)=C_inf_0;


    %% ~~~~~~~~~ calcul avec X_t_prim (Valentin)~~~~~~~~ %%

    X_t_prim = sum(S_vec,1)/(n+1);
    %X_t_prim = mean(vecX_t_prim);
    C_N_j = X_t_prim - K .* ( X_t_prim - K >= 0 );

    % C_N * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_N]= C_N(S_0)
    C_N_0 = exp(-r*T)*C_N_j;

    C_val(j)=C_N_0;


    %% ~~~~~~~~~~ calcul avec X_t_prim (Matthias)~~~~~~~ %%

    %1/N * sum_1^N S_{kT/N}
    % => kT n'est pas un numero entier, il faut arrondir

    index = fliplr(1:n);
    warn_id = 'MATLAB:colon:nonIntegerIndex';
    warning('off', warn_id);
    % ^supprime Warning a cause de arrondir:
    index = index(1:(n/Nd):end); 
    X_t_matthias = sum(S_vec(index,:),1)/Nd;

    C_N_j = X_t_matthias - K .* ( X_t_matthias - K >= 0 );

    % C_N * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_N]= C_N(S_0)
    C_N_0 = exp(-r*T)*C_N_j;

    C_mat(j)=C_N_0;

end

%% ~~~~~~~~~~~~ affichage des estimateurs ~~~~~~~~~~~~~~ %%

% C_inf
C_inf_est = mean(C_inf);
C_inf_est_var = var(C_inf);

fprintf('L''estimateur du C_inf a t0 = %0.5g\n', ...
 C_inf_est);
fprintf('Son ecart type = %0.5g\n', sqrt(C_inf_est_var));

fprintf('\n methodes differentes pour C_N, \n premier Valentin, puis Matthias \n')
% Valentin C_inf
C_N_est_val = mean(C_val);
C_N_est_var_val = var(C_val);

fprintf('L''estimateur du C_N a t0 = %0.5g\n', ...
 C_N_est_val);
fprintf('Son ecart type = %0.5g\n', sqrt(C_N_est_var_val));

% Matthias C_inf
C_N_est_mat = mean(C_mat);
C_N_est_var_mat = var(C_mat);

fprintf('L''estimateur du C_N a t0 = %0.5g\n', ...
 C_N_est_mat);
fprintf('Son ecart type = %0.5g\n', sqrt(C_N_est_var_mat));


duree= toc;
fprintf('\n')
fprintf('%d trajectoires simules\n', nt);
fprintf('Avec S0 = %d, K = %0.5g \n', S0, K);
fprintf('L''integrale par (t_0 - T) de X_t, ');
    fprintf('le prix estime C(T) = %0.5g \n', C_inf_0);
fprintf('La moyenne des X_t: C(T) = %0.5g \n', C_N_0);
fprintf('Fini en %0.5g\n', duree);


%% ~~~~~~~~~~~~~~~~~~~~~ graphes ~~~~~~~~~~~~~~~~~~~~~~~ %%

% 1:   graphe de S; 
% 2-3: histogrammes de C_inf et C_N; 
% 4-5:   boxplot des estimateurs

fprintf('\n 1: graphe de S \n')
input('Tapez [Enter] pour afficher le graphe\n')

%axis([0 T 0.8*min(min(S)) 1.5*max(max(S))]) %x-axe limits

plot(t, S)

% pour comparison, si j'epargne pour le taux r:
%plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
%1% fplot(obligation, [t0 T], "-k"); 

legend("les prix S_t des actions", "sans risque","Location","northwest");


fprintf('\n 2: histogramme de C_inf \n')
input('Tapez [Enter] pour afficher le graphe\n')

% E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum {C(T)}
histogram( C_inf );
title("Histogramm des C(T) pour X_{infinie}");


fprintf('\n 3: histogramme de C_N \n')
input('Tapez [Enter] pour afficher le graphe\n')

histogram( C_mat );
title("Histogramm des C(T) pour X_{N}");


fprintf('\n 4: boxplot des estimateurs \n')
input('Tapez [Enter] pour afficher le graphe\n')

%tiledlayout(1,2)
%nexttile
%hold on 

boxplot( C_inf );
title('boxplot de C_{infinie} a T')
ylabel('C_T, valeurs actualisees')

%hold off
%nexttile
%hold on

boxplot ( C_mat );
title('boxplot de C_{N} a T')

%hold off
