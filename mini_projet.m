% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Valentin DE CRESPIN DE BILLY                      UTF-8 %
% Matthias LANG                                30.11.2021 %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ~~~~~ Mathematiques financieres: Mini-projet 1 ~~~~~~~~ %

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
   

    %% ~~~~~~~~~~~~ C_inf: calcul avec X_T ~~~~~~~~~~~~~ %%

    % integral: l'aire de t0 a T sous S
    X_T = 0.5*S0 + sum(S_vec(2:n,:),1) + 0.5*S_vec(n+1,:);
    X_T = X_T/n; %ou (n+1)?

    C_inf_j = (X_T - K) .* ( X_T - K >= 0 );
    C_inf_0 = exp(-r*T)*C_inf_j;

    % ~ Estimateur ~
    % C_inf * exp(-rT) est une martingale donc 
    % E[exp(-rT)*C_inf]= C_inf(S_0)
    
    C_inf(j)=C_inf_0;


    %% ~~~~~~~~~~ C_N: calcul avec X_T_prim ~~~~~~~~~~~~ %%

    %1/N * sum_1^N S_{kT/N}
    % => kT n'est pas un numero entier, il faut arrondir

    index = fliplr(1:n);
    warn_id = 'MATLAB:colon:nonIntegerIndex';
    warning('off', warn_id);
    % ^supprime Warning a cause de arrondir:
    index = index(1:(n/Nd):end); 
    X_T_prim = sum(S_vec(index,:),1)/Nd;

    C_N_j = (X_T_prim - K) .* ( X_T_prim - K >= 0 );

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

% 1:   graphes de S; 
% 2-3: histogrammes de C_inf et C_N; 
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
        fprintf('< 1: graphes de S >\n')
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
        fprintf('< 2: histogramme de C_inf >\n')
        figure(1)
        % E_\pi (e^-rT (X_T - K)^+ / F_O) ~ 1/nt \sum{C(T)}
        histogram( C_inf );
        title("Histogramm des C(T) pour X_{infinie}");
        P=P+1; input('\n');

    case 3
        fprintf('< 3: histogramme de C_N >\n')
        figure(1)
        histogram( C_N );
        title("Histogramm des C(T) pour X_{N}");
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
