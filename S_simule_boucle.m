S0 = 40;                % Prix initial du sous jacent
K = 41;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 1000; %2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode
Nd = 8;                 % Nombre des sous-intervalles 

nt = 100000;              % Nombre de trajectoires



starttime = datetime('now');
fprintf('\n ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n');
fprintf('La programme a demarre a %s \n', starttime);
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
   
    % Methode de Euler
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
    X_T = X_T/n;

    C_inf_j = (X_T - K) * logical( X_T - K >= 0 );
    C_inf_0 = exp(-r*T)*C_inf_j;
    
    X_inf_v(j)=X_T;
    C_inf(j)=C_inf_0;


    %% ~~~~~~~~~~ C_N: calcul avec X_T_prim ~~~~~~~~~~~~ %%

    %1/N * sum_1^N S_{kT/N}
    % => kT n'est pas un numero entier, il faut arrondir

    index = fliplr(1:n);
    warn_id = 'MATLAB:colon:nonIntegerIndex';
    warning('off', warn_id);
    % ^supprime Warning a cause de arrondir:
    % l'erreur relatif à S_t dépend que de dt, alors avec dt desormais
    % petit pas de souci (d'ailleurs: interpolation lineaire, simuler,
    % choisir N n/a a partie de G
    index = index(1:(n/Nd):end); 
    X_T_prim = sum(S_vec(index))/Nd;

    C_N_j = (X_T_prim - K) * logical( X_T_prim - K >= 0 );
    C_N_0 = exp(-r*T)*C_N_j;

    X_N_v(j)=X_T_prim;
    C_N(j)=C_N_0;

end
endtime = datetime('now');
toc