
  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Valentin DE CRESPIN DE BILLY                      UTF-8 %
% Matthias LANG                                30.11.2021 %
% exige:                                                  %
% - Statistics and Machine Learning Toolbox               %
% - Symbolic Math Toolbox                                 %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% ~~~~~~ Mathematiques financieres: Mini-projet 1 ~~~~~~~ %

%%%%%  PARAMETRES SOUS-JACENT %%%%%
S0 = 40;                % Prix initial du sous jacent
K = 41;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Volatilité du sous jacent

t0 = 0;                 % Debut de la periode
T = 2;                  % Fin de la periode

n = 1000;                % Nombre de intervalles

%%%%%  PARAMETRES APPLICATION NUMERIQUE %%%%%

nlambda = 500;          % nombre de sous jacent simulé pour le calcul du lambda
n_simu = 1000000;       % nombre de sous jacent simulé

nb = 700;               % nombre de valeur conservé pour la calcul du cours moyen de manière discrète


%%%%       /!\ ne pas modifier le code en dessous /!\    

%%%%%  INITILISATION %%%%%

dt = (T-t0)/n;
t = t0:dt:T;


%% ~~~~~~~~~~~~ Simulation de 15 trajectoires ~~~~~~~~~~~~~~ %%

S = zeros(15,n+1);
S(:,1) = S0;

% Simulation pas a pas
for i = 2:(n+1)
    dW_t = normrnd(zeros(15,1),sqrt(dt));
    dSi = S(:,i-1).*( r*dt + sigma*sqrt(S(:,i-1)).*dW_t );
    S(:,i) = S(:,i-1) + dSi;
end

plot(t, S)





%% ~~~~~~~~~~~~ Méthode de Monte Carlo ~~~~~~~~~~~~~~ %%
nb_trajectoires = n_simu;

tic;
%t_d = datetime('now');
S = S0 * ones(nb_trajectoires,1);
X = S/2;

% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0,sqrt(dt), nb_trajectoires,1);
    S = S .* (1 + r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
end

X = X - S/2;

XT = X / n;

XK = max(XT-K,0);
MC = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
MC_IC_inf = MC - 1.96 * sd / sqrt(nb_trajectoires);
MC_IC_sup = MC + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'MC   : '},...
{' P = '},sprintf('%05.3f',MC),...
{' IC = ['},sprintf('%05.3f',MC_IC_inf),...
{' , '},sprintf('%05.3f',MC_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',MC_IC_sup - MC_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(MC_IC_sup - MC_IC_inf)  *sqrt(tps))...
));


%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec réduction de la variance - variables antithétiques ~~~~~~~~~~~~~~ %%

tic
%t_d = datetime('now');
S = S0 * ones(nb_trajectoires,1);
Sa = S;
X = S/2;
Xa = X;

% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0,sqrt(dt), nb_trajectoires,1);
    S = S .* (1 + r*dt + sigma*sqrt(abs(S)).*dW_t );
    Sa = Sa .* (1 + r*dt - sigma*sqrt(abs(Sa)).*dW_t );
    X = X+S;
    Xa = Xa + Sa;
end

X = X - S/2;
Xa = Xa - Sa/2;

XT = X / n;
XTa = Xa/n;

XK = max(XT-K,0);
XKa = max(XTa - K,0);
XK = (XK + XKa)/2;
ANT = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
A_IC_inf = ANT - 1.96 * sd / sqrt(nb_trajectoires);
A_IC_sup = ANT + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'Ant  : '},...
{' P = '},sprintf('%05.3f',ANT),...
{' IC = ['},sprintf('%05.3f',A_IC_inf),...
{' , '},sprintf('%05.3f',A_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',A_IC_sup - A_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(A_IC_sup - A_IC_inf)  *sqrt(tps))...
));


%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec reduction de la variance - variables de contrôle 1 ~~~~~~~~~~~~~~ %%

tic;
%t_d = datetime('now');
nb_trajectoires = nlambda;
S = S0 * ones(nb_trajectoires,1);
X = S/2;
VC = zeros(nb_trajectoires, 1);


% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0,sqrt(dt), nb_trajectoires,1);
    VC = VC + dW_t; % succesion des petits accroissements du brownien => VC et S sont correles
    S = S .* (1 + r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
end

X = X - S/2;

XT = X / n;

XK = max(XT-K,0);
A = cov(XK, VC);
lambda = A(1,2)/A(2,2);

nb_trajectoires = n_simu;
S = S0 * ones(nb_trajectoires,1);
X = S/2;
VC = zeros(nb_trajectoires, 1);


% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0,sqrt(dt), nb_trajectoires,1);
    VC = VC + dW_t; % succesion des petits accroissements du brownien => VC et S sont correles
    S = S .* (1 + r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
end

X = X - S/2;
XT = X / n;

XK = max(XT-K,0) - lambda * (VC - 0);
VC1 = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
VC1_IC_inf = VC1 - 1.96 * sd / sqrt(nb_trajectoires);
VC1_IC_sup = VC1 + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'VC1  : '},...
{' P = '},sprintf('%05.3f',VC1),...
{' IC = ['},sprintf('%05.3f',VC1_IC_inf),...
{' , '},sprintf('%05.3f',VC1_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',VC1_IC_sup - VC1_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(VC1_IC_sup - VC1_IC_inf)  *sqrt(tps))...
));


%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec reduction de la variance - variables de contrôle 2 ~~~~~~~~~~~~~~ %%

tic;
nb_trajectoires = nlambda;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;
VC_aire = S/2;
for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t;
    VC_aire = VC_aire + VC;
    % succesion des petits accroissements du brownien 
    % => VC et S sont correles
end

X = X - S/2;
XT = X / n;

VC_aire = VC_aire - VC/2;
VC_aire_T = VC_aire / n;

XK = max(XT-K,0);
A = cov(XK, VC_aire_T);
lambda = A(1,2)/A(2,2);



nb_trajectoires = n_simu;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;
VC_aire = S/2;
for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t;
    VC_aire = VC_aire + VC;
    % succesion des petits accroissements du brownien 
    % => VC et S sont correles
end

X = X - S/2;
XT = X / n;

VC_aire = VC_aire - VC/2;
VC_aire_T = VC_aire / n;

XK = max(XT-K,0)-lambda * (VC_aire_T - S0);
VC2 = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
VC2_IC_inf = VC2 - 1.96 * sd / sqrt(nb_trajectoires);
VC2_IC_sup = VC2 + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'VC2  : '},...
{' P = '},sprintf('%05.3f',VC2),...
{' IC = ['},sprintf('%05.3f',VC2_IC_inf),...
{' , '},sprintf('%05.3f',VC2_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',VC2_IC_sup - VC2_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(VC2_IC_sup - VC2_IC_inf)  *sqrt(tps))...
));


%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec reduction de la variance - variables de contrôle 3 ~~~~~~~~~~~~~~ %%
tic;
nb_trajectoires = nlambda;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;

for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC .* (1 +r*dt + dW_t/10 );
end

X = X - S/2;
XT = X / n;

XK = max(XT-K,0);
A = cov(XK, VC);
lambda = A(1,2)/A(2,2);



nb_trajectoires = n_simu;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;

for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC .* (1 +r*dt + dW_t/10 );
end

X = X - S/2;
XT = X / n;

E_VC = S0 * exp(r*T);

XK = max(XT-K,0)-lambda * (VC - E_VC);
VC3 = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
VC3_IC_inf = VC3 - 1.96 * sd / sqrt(nb_trajectoires);
VC3_IC_sup = VC3 + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'VC3  : '},...
{' P = '},sprintf('%05.3f',VC3),...
{' IC = ['},sprintf('%05.3f',VC3_IC_inf),...
{' , '},sprintf('%05.3f',VC3_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',VC3_IC_sup - VC3_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(VC3_IC_sup - VC3_IC_inf)  *sqrt(tps))...
));




%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec reduction de la variance - variables de contrôle 4 ~~~~~~~~~~~~~~ %%
tic;
nb_trajectoires = nlambda;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;

for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC .* (1 +r*dt + sigma * sqrt(S0) * dW_t );
end

X = X - S/2;
XT = X / n;

VC = max(VC-K,0);
XK = max(XT-K,0);
A = cov(XK, VC);
lambda = A(1,2)/A(2,2);



nb_trajectoires = n_simu;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;

for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC .* (1 +r*dt + sigma * sqrt(S0) * dW_t );
end

X = X - S/2;
XT = X / n;

z0 = (log(K/S0)-T*(r-0.5*(sigma*sqrt(S0))^2))/((sigma*sqrt(S0))*sqrt(T));
d1 = (sigma*sqrt(S0)) * T^(1/2)-z0;
d2 = -z0;
E_VC = S0 * exp(T*(r-0.5*(sigma*sqrt(S0))^2)+0.5*T*(sigma*sqrt(S0))^2)*normcdf(d1) - K * normcdf(d2);

VC = max(VC-K,0);
XK = max(XT-K,0)-lambda * (VC - E_VC);
VC4 = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
VC4_IC_inf = VC4 - 1.96 * sd / sqrt(nb_trajectoires);
VC4_IC_sup = VC4 + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'VC4  : '},...
{' P = '},sprintf('%05.3f',VC4),...
{' IC = ['},sprintf('%05.3f',VC4_IC_inf),...
{' , '},sprintf('%05.3f',VC4_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',VC4_IC_sup - VC4_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(VC4_IC_sup - VC4_IC_inf)  *sqrt(tps))...
));

%% ~~~~~~~~~~~~ Méthode de Monte Carlo avec reduction de la variance - variables de contrôle 5 ~~~~~~~~~~~~~~ %%

tic;
nb_trajectoires = nlambda;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;
for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t;
    % succesion des petits accroissements du brownien 
    % => VC et S sont correles
end

X = X - S/2;
XT = X / n;

VC = exp(VC/S0);

XK = max(XT-K,0);
A = cov(XK, VC);
lambda = A(1,2)/A(2,2);



nb_trajectoires = n_simu;
S = S0 * ones(nb_trajectoires,1);
VC = S0 * ones(nb_trajectoires, 1);
X = S/2;
for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    X = X+S;
    VC = VC + dW_t;
    % succesion des petits accroissements du brownien 
    % => VC et S sont correles
end

X = X - S/2;
XT = X / n;

VC = exp(VC/S0);

XK = max(XT-K,0)-lambda * (VC - exp(T/(2*S0^2)+1));
VC5 = exp(-r * T ) * mean(XK);
sd = exp(-r * T) * std(XK);
VC5_IC_inf = VC5 - 1.96 * sd / sqrt(nb_trajectoires);
VC5_IC_sup = VC5 + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'VC5  : '},...
{' P = '},sprintf('%05.3f',VC5),...
{' IC = ['},sprintf('%05.3f',VC5_IC_inf),...
{' , '},sprintf('%05.3f',VC5_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',VC5_IC_sup - VC5_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(VC5_IC_sup - VC5_IC_inf)  *sqrt(tps))...
));


%% ~~~~~~~~~~~~ Méthode de Monte Carlo sans prendre toutes les valeurs  ~~~~~~~~~~~~~~ %%

tic

S = S0 * ones(nb_trajectoires, 1);
X = zeros(nb_trajectoires, 1);
l=1;

a = n/nb;


for i = 1:n
    dW_t = normrnd(0, sqrt(dt), nb_trajectoires, 1);
    S = S .*(1 +r*dt + sigma*sqrt(abs(S)).*dW_t );
    if (i/n) >= (l/nb)
        X = X+S;
        l = l+1;
    end
end
XK = (X)/nb;
XK = max(XK-K,0);

C = exp(-r*T) * mean(XK);
sd = exp(-r * T) * std(XK);
C_IC_inf = C - 1.96 * sd / sqrt(nb_trajectoires);
C_IC_sup = C + 1.96 * sd / sqrt(nb_trajectoires);

tps = toc;

% fonction d'affichage 
disp(strcat({'C    : '},...
{' P = '},sprintf('%05.3f',C),...
{' IC = ['},sprintf('%05.3f',C_IC_inf),...
{' , '},sprintf('%05.3f',C_IC_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',C_IC_sup - C_IC_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(C_IC_sup - C_IC_inf)  *sqrt(tps))...
));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 8; % nombre des IC affichees
    %fprintf('< 8: ICs (normales) >')
    plot([MC ANT VC1 VC2 VC3 VC4 VC5 C], ...
        1:w, 'x')
    %line([K K],[0 5],'Color','green','LineStyle','--')
   
    line([MC_IC_inf MC_IC_sup],[1 1])
    line([A_IC_inf A_IC_sup],[2 2])
    line([VC1_IC_inf VC1_IC_sup],[3 3])
    line([VC2_IC_inf VC2_IC_sup],[4 4])
    line([VC3_IC_inf VC3_IC_sup],[5 5])
    line([VC4_IC_inf VC4_IC_sup],[6 6])
    line([VC5_IC_inf VC5_IC_sup],[7 7])
    line([C_IC_inf C_IC_sup],[8 8])
 
    grid on

    %line([1.362 1.362],[0.5 7.5])

    legend("estimateurs",'IC')
%     L1 = C_IC_sup - C_IC_inf;
%     L2 = C_est - K;
%     limf = [C_IC_inf C_IC_sup] + max(L1,L2)*[-1 1];
%     xlim(limf)
      ylim([0 w+2])
      yticks(1:w)

    yticklabels({'MC','C_a','C_{VC1}','C_{VC2}',...
                 'C_{VC3}','C_{VC4}','C_{VC5}', 'C'})
    %title('Intervalles de confiance (sauf C)')
    
    %P=P+1;
    G="q";
