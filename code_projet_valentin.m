
S0 = 40;                % Prix initial du sous jacent
K = 41;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode

nb_trajectoires = 100;

dt = (T-t0)/n;
t = t0:dt:T;


tic

%t_d = datetime('now');
S = S0 * ones(nb_trajectoires,1);
X = S/2;

% Simulation pas a pas
for i = 1:n
    dW_t = normrnd(0,sqrt(dt), nb_trajectoires,1);
    S = S .* (1 + r*dt + sigma*sqrt(S).*dW_t );
    X = X+S;
end

X = X - S/2;

XT = X / n;
XK = max(XT-K,0);
P = exp(-r * T ) * mean(XK);
