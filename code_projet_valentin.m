
S0 = 40;                % Prix initial du sous jacent
K = 41;                 % Prix d'exercice de l'option

r = 0.05;               % Taux d'interet sous risque neutre
sigma = 0.01;           % Variance partie fixe

t0 = 0;                 % Debut de la periode
n = 2^9;                % Nombre de intervalles
T = 1;                  % Fin de la periode

nb_trajectoires = 1000;

dt = (T-t0)/n;
t = t0:dt:T;

S = zeros(15,n+1);
S(:,1) = S0;

% Simulation pas a pas
for i = 2:(n+1)
    dW_t = normrnd(zeros(15,1),sqrt(dt));
    dSi = S(:,i-1).*( r*dt + sigma*sqrt(S(:,i-1)).*dW_t );
    S(:,i) = S(:,i-1) + dSi;
end

 plot([t0 T],[K K], ':k', 'LineWidth',2)
        hold on
        plot(t, S)
        plot([t0 T],[K K], ':k', 'LineWidth',2)
        hold off



%% ~~~~~~~~~~~~ Méthode de Monte Carlo ~~~~~~~~~~~~~~ %%

tic;
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
sd = exp(-r * T) * std(XK);
i_inf = P - 1.96 * sd;
i_sup = P +1.96 * sd;

tps = toc;

% fonction d'affichage 
disp(strcat({'MC : '},...
{' P = '},sprintf('%05.3f',P),...
{' IC = ['},sprintf('%05.3f',i_inf),...
{' , '},sprintf('%05.3f',i_sup),...
{'] '},...
{' largeur = '},sprintf('%05.3f',i_sup - i_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(i_sup - i_inf)  *sqrt(tps))...
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
    S = S .* (1 + r*dt + sigma*sqrt(S).*dW_t );
    Sa = Sa .* (1 + r*dt - sigma*sqrt(Sa).*dW_t );
    X = X+S;
    Xa = Xa + Sa;
end

X = X - S/2;
Xa = Xa - Sa/2;

XT = X / n;
XTa = Xa/n;

XK = max(XT-K,0);
XKa = max(XTa - K,0);
P = exp(-r * T ) * mean(XK/2 + XK/2);
sd = exp(-r * T) * std(XK);
i_inf = P - 1.96 * sd;
i_sup = P +1.96 * sd;

tps = toc;

% fonction d'affichage 
disp(strcat({'ANT : '},...
{' P = '},sprintf('%05.3f',P),...
{' IC = ['},sprintf('%05.3f',i_inf),...
{' , '},sprintf('%05.3f',i_sup),... 
{'] '},...
{' largeur = '},sprintf('%05.3f',i_sup - i_inf),...
{' t = '},sprintf('%05.3f',tps),...
{' eff = '},sprintf('%05.3f',(i_sup - i_inf)  *sqrt(tps))...
)); 
