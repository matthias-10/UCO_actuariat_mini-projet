

r = 0.02;
sigma = 0.02;
S0 = 40;

T = 100;
t0=0;
n=2^8;
dt = ((T-t0)/n);
t = t0:dt:T;

t_mm = n;
nt = 10; % nombre des trajectoires



tic

% simulation en forme matrice:

S_sim = zeros(length(t),nt);
for i = 1:nt
    S_sim(:,i) = brownmo(S0, r, sigma ,t0, T, n); %brownmo est definie en bas
end

duree= toc; 
% Matthias HPPavilion: ~ 0.22, le script R: ~ 0.09 secs; nt=30 n=2^8
fprintf('%d trajectoires plotees\n', nt);
fprintf('Fini en %0.5g\n', duree);


% plot
tiledlayout(3,1)

x_ax = t; % x-axe

axis([0 T 0.9*min(min(S_sim)) 1.1*max(max(S_sim))]) %x-axe limits
%axis([0 T S0-sigma*n S0*(1+r)^(T-t0)+sigma*n]); 

nexttile
hold on
% pour comparison, si j'epargne pour le taux r:
plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
plot(x_ax,S_sim);


% calculer X_T
 % cumsum au lieu de l'integrale

 %X = [0 1 2; 3 4 5];cumsum(X,1)

% X_T = cumsum(S_sim,1);
% X_T = X_T(n+1,:);
% X_T = X_T/T;

% calculer X_t
hold off
nexttile
hold on

plot([t0 T], [S0 S0*(1+r)^(T-t0)], "--k"); % obligation
plot(x_ax(1:(t_mm+1)),MoyMob(S_sim, t_mm), ':');
plot(x_ax,S_sim);

%hold off

%MoyMob(S_sim, t_mm)
%size(MoyMob(S_sim, t_mm))

% calculer K = S_0 * exp(rT)

K = S0*exp(r*T);

% calculer le prix de l'option C_t = max(0, X_t-S_t)
Xn = MoyMob(S_sim, n);

hold off
nexttile
hold on

plot(x_ax, Xn-S_sim);
plot([t0 T], [0 0], ":k"); % zero

hold off
nexttile

C = ( (Xn - S_sim) >= 0 ) .* (Xn - S_sim);


plot(x_ax, C); % Il y a trop de valeurs >0 ? r?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% fonctions %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = brownmo(X0, mu, sigma, t0, T, n) %x0 

  delta = (T-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((T-t0)/n):T;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp( (mu-(sigma^2)/2) * (tseq-t0) + sigma*W );

end


function X_t = MoyMob(M, t_m)
    X_t = cumsum(M,1);
    t_m2= t_m+1;
    X_t = X_t(1:t_m2,:);
    for i = 1:t_m2
        X_t(i,:) = X_t(i,:)/i; %vectoriel?
    end
end








