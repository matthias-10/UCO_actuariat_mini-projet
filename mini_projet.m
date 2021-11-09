%% Parametres

r = 0.05;
sigma = 0.5;
T = 1;
n=2^8;
nt = 10; % nombre de trajectoires
N = 5; % pour le cas discret, verifier que N << n
S0 = 40;
t0=0;

% calculer K

K = S0*exp(r*T);
%K = int(obligation,t0,T)/(T-t0);

%% Simulation
% simulation en forme matrice:

dt = ((T-t0)/n);
t = t0:dt:T;

%simulation de nt trajectoires
S = zeros(length(t),nt);
for i = 1:nt
    S(:,i) = brownmo(S0, r, sigma ,t0, T, n); %brownmo est definie en bas
end
plot(S)

%calcul de X_t
vecX_t = zeros(1,nt);
for i = 1:nt
    for j = (1:n)
        vecX_t(i) = vecX_t(i) + (S(j,i)+S(j+1,i))/2;
    end
end
vecX_t = vecX_t/n;
X_t = mean(vecX_t(:,1));
vecC_inf = vecX_t-K;
for i = 1:nt
    vecC_inf(i) = max(vecC_inf(i),0);
end
C_inf = mean(vecC_inf);
%C_inf * exp(-rT) est une martingale donc E[exp(-rT)*C_inf]= C_inf(S_0)
C_inf_0 = exp(-r*T)*C_inf;

%calcul de X_t_prim
vecX_t_prim = S(1,:);
for i = 2:(n+1)
    vecX_t_prim = vecX_t_prim + S(i,:);
end
vecX_t_prim = vecX_t_prim/(n+1);
X_t_prim = mean(vecX_t_prim);
vecC_N = vecX_t_prim-K;
for i = 1:nt
    vecC_N(i) = max(vecC_N(i),0);
end
C_N = mean(vecC_N);
%C_N * exp(-rT) est une martingale donc E[exp(-rT)*C_N]= C_N(S_0)
C_N_0 = exp(-r*T)*C_N;




function S = brownmo(X0, mu, sigma, t0, T, n) %x0 

  delta = (T-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((T-t0)/n):T;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp( (mu-(sigma^2)/2) * (tseq-t0) + sigma*W );

end

