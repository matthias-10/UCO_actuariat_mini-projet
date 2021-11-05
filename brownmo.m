tic

r = 2;
sigma = 0.2;
P0 = 40;

T = 1/12;
t0=0;
n=2^5;
dt = ((T-t0)/n);
t = t0:dt:T;

nt = 1;


% n trajectoires
% plot(P0); %pour effacer le plot
% hold on %pour faire plusieurs plots
% x_ax = t; % x-axe
% 
% axis([0 T P0+r*(-2) P0+r*8]); %x-axe limits
% 
% for i = 1:nt
%     plot(x_ax, brownmo(P0, r, sigma ,t0, T, n)) %brownmo est definie en bas
% end




duree= toc; % Matthias HPPavilion: ~ 0.22, le script R: ~ 0.09 secs; nt=30 n=2^8
disp(sprintf('%d trajectoires plote', nt));
disp(sprintf('Fini en %0.5g', duree));

% nouvelle simulation en forme matrice:

S_sim = zeros(length(t),nt);
for i = 1:nt
    S_sim(:,i) = brownmo(P0, r, sigma ,t0, T, n); %brownmo est definie en bas
end

% plot
plot(P0); %pour effacer le plot
hold on %pour faire plusieurs plots
x_ax = t; % x-axe

axis([0 T P0+r*(-2) P0+r*8]); %x-axe limits

for i = 1:nt
    plot(x_ax,S_sim(:,i));
    %plot(x_ax, brownmo(P0, r, sigma ,t0, T, n)) 
end

% calculer X_T
 % au lieu de l'integrale cumsum

 %X = [0 1 2; 3 4 5];cumsum(X,1)

X_T = cumsum(S_sim,1);
X_T = X_T(n+1,:);
X_T = X_T/T;

% calculer X_t


plot(x_ax(1:(20+1)),MoyMob(S_sim, 20), ':');
hold off
% calculer le prix des options
MoyMob(S_sim, 20)
size(MoyMob(S_sim, 20))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% fonctions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = brownmo(X0, mu, sigma, t0, t, n) %x0 

  delta = (t-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((t-t0)/n):t;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp((mu-(sigma^2)/2)*(tseq-t0)+sigma*W);

end


function X_t = MoyMob(M, t_m)
    X_t = cumsum(M,1);
    t_m2= t_m+1;
    X_t = X_t(1:t_m2,:);
    for i = 1:t_m2
        X_t(i,:) = X_t(i,:)/i;
    end
end
