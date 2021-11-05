tic

mu = 2;
sigma = 0.2;
P0 = 40;

T = 1/12;
t0=0;
n=2^8;
dt = ((T-t0)/n);
t = t0:dt:T;

nt = 30;


% n trajectoires
plot(P0); %pour effacer le plot
hold on %pour faire plusieurs plots
x.ax = t; % x-axe

axis([0 T P0+mu*(-2) P0+mu*8]); %x-axe limits

for i = 1:nt
    plot(x.ax, brownmo(P0, mu, sigma ,t0, T, n)) %brownmo est definie en bas
end

hold off

% tic toc, comparer avec r

duree= toc;
disp(sprintf('%d trajectoires', nt));
disp(sprintf('Fini en %0.5g', duree));



function S = brownmo(X0, mu, sigma, t0, t, n) %x0 

  delta = (t-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((t-t0)/n):t;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp((mu-(sigma^2)/2)*(tseq-t0)+sigma*W);

end
