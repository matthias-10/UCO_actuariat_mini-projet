% n trajectoires

hold on %pour faire plusieurs plots

for i = 1:nt
    plot(brownmo(P0, mu, sigma ,t0, T, n))
end

hold off


disp("fini")


function S = brownmo(X0, mu, sigma, t0, t, n) %x0 

  delta = (t-t0)/n;
  W = zeros(1,n+1);
  tseq = t0:((t-t0)/n):t;
  for i = 2:(n+1)
    W(i) = W(i-1)+normrnd(0,1)*sqrt(delta);
  end
  S = X0 * exp((mu-(sigma^2)/2)*(tseq-t0)+sigma*W);
end
