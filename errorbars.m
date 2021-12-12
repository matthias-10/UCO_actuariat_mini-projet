randn('state', 100);
S=10;   E=9;   sigma=0.1;   r=0.06;   T=1;   Dt=1e-3;   N=T/Dt; 
M = [2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14,2^15,2^16,2^17];
hold on;
for k=1:numel(M)
    %No need of loop here. Generate all random values in one go 
    Sfinal = S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*randn(M(k),1));
    V = exp(-r*T)*max(Sfinal-E,0);
    aM = mean(V);   bM = std(V);

    plot(M(k),aM,'x');
    errorbar(M(k), aM, 1.96*bM/sqrt(M(k)));    
end
chvar = repmat(ch08(10,9,0.06,0.1,1),1,numel(M));  %<----Notice this
plot(M, chvar,'--.k'); 

%Other figure cosmetics
%These commands shouldn't be inside the loop to avoid unnecessary computations


plot([1 2],[X_ab_mu X_mu],'x');
hold on
errorbar([1 2], [X_ab_mu X_mu] , [2 2]);  
hold off
titre = ['Monte Carlo estimateurs de X avec '...
          num2str(nt) ' trajectoires simulees' ];
title(titre);
legend("test2","","test3","");
xlabel(''); % x-axis label
ylabel('Option value approximation'); % y-axis label

boxplot([1 2 ;  5 6; 7 8], 'whisker', 0, ...
        'widths', 0.001,'orientation','horizontal', ...
        'medianstyle', 'target', ...
        'labels' , ['test1'; 'test2'])

plot([1.2 1.1 1.6 1.5],[1 2 3 4], 'x')
xlim([0 3])
ylim([0 5])
line([1 2],[2 2])

function [C, Cdelta, P, Pdelta] = ch08(S,E,r,sigma,tau)
% Program for Chapter 8
% This is a MATLAB function
%
% Input arguments: S = asset price at time t
%                  E = Exercise price
%                  r = interest rate
%                  sigma = volatility
%                  tau = time to expiry (T-t) 
%
% Output arguments: C = call value, Cdelta = delta value of call 
%                   P = Put value, Pdelta = delta value of put
%
%   function [C, Cdelta, P, Pdelta] = ch08(S,E,r,sigma,tau)

if tau > 0
   d1 = (log(S/E) + (r + 0.5*sigma^2)*(tau))/(sigma*sqrt(tau));
   d2 = d1 - sigma*sqrt(tau);
   N1 = 0.5*(1+erf(d1/sqrt(2)));
   N2 = 0.5*(1+erf(d2/sqrt(2)));
   C = S*N1-E*exp(-r*(tau))*N2;
   Cdelta = N1;
   P = C + E*exp(-r*tau) - S;
   Pdelta = Cdelta - 1;
else 
   C = max(S-E,0);
   Cdelta =  0.5*(sign(S-E) + 1);
   P = max(E-S,0);
   Pdelta = Cdelta - 1;
end

end