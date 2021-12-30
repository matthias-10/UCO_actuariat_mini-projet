S = zeros(nt, n+1);

% Methode de Euler

S(:, 1) = S0;
S_anti = S;
VC = S;
for i = 2:(n+1)

    % ~~~~~~~~~~~~~~ Simulation des prix ~~~~~~~~~~~~~~~~ %
    
    dWt = normrnd(zeros(nt,1),sqrt(dt));
    dSi = S(:,i-1).* ...
          ( r*dt + sigma*sqrt(S(:,i-1)).*dWt );
    S(:,i) = S(:,i-1) + dSi;
    S(:,i) = S(:,i) .* (S(:,i) >= 0);    

    % ~~~~~~~~~~~~~ variable antithetique ~~~~~~~~~~~~~~~ %

    dWt_a = -1*dWt;
    %dS_anti = S_anti(:,i-1).* ...
    %         ( r*dt + sigma*sqrt(S_anti(:,i-1)).*dWt_a );
    % sans biais (?):    
    dS_anti = S_anti(:,i-1).* ...
             ( r*dt + sigma*sqrt(S(:,i-1)).*dWt_a );
    S_anti(:,i) = S_anti(:,i-1) + dS_anti;
    S_anti(:,i) = S_anti(:,i) .* (S_anti(:,i) >= 0);

    % ~~~~~~~~~~~~~ variable de controle ~~~~~~~~~~~~~~~~ %

    dWt_vc = normrnd(zeros(nt,1),sqrt(dt));
    dWt_vc = dWt_vc + dWt_vc ...
              .*sign(dWt).*(sign(dWt_vc) - sign(dWt));
    dSi = VC(:,i-1).* ...
          ( r*dt + sigma*sqrt(S(:,i-1)).*dWt_vc );
    VC(:,i)      =    VC(:,i-1) + dSi;
    VC(:,i) = VC(:,i) .* (VC(:,i) >= 0);
end