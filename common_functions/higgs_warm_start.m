function [sigma2xi0,Sigmajj0,llh0] = higgs_warm_start(Svv,Lvj,Ljv,Thetajj,Sigmajj,param)
%% Warm-start by cross-validating diagonal initial values
p               = param.p;
grid_sigma2xi   = (-5:0.25:5);
grid_sigma2j    = (-5:0.25:5);
llh          = zeros(length(grid_sigma2xi),length(grid_sigma2j));
for cont_sigma2xi = 1:length(grid_sigma2xi)
    for cont_sigma2j = 1:length(grid_sigma2j)
        sigma2xi0        = 10^(grid_sigma2xi(cont_sigma2xi));
        sigma2xi0        = sigma2xi0*ones(p,1);
        sigma2j0         = 10^(grid_sigma2j(cont_sigma2j));
        Sigmajj0         = sigma2j0*Sigmajj;
        Thetajj0         = (1/sigma2j0)*Thetajj;
        % Expectation
        [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
        % Compute likelihood
        llh(cont_sigma2xi,cont_sigma2j)       = higgs_likelihood(Sxixi,sigma2xi0,Sjj,Thetajj0,Sigmajj_pst,param);
    end
end
%% optimal value
[val]           = max(llh(:));
[id_xi,id_j]    = find(llh == val);
sigma2xi        = 10^(grid_sigma2xi(id_xi));
sigma2xi        = sigma2xi*ones(p,1);
sigma2j         = 10^(grid_sigma2j(id_j));
Sigmajj         = sigma2j*Sigmajj;
Thetajj         = (1/sigma2j)*Thetajj;
%% Recompute likelihood
% Expectation
[Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi,Sigmajj,param);
% Compute likelihood
llh0      = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_pst,param);
sigma2xi0 = sigma2xi;
Sigmajj0  = Sigmajj;
end