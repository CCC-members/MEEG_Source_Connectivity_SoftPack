function [rth] = higgs_lasso_stabilizer(Thetajj_unb,Thetajj_var,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
rth1      = param.rth1;
rth2      = param.rth2;
m         = param.m;
ntry      = param.ntry;
%% try likelihood evolution for the mask given by all thresholds
rth_grid  = rth1:0.1:rth2;
llh_grid  = zeros(length(rth_grid),ntry);
for rth_count = 1:length(rth_grid)
    %% Zeros mask
    rth                = rth_grid(rth_count);
    mask               = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
    %% Perform ntry iterations
    [llh_grid(rth_count,:)] = try_likelihood(mask,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
end
[rth]     = check_likelihood_trend(llh0,llh_grid,rth_grid,ntry);
end

%%
function [llh] = try_likelihood(mask,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
%% Perform ntry iterations
ntry          = param.ntry;
param.penalty = 2;
llh           = zeros(1,ntry);
for try_iter = 1:ntry
        %% Expectation
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
    %%
    %% Maximization
    [sigma2xi,Sigmajj,Thetajj]            = higgs_maximization(Psixixi,Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
    Thetajj(mask)                         = 0;
    [U,D]                                 = eig(Thetajj);
    d                                     = real(diag(D));
    Sigmajj                               = U*diag(1./d)*U';
    Sigmajj                               = (Sigmajj + Sigmajj')/2;
    %% Compute likelihood
    [llh(try_iter)]                       = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_pst,param);
    llh0                                  = llh(try_iter);
    sigma2xi0                             = sigma2xi;
    Sigmajj0                              = Sigmajj;
end
end