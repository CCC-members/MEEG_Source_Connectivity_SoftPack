function [sigma2xi,Sigmajj,Thetajj] = higgs_maximization(Psixixi,Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
p       = param.p;
axi     = param.axi;
aj      = param.aj;
penalty = param.penalty;
%% Source precision/covariance matrix estimator
if penalty == 0 % naive
    [U,D]                 = eig(Psijj);
    d                     = real(diag(D));
    Thetajj               = U*diag(1./d)*U';
    Thetajj               = (Thetajj + Thetajj')/2;
    Sigmajj               = Psijj;
elseif penalty == 1 % lasso
    [Thetajj,Sigmajj]     = higgs_lasso_caller(Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
elseif penalty == 2 % ridge
    [U,D]                 = eig(Psijj);
    d                     = real(diag(D));
    d_frob                = (sqrt(d.^2 + 4*aj^2) - d)/(2*aj^2);
    Thetajj               = U*diag(d_frob)*U';
    Thetajj               = (Thetajj + Thetajj')/2;
    Sigmajj               = U*diag(1./d_frob)*U';
    Sigmajj               = (Sigmajj + Sigmajj')/2;
end
%%
%% Residual precision/variance estimator
% sigma2xi          = (sum(abs(diag(Psixixi)))/p + axi)*ones(p,1);
sigma2xi          = abs(diag(Psixixi)) + axi;
end