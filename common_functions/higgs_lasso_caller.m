function [Thetajj,Sigmajj,llh_inner] = higgs_lasso_caller(Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
%% Effectutates the higgs-lasso maximization step
%  attempts to solve two ambigities in the higgs-lasso maximization 1st- the hg-lasso sparsity biases the results
%  in local and global computations of Thetajj,

m              = param.m;
aj             = param.aj;
Ajj            = param.Ajj;
nu             = param.nu;
rth1           = param.rth1;
rth2           = param.rth2;
maxiter_inner  = param.maxiter_inner;
ntry           = param.ntry;
%% Compute hg_lasso_lqa
%  it first corrects the eigenvalues of the effective empirical covariance Psijj and normalize it by the infimun
%  norm Psijjfix, the scaling guarantees more stable computations of the hg-lasso-lqa, continues with performing
%  hg-lasso-lqa Thetajj_lasso and the unbiased statistics Thetajj_unb and Thetajj_var

[U,D]          = eig(Psijj);
d              = real(diag(D));
dfix           = d + max(d)*1E-4;
Psijjfix       = U*diag(dfix)*U';
Psijjfix       = Psijjfix/sqrt(sum(abs(diag(Psijjfix*Psijjfix')))/length(Psijjfix));
[Thetajj_lasso,llh_inner] = hg_lasso_lqa1(Psijjfix,m,Ajj,aj,nu,maxiter_inner);
Thetajj_unb    = 2*Thetajj_lasso - Thetajj_lasso*Psijjfix*Thetajj_lasso;
Thetajj_unb    = (Thetajj_unb + Thetajj_unb')/2;
Thetajj_var    = sqrt(abs(diag(Thetajj_lasso))*abs(diag(Thetajj_lasso))' + abs(Thetajj_lasso).^2);
Thetajj_var    = (Thetajj_var + Thetajj_var')/2;
%% Precompute hg_lasso_lqa equivalent to ridge solution
%  this is the closest solution to the hermitian graphical lasso with parameter aj and equivalent to the first
%  iteration of the lqa with gamma == 1, if for a given iteration the lasso model biases the computations and
%  the global hyperparameters posterior likelihood does not increases in ntry iterations higgs-lasso still goes
%  on with higgs-ridge prewarming, this approximation is a result of the empirical study of higgs-ridge with
%  parameter aj^2 (aj same as lasso aj = sqrt(log(q)/m)) in terms of stability and similarity to the unbiased
%  statistics given by Jankova conditions,

d_frob         = (sqrt(d.^2 + 4*aj^2) - d)/(2*aj^2);
Thetajj        = U*diag(d_frob)*U';
Thetajj        = (Thetajj + Thetajj')/2;
%% try likelihood evolution for the mask given by all thresholds
%  evaluate performance of Rayleigh thresholds in terms of the global
%  hyperparameters posterior likelihood, the overlapping of the distributions for null
%  hypothesis values (Rayleigh) and the alternative introduces ambiguity in the decision for a given Rayleigh threshold
%  we move rth between the value corresponding to the peak of the Rayleigh
%  distribution (rth1)

if ntry == 0
    %% check only partial likelihood
    rth_grid         = rth1:0.1:rth2;
    llhjj_grid       = zeros(length(rth_grid),1);
    for rth_count = 1:length(rth_grid)
        %% Zeros mask
        rth                   = rth_grid(rth_count);
        mask                  = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
        Thetajj_mask          = Thetajj;
        Thetajj_mask(mask)    = 0;
        [U,D]                 = eig(Thetajj_mask);
        d_mask                = abs(diag(D));
        llhjj_grid(rth_count) = sum(log(abs(d_mask))) - sum(abs(diag(Thetajj*Psijj))) - aj*sum(abs(Thetajj(:)));
    end
    [val,id_rth]     = max(llhjj_grid);
    id_rth           = min(id_rth);
    rth              = rth_grid(id_rth);
else
    %% check global likelihood
    [rth]            = higgs_lasso_stabilizer(Thetajj_unb,Thetajj_var,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
end
Thetajj(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
[U,D]          = eig(Thetajj);
d              = real(diag(D));
Sigmajj        = U*diag(1./d)*U';
Sigmajj        = (Sigmajj + Sigmajj')/2;
end
