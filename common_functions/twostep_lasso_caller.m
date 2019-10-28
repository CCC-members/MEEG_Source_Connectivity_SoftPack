function [Thetajj,Sigmajj] = twostep_lasso_caller(Psijj,param)
m              = param.m;
aj             = param.aj;
Ajj            = param.Ajj;
nu             = param.nu;
rth1           = param.rth1;
rth2           = param.rth2;
maxiter_inner  = param.maxiter_inner;
%% Compute hg_lasso_lqa
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
%% Precompute hg_lasso_lqa equivalent to Frobenious solution
d_frob         = (sqrt(d.^2 + 4*aj^2) - d)/(2*aj^2);
Thetajj        = U*diag(d_frob)*U';
Thetajj        = (Thetajj + Thetajj')/2;
%% try likelihood evolution for the mask given by all thresholds
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
Thetajj(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
[U,D]            = eig(Thetajj);
d                = real(diag(D));
Sigmajj          = U*diag(1./d)*U';
Sigmajj          = (Sigmajj + Sigmajj')/2;
end
