function [Thetajj,Sigmajj,llh_inner] = hggm_solve(Psijj,m,aj,Ajj,Iq,rth,maxiter_inner)
%% Precompute hggm_lasso_ssbl for gamma = ones(q) equivalent to tikhonov regularization (Frobenious norm) solution
[U,D,V]            = svd(Psijj);
d                  = abs(diag(D));
d_frob             = (sqrt(d.^2 + 4*aj^2) - d)/(2*aj^2);
Thetajj_ridge      = U*diag(d_frob)*V';
Thetajj_ridge      = (Thetajj_ridge + Thetajj_ridge')/2;
%% Compute hggm_lasso_ssbl for scaled empirical covariance
dmin               = min(d);
if dmin < 0
    Psijj          = Psijj + abs(dmin)*Iq + Iq;
else
    Psijj          = Psijj - abs(dmin)*Iq + Iq;
end
scale_Psijj        = sqrt(mean(abs(diag(Psijj*Psijj'))));
Psijj              = Psijj/scale_Psijj;
[Thetajj_lasso,llh_inner] = hggm_lasso_ssbl(Psijj,m,Ajj,aj,maxiter_inner);
%% Compute Rayleigh threshold
Thetajj_unb        = 2*Thetajj_lasso - Thetajj_lasso*Psijj*Thetajj_lasso;
Thetajj_unb        = (Thetajj_unb + Thetajj_unb')/2;
Thetajj_var        = sqrt(abs(diag(Thetajj_lasso))*abs(diag(Thetajj_lasso))' + abs(Thetajj_lasso).^2);
Thetajj_var        = (Thetajj_var + Thetajj_var')/2;
%% Apply Rayleigh threshold to precomputed ridge solution
Thetajj_ridge(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
%%
Thetajj            = Thetajj_ridge;
Sigmajj            = Iq/Thetajj;
Sigmajj            = (Sigmajj + Sigmajj')/2;
end