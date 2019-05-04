function [Thetajj,Sjj,gamma_grid,gamma,gcv] = eloreta_hggm(Svv,Lvj,param)
p             = size(Lvj,1);
q             = size(Lvj,2);
Ip            = eye(p);
Iq            = eye(q);
maxiter_inner = param.maxiter_inner;
m             = param.m;
aj            = sqrt(log(q)/m);
rth           = param.rth;
Ajj_diag      = 0;
Ajj_ndiag     = 1;
Ajj           = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
gamma1        = param.gamma1;
gamma2        = param.gamma2;
delta_gamma   = param.delta_gamma;
gamma_grid    = gamma1:delta_gamma:gamma2;
gcv           = zeros(length(gamma_grid),1);
count         = 1;

%%

for gamma = gamma_grid
   [Tjv,Wout]       = mkfilt_eloreta(Lvj,gamma);
    Tjv              = Tjv';
    Txiv             = Ip - Lvj*Tjv;
    gcv(count)       = (1/p)*sum(abs(diag(Txiv*Svv*Txiv')))/((1/p)*sum(abs(diag(Txiv))))^2;
    count             = count + 1;
end
%%
[gcv_opt,idx_gamma]       = min(gcv);
gamma                     = gamma_grid(idx_gamma);
[Tjv,Wout]                = mkfilt_eloreta(Lvj,gamma);
Tjv                       = Tjv';
Sjj                       = Tjv*Svv*Tjv';
Sjj                       = (Sjj + Sjj')/2;
[Thetajj,Sigmajj]         = hggm_solve(Sjj,m,aj,Ajj,Iq,rth,maxiter_inner);

end