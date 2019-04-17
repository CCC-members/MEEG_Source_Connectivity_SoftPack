function [Thetajj,Sjj] = lcmv_hggm(Svv,Lvj,param)
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
gamma         = param.gamma;
%%
[Tjv,T1jv,Wout]   = mkfilt_lcmv(Lvj,Svv,gamma);
Tjv               = Tjv';
%%
Sjj               = Tjv*Svv*Tjv';
Sjj               = (Sjj + Sjj')/2;
[Thetajj,Sigmajj] = hggm_solve(Sjj,m,aj,Ajj,Iq,rth,maxiter_inner);

end