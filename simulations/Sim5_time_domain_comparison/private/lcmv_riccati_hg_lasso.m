function [Thetajj,Sjj,Tjv] = lcmv_riccati_hg_lasso(Svv,Lvj,param)
p             = size(Lvj,1);
scaleLvj      = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj           = Lvj/scaleLvj;
scaleV        = (sum(abs(diag(Svv)))/p);
Svv           = Svv/scaleV;
gamma         = param.gamma;
%%
[Tjv,T1jv,Wout]   = mkfilt_lcmv(Lvj,Svv,gamma);
Tjv               = Tjv';
%%
Sjj               = Tjv*Svv*Tjv';
Sjj               = (Sjj + Sjj')/2;
Sigmajj                       = riccati_leakage_correction(Lvj,Tjv,Sjj,param);
[Thetajj,Sigmajj] = twostep_lasso_caller(Sigmajj,param);

end