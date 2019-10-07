function [Thetajj,Sjj] = lcmv_hg_lasso(Svv,Lvj,param)
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
[Thetajj,Sigmajj] = twostep_lasso_caller(Sjj,param);

end