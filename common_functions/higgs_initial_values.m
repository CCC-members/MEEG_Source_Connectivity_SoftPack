function [Svv,Lvj,Ljv,scale,sigma2xi0,Sigmajj0,llh0,param] = higgs_initial_values(Svv,Lvj,param)
%% Parameters
prew            = param.prew;
m               = param.m;
p               = param.p;
q               = param.q;
Iq              = param.Iq;
aj              = param.aj;
penalty         = param.penalty;
%%
%% Lead Field and Data scaling by frobenious norm equivalent to the largest singular values or the 
% average variance of the observations, helps to treat the problen in a standard scale for every case,
% a projected identity matrix Svv = Lvj*Iq*Ljv under this conditions will generate a signal with 
% (sum(abs(diag(Svv)))/p) = 1 if we fix the value of param.axi = 1E-2 the inferior admisible noise 
% will be 1% of the average signal amplitude

scaleLvj        = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj             = Lvj/scaleLvj;
Ljv             = Lvj';
scaleV          = (sum(abs(diag(Svv)))/p);
Svv             = Svv/scaleV;
%%
%% Decide not-prewarming (prew = 0) or prewarming (prew = 1), prew = 1 is recomended for severe 
%  ill-conditioning, whereas the number of observed variables is of the same order of magnitude 
%  that the number of hidden variables use prew = 0, in both cases the scale of Svv is readjusted 
%  by and initial computation of Sjj making the maximum singular values of Sjj and Sigmajj_pst 
%  identique,this will avoid biasing the EM initial step

if prew == 0 % not-prewarming
    %  - if prew = 0 Sjj is precomputed by the first EM Expectation step with unitary initialization 
    %  Sigmajj0 = Iq sigma2xi0 = 1, it rescales Svv and axi by scaleJ the structure of Thetajj and Sigmajj 
    %  is the identity matrix 
    
    Sigmajj0        = Iq;
    sigma2xi0       = ones(p,1);
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
    scaleJ          = (sum(abs(diag(Sjj)))/q)/(sum(abs(diag(Sigmajj_pst)))/q);
    Svv             = Svv/scaleJ;
    param.axi       = param.axi/scaleJ;
    scale           = scaleLvj^2/(scaleV*scaleJ);
    Thetajj         = Iq;
    Sigmajj         = Iq;
elseif prew == 1 % prewarming
    %  - if prew = 1 we use an initial subspace two-step solution based on the crossspectral enet-ssbl 
    %  (equivalent to an univariate higgs) and hermitian graphical model, this provides a multivariate 
    %  initialization to the EM or prewarming, crossspectral enet-ssbl also adjusts the scale by scaleJ
    %  for every penalty a two-step graphical computation is performed to produce the initial structure
    %  of Thetajj and Sigmajj
    
    disp('prewarming')
    nonovgroups     = [];
    counter         = 1;
    for ii = 1:1:q
        nonovgroups{counter} = ii;
        counter              = counter + 1;
    end
    [sigma2j,sigma2j_pst,Tjv,Svv,Lvj,scaleJ,scaleLvj] = cross_nonovgrouped_enet_ssbl({Svv},{Lvj},m,nonovgroups);
    Lvj             = Lvj{1};
    Ljv             = Lvj';
    Tjv             = Tjv{1};
    Tvj             = Tjv';
    Svv             = Svv{1};
    Sjj             = Tjv*Svv*Tvj;
    Sjj             = (Sjj + Sjj')/2;
    if penalty == 0 % naive
        [U,D]                 = eig(Sjj);
        d                     = real(diag(D));
        dfix                  = d + max(d)*1E-4;
        Thetajj               = U*diag(1./dfix)*U';
        Thetajj               = (Thetajj + Thetajj')/2;
        Sigmajj               = U*diag(dfix)*U';
    elseif penalty == 1 % lasso
        [Thetajj,Sigmajj]     = twostep_lasso_caller(Sjj,param);
    elseif penalty == 2 % ridge
        [U,D]                 = eig(Sjj);
        d                     = real(diag(D));
        dfix                  = d + max(d)*1E-4;
        d_frob                = (sqrt(dfix.^2 + 4*aj^2) - dfix)/(2*aj^2);
        Thetajj               = U*diag(d_frob)*U';
        Thetajj               = (Thetajj + Thetajj')/2;
        Sigmajj               = U*diag(1./d_frob)*U';
        Sigmajj               = (Sigmajj + Sigmajj')/2;
    end
    scaleLvj        = scaleLvj{1};
    scaleJ          = scaleJ{1};
    param.axi       = param.axi/scaleJ;
    scale           = scaleLvj^2/(scaleV*scaleJ);
end
%% General warm-start by cross-validating diagonal initial values
%  search for the optimal initial values sigma2xi0 Sigmajj0 by cross-validating sigma2xi and sigma2j*Sigmajj,
%  Sigmajj is defined above in any of the cases prew = 0 (Sigmajj = Iq) or prew = 1 (Sigmajj is determined by 
%  the twostep_lasso_caller), the optimal criteria is due to the hyperparameters posterior distribution 

disp('warm-start higgs')
[sigma2xi0,Sigmajj0,llh0] = higgs_warm_start(Svv,Lvj,Ljv,Thetajj,Sigmajj,param);
end
