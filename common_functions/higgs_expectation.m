function [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi,Sigmajj,param)
Ip                = param.Ip;
%% Source Posterior Covariance (SPC)
SigmajjLjv        = Sigmajj*Ljv; % SEC*SDTF'
LvjSigmajj        = SigmajjLjv'; % Tranconjugated SEC*SDTF'
Sigmajj_pst       = Sigmajj - SigmajjLjv*(Ip/(Lvj*SigmajjLjv+diag(sigma2xi)))*LvjSigmajj; % SPC by Woodbury formula
Sigmajj_pst       = (Sigmajj_pst + Sigmajj_pst')/2;
%%
%% Data to Source Transfer Function (DSTF)
Tjv               = Sigmajj_pst*Ljv*diag(1./sigma2xi); % DSTF
Tvj               = Tjv'; % Tranconjugated DSTF
%%
%% Source Empirical Covariance (ESEC)
Sjj               = Tjv*Svv*Tvj; % Source Empirical Covariance (SEC)
Sjj               = (Sjj + Sjj')/2;
%%
%% Effective Source Empirical Covariance (ESEC)
Psijj             = Sigmajj_pst + Sjj; % ESEC
%%
%% Data to Residuals Transfer Function (DRTF)
Txiv              = (Ip - Lvj*Tjv); % DRTF
Tvxi              = Txiv'; % Tranconjugated DRTF
%%
%% Residuals Posterior Covariance (SPC)
Sigmaxixi_pst     =  Lvj*Sigmajj_pst*Ljv;
Sigmaxixi_pst     = (Sigmaxixi_pst + Sigmaxixi_pst')/2;
%%
%% Residuals Empirical Covariance (REC)
Sxixi             = Txiv*Svv*Tvxi; % Residual Empirical Covariance (REC)
Sxixi             = (Sxixi + Sxixi')/2;
%% Effective Residual Empirical Covariance (EREC)
Psixixi           = Sigmaxixi_pst + Sxixi; % EREC
end