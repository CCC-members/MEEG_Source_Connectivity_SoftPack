function [Thetajj,Sjj,llh,jj_on,xixi_on] = h_hggm(Svv,Lvj,param)
% Hidden Hermitian Gaussian Grpahical Model (H_HGGM) solver. Computes the Source Empirical Covariance (Sjj) and Source Partial Correlations (Thetajj) by two sequential steps.
% First: Unhides (Expectation) the Type II Likelihood approximated representation which shapes a pair of Hermitian Gaussian Graphical Model (HGGM), one of the state equation
% (with empirical covariance Psijj) and another of the observation equation residuals (with empirical covariance Psixixi). The Hyperparameters are computed by maximum posterior
% analysis (Maximization) regularized with priors. The states (source) HGGM is estimated with hermitian graphical lasso solver controled by Rayleigh threshold
% and the residuals HGGM with exponential prior of the noise presicion controled by nuisance inferior limit (scale parameter).
%
%
% inputs:
%    Svv                 : M/EEG time series cross-spectra
%    Lvj                 : M/EEG Lead Field
%    param.maxiter_outer : number of iterations of the h-hggm unhide and solve (em) loop
%    param.maxiter_inner : number of iterations of the source_hggm_solver loop
%    param.m             : Sample number of MEEG cross-spectra
%    param.penalty       : Penalty type for source partial correlations 0 (penalty free) 1 (hermitian graphical lasso) 2 (hermitian graphical ridge)
%    param.rth           : Rayleigh limit source partial correlations
%    param.axi           : Noise inferior limit
%    param.Axixi         : Noise partial correlations
%
% outputs:
%    Thetajj             : source partial correlations
%    Sjj                 : source empirical covariance
%    llh_outer           : likelihood of the h-hggm unhide and solve (em) loop
%    llh_inner           : likelihood of the source_hggm_solver loop
%    jj_on               : number of iterations in which the source HGGM maximizes the likelihood
%
% Pedro Valdes-Sosa, March 2019
% Deirel Paz Linares, March 2019
% Eduardo Gonzalez-Moreira, March 2019
%**************************************************************************
%% step 1 Initialization EM algorithm
maxiter_outer         = param.maxiter_outer;
maxiter_inner         = param.maxiter_inner;
m                     = param.m;
penalty               = param.penalty;
rth                   = param.rth;
axi                   = param.axi;
sigma2xi              = param.sigma2xi;
Axixi                 = param.Axixi;
%%
[sigma2xi,theta2xi,p,Ip,Axixi_inv,Lvj,Ljv,Svv,Sigmajj,Thetajj,Iq,Ajj,aj,llh] = h_hggm_initial_values(Svv,Lvj,m,maxiter_outer,sigma2xi,Axixi);
%% Outer loop EM algorithm
jj_on   = 0;
xixi_on = 0;
for k_outer = 1:maxiter_outer
    disp(['h-hggm unhide and solve (em) loop # ',num2str(k_outer)])
    %% Source Posterior Covariance (SPC)
    SigmajjLjv        = Sigmajj*Ljv; % SEC*SDTF'
    LvjSigmajj        = SigmajjLjv'; % Tranconjugated SEC*SDTF'
    Sigmajj_pst       = Sigmajj - SigmajjLjv*(Ip/(Lvj*SigmajjLjv+sigma2xi*Axixi_inv))*LvjSigmajj; % SPC by Woodbury formula
    Sigmajj_pst       = (Sigmajj_pst + Sigmajj_pst')/2;
    %%
    %% Data to Source Transfer Function (DSTF)
    Tjv               = Sigmajj_pst*Ljv*Axixi*theta2xi; % DSTF
    Tvj               = Tjv'; % Tranconjugated DSTF
    %%
    %% Source Empirical Covariance (ESEC)
    Sjj               = Tjv*Svv*Tvj; % Source Empirical Covariance (SEC)
    Sjj               = (Sjj + Sjj')/2;
    %%
    %% Effective Source Empirical Covariance (ESEC)
    Psijj             = Sigmajj_pst + Sjj; % ESEC
    %%
    %% Source precision/covariance matrix estimator
    if penalty == 0 % naive
        if k_outer > 1
            llhjj_nupd        = sum(log(abs(d))) - sum(abs(diag(Thetajj*Psijj))); % llh without updating variable
        end
        Sigmajj_upd           = Psijj;
        Thetajj_upd           = Iq/Sigmajj_upd;
        Thetajj_upd           = (Thetajj_upd + Thetajj_upd')/2;
        [U,D,V]               = svd(Thetajj_upd);
        d                     = abs(diag(D));
        llhjj_upd             = sum(log(abs(d))) - sum(abs(diag(Thetajj_upd*Psijj))); % llh updating variable
    elseif penalty == 1 % lasso
        if k_outer > 1
            llhjj_nupd        = sum(log(abs(d))) - sum(abs(diag(Thetajj*Psijj))) - aj*sum(abs(Ajj(:).*Thetajj(:))); % llh without updating variable
        end
        [Thetajj_upd,Sigmajj_upd] = hggm_solve(Psijj,m,aj,Ajj,Iq,rth,maxiter_inner);
        [U,D,V]               = svd(Thetajj_upd);
        d                     = abs(diag(D));
        llhjj_upd             = sum(log(abs(d))) - sum(abs(diag(Thetajj_upd*Psijj))) - aj*sum(abs(Ajj(:).*Thetajj_upd(:))); % llh updating variable
    elseif penalty == 2 % ridge
        if k_outer > 1
            llhjj_nupd        = sum(log(abs(d))) - sum(abs(diag(Thetajj*Psijj))) - (aj^2/2)*sum(abs(diag(Thetajj*Thetajj'))); % llh without updating variable
        end
        [U,D,V]               = svd(Psijj);
        d                     = abs(diag(D));
        d                     = (sqrt(d.^2 + 4*aj^2) - d)/(2*aj^2);
        Thetajj_upd           = U*diag(d)*V';
        Thetajj_upd           = (Thetajj_upd + Thetajj_upd')/2;
        Sigmajj_upd           = V*diag(1./d)*U';
        Sigmajj_upd           = (Sigmajj_upd + Sigmajj_upd')/2;
        llhjj_upd             = sum(log(abs(d))) - sum(abs(diag(Thetajj_upd*Psijj))) - (aj^2/2)*sum(abs(diag(Thetajj_upd*Thetajj_upd'))); % llh updating variable
    end
    %%
    %% Data to Residuals Transfer Function (DRTF)
    Txiv              = (Ip - Lvj*Tjv); % DRTF
    Tvxi              = Txiv'; % Tranconjugated DRTF
    %%
    %% Effective Residual Empirical Covariance (EREC)
    Psixixi           = Txiv*Svv*Tvxi + Lvj*Sigmajj_pst*Ljv; % EREC
    Psixixi           = (Psixixi + Psixixi')/2;
    %%
    %% Residual precision/variance estimator
    llhxixi_nupd      = p*log(theta2xi) - theta2xi*sum(abs(diag(Axixi*Psixixi))) - p*axi*theta2xi; % llh without updating variable
    sigma2xi_upd      = sum(abs(diag(Axixi*Psixixi)))/p + axi;
    theta2xi_upd      = 1/sigma2xi_upd;
    llhxixi_upd       = p*log(theta2xi_upd) - theta2xi_upd*sum(abs(diag(Axixi*Psixixi))) - p*axi*theta2xi_upd; % llh updating variable
    %%
    %% Check llh trend and update source precision/covariance
    if k_outer == 1
        jj_on         = jj_on + 1;
        Thetajj       = Thetajj_upd;
        Sigmajj       = Sigmajj_upd;
        llhjj         = llhjj_upd;
    elseif llhjj_upd + llhxixi_nupd > llhjj + llhxixi
        jj_on         = jj_on + 1;
        Thetajj       = Thetajj_upd;
        Sigmajj       = Sigmajj_upd;
        llhjj         = llhjj_upd;        
    end
    %% Check llh trend and update residuals precision/covariance
    if k_outer == 1
        xixi_on       = xixi_on + 1;
        sigma2xi      = sigma2xi_upd;
        theta2xi      = theta2xi_upd;
        llhxixi       = llhxixi_upd;       
    elseif llhjj_nupd + llhxixi_upd > llhjj + llhxixi
        xixi_on       = xixi_on + 1;
        sigma2xi      = sigma2xi_upd;
        theta2xi      = theta2xi_upd;
        llhxixi       = llhxixi_upd;
    end
    %%
    %% Likelihood
    llh(k_outer)      = llhjj + llhxixi;
    if (k_outer > 1) && (abs(llh(k_outer) - llh(k_outer-1))/abs(llh(k_outer-1)) < 1E-8)
        llh(k_outer + 1:end) = llh(k_outer);
        break
    end
    %%
end %iterations outer loop
%%
end