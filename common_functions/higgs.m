function [Thetajj,Sjj,llh] = higgs(Svv,Lvj,param)
% Hidden Gaussian Grpahical Source-Model (HIGGS) solver. Computes the Source Empirical Covariance (Sjj) and Source Partial Correlations (Thetajj) by two sequential steps.
% First: Unhides (Expectation) the Type II Likelihood approximated representation which shapes a pair of Hermitian Gaussian Graphical Model (HGGM), one of the state equation
% (with empirical covariance Psijj) and another of the observation equation residuals (with empirical covariance Psixixi). The Hyperparameters are computed by maximum posterior
% analysis (Maximization) regularized with priors. The states (source) HGGM is estimated with hermitian graphical lasso (HG-LASSO) solver controled by Rayleigh threshold
% and the residuals HGGM with exponential prior of the noise presicion controled by nuisance inferior limit (scale parameter).
%
% inputs:
%    Svv                 : M/EEG time series cross-spectra
%    Lvj                 : M/EEG Lead Field
%    param               : set of parameters
%
% outputs:
%    Thetajj             : source partial coherence
%    Sjj                 : source empirical covariance
%    llh                 : likelihood of the h-hggm unhide and solve (em) loop
%
% Pedro Valdes-Sosa, March 2019
% Deirel Paz Linares, March 2019
% Eduardo Gonzalez-Moreira, March 2019
%************************************************************************** 
%% Initialization EM algorithm
[Svv,Lvj,Ljv,scale,sigma2xi0,Sigmajj0,llh0,param] = higgs_initial_values(Svv,Lvj,param);
%% Outer loop EM algorithm
maxiter_outer         = param.maxiter_outer;
llh                   = zeros(maxiter_outer,1);
for k_outer = 1:maxiter_outer
    disp(['higgs expectation-maximization loop # ',num2str(k_outer)])
    %% Expectation
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
    %%
    %% Maximization
    [sigma2xi,Sigmajj,Thetajj]            = higgs_maximization(Psixixi,Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
    %% Compute the global hyperparamters posterior log-likelihood
    [llh(k_outer)]                        = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_pst,param);
    llh0                                  = llh(k_outer);
    sigma2xi0                             = sigma2xi;
    Sigmajj0                              = Sigmajj;
    %% Stopping criteria
    if (k_outer > 1) && ((abs(llh(k_outer) - llh(k_outer-1))/abs(llh(k_outer-1)) < 1E-4) || (llh(k_outer) < llh(k_outer-1))) 
        llh(k_outer:end) = llh(k_outer-1);
        break
    end
end %iterations outer loop
%%
Sjj     = Sjj/scale;
Thetajj = Thetajj*scale;
end