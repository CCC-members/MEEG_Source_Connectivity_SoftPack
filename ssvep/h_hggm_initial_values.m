function [sigma2xi,theta2xi,p,Ip,Axixi_inv,Lvj,Ljv,Svv,Sigmajj,Thetajj,Iq,Ajj,aj,llh] = h_hggm_initial_values(Svv,Lvj,m,maxiter_outer,Axixi)
%% Initialization of variables and tunning parameters
[p,q]           = size(Lvj);
llh             = zeros(maxiter_outer,1);                                    % Likelihood h_hggm
sigma2xi        = 1;                                                         % residuals (noise) variance
theta2xi        = 1/sigma2xi;                                                % residuals (noise) precision
Ip              = eye(p);                                                    % Identity matrix on the sensors space
Iq              = eye(q);                                                    % Identity matrix on the cortical generators space
aj              = sqrt(log(q)/m);                                            % Jankova regularization parameter
Ajj_diag        = 0;                                                         % Regularization mask diagonal parameter
Ajj_ndiag       = 1;                                                         % Regularization mask nondiagonal parameter
Ajj             = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));                % Regularization mask
Axixi_inv       = Ip/Axixi;                                                  % Residuals (noise) Covariance
%% Lead Field scaling
scale_Lvj       = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj             = Lvj/scale_Lvj;
Ljv             = Lvj';
%% Transfer Function and Posterior Covariance
Sigmajj         = Iq;
Thetajj         = Iq/Sigmajj;
SigmajjLjv      = Sigmajj*Ljv;
LvjSigmajj      = Lvj*Sigmajj;
Sigmajj_pst     = Sigmajj - SigmajjLjv*(Ip/(Lvj*SigmajjLjv+sigma2xi*Axixi_inv))*LvjSigmajj; 
Sigmajj_pst     = (Sigmajj_pst + Sigmajj_pst')/2;
Tjv             = (Sigmajj_pst*Ljv)*Axixi*theta2xi;
Tvj             = Tjv';
%% Computation of sources activity calibration empirical covariance
Sjj             = Tjv*Svv*Tvj; 
Sjj             = (Sjj + Sjj')/2;
%% Scaling
scale           = (sum(abs(diag(Sjj)))/q)/(max(abs(diag(Sigmajj_pst))));
Svv             = Svv/scale;
end