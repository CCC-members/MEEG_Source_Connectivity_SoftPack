function [sol_h_hggm] = InverseSolver_higgs(V_sim,LeadFields,SeedsIdx,Nsegments,sens_system,Fs)
Nsim                = size(V_sim,2);
Nseed               = size(SeedsIdx,1);
Ntpoints            = size(V_sim{1,1}{1,1},2);
sol_h_hggm          = cell(4,Nsim);
penalty             = [1 2 0]; % 1 (lasso) 2 (frobenious) 0 (naive)
param.maxiter_outer = 60;
param.maxiter_inner = 30;
p                   = size(V_sim{1,1}{1,1},1);%sensor number
q                   = Nseed;
param.p             = p;
param.q             = q;
param.Ip            = eye(p);
param.Iq            = eye(q);
param.m             = Nsegments;
aj                  = sqrt(log(q)/Nsegments);                                           
Ajj_diag            = 0;                                                      
Ajj_ndiag           = 1;                                               
Ajj                 = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj            = aj;
param.Ajj           = Ajj;
param.axi           = 1E-4;
param.Axixi         = eye(p);
param.Axixi_inv     = eye(p);
param.ntry          = 0;
param.prew          = 0;
param.nu            = Nsegments;
param.rth1          = 0.7;
param.rth2          = 3.16;
llh                 = cell(1,length(penalty));
%%

process_waitbar = waitbar(0,'Please wait...');
for k_sim = 1:Nsim
    waitbar((k_sim)/(Nsim),process_waitbar,strcat('higgs & two-step solution for Sim #',num2str(k_sim),' - SS:',sens_system));
    disp(['higgs & two-step solution for simulation # ',num2str(k_sim),' '])
    %% Preprocessing
    alphaBand=[9.5,10.5];
    V_filt{k_sim}{1} =calc_band_pass(V_sim{k_sim}{1},alphaBand,Fs);
    V_env{k_sim}{1}  =calc_envelope(V_filt{k_sim}{1},Fs);
    Svv_sim{k_sim}{1}=calc_covariance(V_env{k_sim}{1});
%     V_filt{k_sim}{1} =calc_band_pass(V_sim{k_sim}{1});
%     Svv_sim{k_sim}{1}=calc_covariance(V_filt{k_sim}{1});
    %% connectivity leakage module
    Thetajj_est             = zeros(Nseed,Nseed,length(penalty) + 2);
    Sjj_est                 = zeros(Nseed,Nseed,length(penalty) + 2);
    %% h-hggm
    for k_penalty = 1:length(penalty)
        param.penalty  = penalty(k_penalty);
        [Thetajj_est(:,:,k_penalty),Sjj_est(:,:,k_penalty),llh{k_penalty}] = higgs(Svv_sim{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    end
    %% eloreta + hggm
    param.gamma1        = 0.001;
    param.gamma2        = 0.05;
    param.delta_gamma   = 0.001;
    [Thetajj_est(:,:,4),Sjj_est(:,:,4),gamma_grid,gamma,gcv] = eloreta_hg_lasso(Svv_sim{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    %%
    
      
    %% lcmv + hggm
    param.gamma         = sum(abs(diag(Svv_sim{k_sim}{1})))/(length(Svv_sim{k_sim}{1})*100);
    [Thetajj_est(:,:,5),Sjj_est(:,:,5)] = lcmv_hg_lasso(Svv_sim{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);    
        
    
    %%
    sol_h_hggm{1,k_sim}       = SeedsIdx;
    %%
    sol_h_hggm{2,k_sim}       = Sjj_est;
    %%
    sol_h_hggm{3,k_sim}       = Thetajj_est;
    %%
    sol_h_hggm{4,k_sim}       = {llh,gamma_grid,gcv,gamma,param.gamma};
end
delete(process_waitbar);