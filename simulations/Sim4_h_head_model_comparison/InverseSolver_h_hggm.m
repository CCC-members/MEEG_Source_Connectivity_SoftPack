function [sol_h_hggm] = InverseSolver_h_hggm(Svv_sim,LeadFields,Seeders_sim,Nsamp,sens_system)
Nsim                = size(Svv_sim,2);
Nseed               = size(Seeders_sim,1);
sol_h_hggm          = cell(4,Nsim);
penalty             = [1 2 0]; % 1 (lasso) 2 (frobenious) 0 (naive)
param.maxiter_outer = 60;
param.maxiter_inner = 30;
param.m             = Nsamp;
param.rth           = 3.16;
param.axi           = 1E-5;
param.Axixi         = eye(length(Svv_sim{1,1}{1,1}));
llh                 = cell(1,length(penalty));
%%

process_waitbar = waitbar(0,'Please wait...');
for k_sim = 1:Nsim
    waitbar((k_sim)/(Nsim),process_waitbar,strcat('h-hggm & eloreta-hggm solution for Sim #',num2str(k_sim),' - SS:',sens_system));
    disp(['h-hggm and eloreta-hggm solution for simulation # ',num2str(k_sim),' '])
    %% connectivity leakage module
    Thetajj_est             = zeros(Nseed,Nseed,length(penalty) + 2);
    Sjj_est                 = zeros(Nseed,Nseed,length(penalty) + 2);
    %% h-hggm
    for k_penalty = 1:length(penalty)
        param.penalty  = penalty(k_penalty);
        [Thetajj_est(:,:,k_penalty),Sjj_est(:,:,k_penalty),llh{k_penalty},jj_on,xixi_on] = h_hggm(Svv_sim{k_sim}{1},LeadFields{1}(:,Seeders_sim(:,k_sim)),param);
    end
    %% eloreta + hggm
    param.gamma1        = 0.001;
    param.gamma2        = 0.05;
    param.delta_gamma   = 0.001;
    [Thetajj_est(:,:,4),Sjj_est(:,:,4),gamma_grid,gamma,gcv] = eloreta_hggm(Svv_sim{k_sim}{1},LeadFields{1}(:,Seeders_sim(:,k_sim)),param);
    %%
    
      
    %% lcmv + hggm
    param.gamma         = sum(abs(diag(Svv_sim{k_sim}{1})))/(length(Svv_sim{k_sim}{1})*100);
    [Thetajj_est(:,:,5),Sjj_est(:,:,5)] = lcmv_hggm(Svv_sim{k_sim}{1},LeadFields{1}(:,Seeders_sim(:,k_sim)),param);    
        
    
    %%
    sol_h_hggm{1,k_sim}       = Seeders_sim;
    %%
    sol_h_hggm{2,k_sim}       = Sjj_est;
    %%
    sol_h_hggm{3,k_sim}       = Thetajj_est;
    %%
    sol_h_hggm{4,k_sim}       = {llh,gamma_grid,gcv,gamma};
end
delete(process_waitbar);