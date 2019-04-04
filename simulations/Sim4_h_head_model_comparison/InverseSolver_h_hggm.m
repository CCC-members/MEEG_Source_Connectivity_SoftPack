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
     waitbar((k_sim)/(Nsim),process_waitbar,strcat('h-hggm solution for simulation # ',num2str(k_sim),'-  sens-system: ',sens_system));
    disp(['h_hggm solution for simulation # ',num2str(k_sim)])
    %% connectivity leakage module
    Thetajj_est             = zeros(Nseed,Nseed,length(penalty));
    Sjj_est                 = zeros(Nseed,Nseed,length(penalty));
    for k_penalty = 1:length(penalty)
        param.penalty  = penalty(k_penalty);
        [Thetajj_est(:,:,k_penalty),Sjj_est(:,:,k_penalty),llh{k_penalty},jj_on,xixi_on] = h_hggm(Svv_sim{k_sim}{1},LeadFields{1}(:,Seeders_sim(:,k_sim)),param);
    end
    %%
    sol_h_hggm{1,k_sim}       = Seeders_sim;
    %%
    sol_h_hggm{2,k_sim}       = Sjj_est;
    %%
    sol_h_hggm{3,k_sim}       = Thetajj_est;
    %%
    sol_h_hggm{4,k_sim}       = llh;
end
delete(process_waitbar);