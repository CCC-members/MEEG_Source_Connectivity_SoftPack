function [] = InverseSolvers(output_sourse)

%% Run Inverse Solvers
%% Loading Simulation Substrate or Real Data.

load('Pseudorand_Net.mat')
%%
if strcmp(sens_system,'large') == 1
    load ('LeadFields_large.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
elseif strcmp(sens_system,'small') == 1
    load ('LeadFields_small.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
elseif strcmp(sens_system,'pseudo') == 1
    load ('LeadFields_pseudo.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
end
%%
subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
LeadFields = LeadFields(1,subject);
% 'Vsim' is a cell array containing Simulated or Real Data at different conditions,
% every cell contains a tensor of (Number of Sensors, Time points or Frequencies, Subjects X Samples)
%% Run h-hggm
sol_h_hggm = InverseSolver_h_hggm(Svv_sim,LeadFields,Seeders_sim,Nsamp,sens_system);

save(strcat('simulations',filesep,'Sim4_h_head_model_comparison',filesep,'Solutions_h_hggm'),'sol_h_hggm', '-v7.3')

end