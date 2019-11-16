function result = Main_HIGGS_TimeDomain(output_source)
%%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: April 4, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: April 4, 2019

rng('default')

% addpath([fileparts(mfilename('fullpath')),filesep,'function'])
addpath(genpath([fileparts(mfilename('fullpath')),filesep,'external']))
if nargin<1 || isempty(output_source)
    output_source='.\result';
end
output_source = strcat(output_source, filesep,'HIGGS_TimeDomain');
if(~exist(output_source))
    mkdir(output_source);
end

%% Input Lead Field model

for i = 2:2
    if(i == 1)
        sens_system = 'pseudo';
    elseif(i == 2)
        sens_system = 'small';
    else
        sens_system = 'large';
    end
    %% Simulation
%     Simulation(sens_system, output_source);
    %% Inverse Solvers
    InverseSolvers(output_source);
    %% Results
    Results(output_source);
end
