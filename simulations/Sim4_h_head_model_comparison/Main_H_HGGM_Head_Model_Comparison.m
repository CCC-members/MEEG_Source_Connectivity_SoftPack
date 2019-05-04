function result = Main_H_HGGM_Head_Model_Comparison(output_sourse)


%%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: April 4, 2019


% Updates
% - Ariosky Areces Gonzalez

% Date: April 4, 2019


output_sourse = strcat(output_sourse, filesep,'H_HGGM_Head_Model_Comparison');
if(~isfolder(output_sourse))
mkdir(output_sourse);
end

%% Input Lead Field model

for i = 1:3
    
    if(i == 1)
        sens_system = 'pseudo';
    elseif(i == 2)
        sens_system = 'small';
    else
        sens_system = 'large';
    end
   
    
    %% Simulation
    Simulation(sens_system, output_sourse);
    

%% Inverse Solvers

InverseSolvers(output_sourse);
%%
%% Results

Results(output_sourse);

end
