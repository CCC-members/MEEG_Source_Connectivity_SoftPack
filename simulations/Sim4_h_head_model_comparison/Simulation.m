function [] = Simulation(sens_system,output_sourse)
%% initialization
[LeadFields,Nv,Nsubj,Nsamp,Nnoise,db_source,db_sens,Nsim,Nseed,options,index_seed,index_full,vertices,faces,elec_pos] = Initialize_simulation(sens_system);
%% Run Simulation for specific configurations of sources and subjects
J_sim       = cell(1,Nsim);
Svv0_sim    = cell(1,Nsim);
Svv_sim     = cell(1,Nsim);
Theta_sim   = cell(1,Nsim);
Seeders_sim = zeros(Nseed,Nsim);

process_waitbar1 = waitbar(0,'Please wait...');

for sim = 1:Nsim
    waitbar((sim)/(Nsim),process_waitbar1,strcat('simulation # ',num2str(sim), '  to sens-system: ',sens_system));
    disp(['simulation # ',num2str(sim)])
    % Setting correlation structure between points
%     [S,Data,X]              = gen_hggm1(Nsamp,sum(options.extensions),2,options);
%     Data                    = transpose(Data);
    [S,Data,X]              = gen_hggm2(Nsamp,sum(options.extensions),options);
    %% Simulating Sensors Time Series Data
    V0                      = cell(1,Nsubj);
    for cont_seed = 1:Nseed        
        Seeders_sim(cont_seed,sim) = index_seed{cont_seed}(randi(length(index_seed{cont_seed})));
    end
     Seeders_sim(:,sim) = Seeders_sim(randperm(Nseed),sim); 
    for cont4 = 1:Nsubj
        K  = LeadFields{cont4};
        V0{cont4}  = K(:,Seeders_sim(:,sim))*Data;
    end
    %% Computing Sources Variances
    J_tmp                   = sum(abs(Data).^2,2);
    J                       = zeros(Nv,1);
    J(Seeders_sim(:,sim),:) = J_tmp/Nsamp; % 3 source
    %% Plot
%     figure
%     patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',J/max(J),'FaceColor','interp','EdgeColor','b','FaceAlpha',.95);
%     colormap('hot');
%     axis off
%     set(gcf,'Color','k');
%     caxis([0 1]);
    %% Simulating Biological Noise
    rs0                     = randn(Nnoise,Nsamp,Nsubj) + 1i*randn(Nnoise,Nsamp,Nsubj);
    %% Projecting Biological Noise to Sensors
    noisesources            = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        K                   = LeadFields{cont4};
        Ne                  = size(K,1);
        rs_tmp              = K(:,index_full)*squeeze(rs0(:,:,cont4));
        noisesources{cont4} = db_source*sum(abs(V0{cont4}(:)).^2)^(1/2)*rs_tmp/sum(abs(rs_tmp(:)).^2)^(1/2);
    end
    %% Simulating Sensors Noise
    noisesensors            = cell(1,Nsubj);
    for cont4 = 1:Nsubj
        K                   = LeadFields{cont4};
        Ne                  = size(K,1);
        rs_tmp              = randn(Ne,Nsamp) + 1i*randn(Ne,Nsamp);
        noisesensors{cont4} = db_sens*sum(abs(V0{cont4}(:)).^2)^(1/2)*rs_tmp/sum(abs(rs_tmp(:)).^2)^(1/2);
    end
    %% Simulating Data by K*J + SensorsNoise + BiologicalNoise
    V = cell(1,Nsubj);
    for cont4 = 1:Nsubj
    V{cont4} = V0{cont4} + noisesources{cont4} + noisesensors{cont4};
    end
    %% Computing Covariance across Data samples
    Svv  = cell(1,Nsubj);
    Svv0 = cell(1,Nsubj);
    for cont4 = 1:Nsubj
            Svv0{cont4}     = squeeze(V0{cont4})*squeeze(V0{cont4})'/Nsamp;
            Svv{cont4}      = squeeze(V{cont4})*squeeze(V{cont4})'/Nsamp;
    end
    %%
    J_sim{:,sim}            = J;
    Svv0_sim{:,sim}         = Svv0;
    Svv_sim{:,sim}          = Svv;
    Theta_sim{:,sim}        = X;
end
delete(process_waitbar1);
save(strcat('simulations',filesep,'Sim4_h_head_model_comparison',filesep,'Pseudorand_Net'),'J_sim', 'Seeders_sim', 'index_full', 'Svv0_sim', 'Svv_sim' ,'Theta_sim', 'vertices', 'faces', 'elec_pos', 'LeadFields', 'Nsubj', 'Nsamp', 'sens_system', 'options')

end