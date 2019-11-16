function [] = Simulation(sens_system,output_source)
FLAG_SAVESOURCE = 0;
%% initialization
[LeadFields,Nsensor,Nsource,Nsubj,Nsim,Nsegments,Nseed,Nnoise,Nfreqs,Ntpoints,db_source,db_sens,options,index_seed,index_full] = Initialize_simulation(sens_system);
%% Run Simulation for specific configurations of sources and subjects
SeedsIdx                        = zeros(Nseed,Nsim);
SeedsPos                        = zeros(Nseed,size(options.vertices,2),Nsim);
Jvar_sim                        = cell(1,Nsim);
Svv0_sim                        = cell(1,Nsim);
Svv_sim                         = cell(1,Nsim);
V0_sim                          = cell(1,Nsim);
V_sim                           = cell(1,Nsim);
Jseed_sim                       = cell(1,Nsim);
sigma_sim                       = cell(1,Nsim);
theta_sim                       = cell(1,Nsim);
Jwhole_sim                      = cell(1,Nsim);
for sim=1:Nsim
    Svv0_sim{:,sim}             = cell(1,Nsubj); % Sensor covariances without noise
    Svv_sim{:,sim}              = cell(1,Nsubj);  % Sensor covariances with noise
    V0_sim{:,sim}               = cell(1,Nsubj);% Sensor space without noise
    V_sim{:,sim}                = cell(1,Nsubj); % Sensor space with noise
end
for sim=1:Nsim
    Jvar_sim{:,sim}              = zeros(Nsource,1);    % Sources Variances
    Jseed_sim{:,sim}             = zeros(Nseed,Ntpoints,Nsegments);       % Activate source signal
    theta_sim{:,sim}             = zeros(Nseed,Nseed,Nfreqs);     % Precision matrix
    sigma_sim{:,sim}             = zeros(Nseed,Nseed,Nfreqs);  % Covariance matrix
    if FLAG_SAVESOURCE==1
        Jwhole_sim{:,sim} = ndSparse.build([Nsource,Ntpoints,Nsegments]);    % Whole brain source signal without noise
    end
    for subj=1:Nsubj
        
        Svv0_sim{:,sim}{:,subj}  = zeros(Nsensor,Nsensor); % Sensor covariances without noise
        Svv_sim{:,sim}{:,subj}   = zeros(Nsensor,Nsensor);  % Sensor covariances with noise
        V0_sim{:,sim}{:,subj}    = zeros(Nsensor,Ntpoints);% Sensor space without noise
        V_sim{:,sim}{:,subj}     = zeros(Nsensor,Ntpoints); % Sensor space with noise
    end
end
% process_waitbar1 = waitbar(0,'Please wait...');
%%
for sim = 1:Nsim
    % waitbar((sim)/(Nsim),process_waitbar1,strcat('simulation # ',num2str(sim), '  to sens-system: ',sens_system));
    disp(['simulation # ',num2str(sim)])
    % Seeds
    for cont_seed = 1:Nseed
        SeedsIdx(cont_seed,sim) = index_seed{cont_seed}(randi(length(index_seed{cont_seed})));
    end
    SeedsIdx(:,sim)   = SeedsIdx(randperm(Nseed),sim);
    SeedsPos(:,:,sim) = options.vertices(SeedsIdx(:,sim),:);
    % Setting correlation structure between points
    options.SeedsIdx=SeedsIdx;
    options.SeedsPos=SeedsPos;
    %% Simulating Alpha signal
    [sigma,theta,process,process_eneger] = gen_processes(Nsegments,Nseed,options);
    %% Simulating Biological Noise
    [rs0,noise_eneger] = gen_process_xi(Nnoise,Ntpoints,Nsegments,Nsubj,index_full,options);
    %%
    for segment = 1:Nsegments
        %% Simulating Sensors Time Series Data
        V0       = cell(1,Nsubj);
        for subj = 1:Nsubj
            K        = LeadFields{subj};
            V0{subj} = K(:,SeedsIdx(:,sim))*process(:,:,segment);
        end
        %% Projecting Biological Noise to Sensors
        noisesources            = cell(1,Nsubj);
        for subj = 1:Nsubj
            K                   = LeadFields{subj};
            rs_tmp              = K(:,index_full)*squeeze(rs0(:,:,segment,subj));
            %                         noisesources{subj}  = db_source*sum(abs(V0{subj}(:)).^2)^(1/2)*rs_tmp/sum(abs(rs_tmp(:)).^2)^(1/2);
%             SourceNoiseRatio1=db_source*sum(abs(V0{subj}(:)).^2)^(1/2)/sum(abs(rs_tmp(:)).^2)^(1/2);
            SourceNoiseRatio=db_source* sum(abs(K(:,SeedsIdx(:,sim))*process_eneger).^2,[1,2])^(1/2) / sum(abs(K(:,index_full)*noise_eneger).^2,[1,2])^(1/2);

%             SourceNoiseRatio1/SourceNoiseRatio;
            noisesources{subj}  = SourceNoiseRatio*rs_tmp;
        end
        %% Simulating Sensors Noise
        noisesensors            = cell(1,Nsubj);
        for subj = 1:Nsubj
            K                   = LeadFields{subj};
            Ne                  = size(K,1);
            rs_tmp              = randn(Ne,Ntpoints);
            noisesensors{subj}  = db_sens*sum(abs(V0{subj}(:)).^2)^(1/2)*rs_tmp/sum(abs(rs_tmp(:)).^2)^(1/2);
        end
        %% Simulating Data by K*J + SensorsNoise + BiologicalNoise
        V = cell(1,Nsubj);
        for subj = 1:Nsubj
            V{subj} = V0{subj} + noisesources{subj} + noisesensors{subj};
        end

        %% Computing Covariance across Data samples
        Svv  = cell(1,Nsubj);
        Svv0 = cell(1,Nsubj);
        for subj = 1:Nsubj
            Svv0{subj}     = squeeze(V0{subj})*squeeze(V0{subj})'/Ntpoints;
            Svv{subj}      = squeeze(V{subj})*squeeze(V{subj})'/Ntpoints;
        end
        %%
        if isa(process,'gpuArray')
            for subj = 1:Nsubj
                Svv0_sim{:,sim}{subj}(:,:,segment)         = gather(Svv0{subj}); % Sensor covariances without noise
                Svv_sim{:,sim}{subj}(:,:,segment)          = gather(Svv{subj});  % Sensor covariances with noise
                V0_sim{:,sim}{subj}(:,:,segment)           = gather(V0{subj});   % Sensor space without noise
                V_sim{:,sim}{subj}(:,:,segment)            = gather(V{subj});    % Sensor space with noise
            end
        else
            for subj = 1:Nsubj
                Svv0_sim{:,sim}{subj}(:,:,segment)         = Svv0{subj}; % Sensor covariances without noise
                Svv_sim{:,sim}{subj}(:,:,segment)          = Svv{subj};  % Sensor covariances with noise
                V0_sim{:,sim}{subj}(:,:,segment)           = V0{subj};   % Sensor space without noise
                V_sim{:,sim}{subj}(:,:,segment)            = V{subj};    % Sensor space with noise
            end
        end
    end %for segment
    %% Computing Sources Variances
    Jvar_tmp                      = sqrt(sum(squeeze(sum(process.^2,3)),2));
    Jvar                          = zeros(Nsource,1);
    if test_gpu
        Jvar(SeedsIdx(:,sim),:)   = gather(Jvar_tmp/(Ntpoints*Nsegments)); % 3 source % Variances
    else
        Jvar(SeedsIdx(:,sim),:)   = Jvar_tmp/(Ntpoints*Nsegments); % 3 source % Variances
    end
    %%
    if isa(process,'gpuArray')
        Jvar_sim{:,sim}                   = gather(Jvar);    % Sources Variances
        Jseed_sim{:,sim}                  = gather(process);              % Seeds signal
        theta_sim{:,sim}                  = gather(theta);          % Precision matrix
        sigma_sim{:,sim}                  = gather(sigma);
        if FLAG_SAVESOURCE==1
            Jwhole_sim{:,sim}(SeedsIdx(:,sim),:,:)             = gather(process);            % Source signal without noise
        else
            Jwhole_sim{:,sim}(SeedsIdx(:,sim),:,:)             =nan;
        end
    else
        Jvar_sim{:,sim}                   = Jvar;    % Sources Variances
        Jseed_sim{:,sim}                  = process;              % Seeds signal
        theta_sim{:,sim}                  = theta;          % Precision matrix
        sigma_sim{:,sim}                  = sigma;
        if FLAG_SAVESOURCE==1
            Jwhole_sim{:,sim}(SeedsIdx(:,sim),:,:)        = process;            % Source signal without noise
        else
            Jwhole_sim{:,sim}(SeedsIdx(:,sim),:,:)        =nan;
        end
    end
end %for sim
% % delete(process_waitbar1);

save([output_source,filesep,'Pseudorand_Net.mat'],'Jvar_sim', 'index_full', 'Svv0_sim', 'Svv_sim' ,'LeadFields', ...
    'sens_system', 'options','V0_sim','V_sim','K','Jwhole_sim','sigma_sim','theta_sim','Jseed_sim','options','-v7.3')
end
