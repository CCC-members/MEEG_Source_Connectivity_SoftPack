%% Time domain simulation
% Authors:
% - Ying Wang
% - Ariosky
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa
% Date: Dec 17, 2019
function result = Main_HIGGS_TimeDomain_Demo(output_source)
%%
rng('default')
addpath(genpath([fileparts(mfilename('fullpath')),filesep,'external']))
if nargin<1 || isempty(output_source)
    output_source='./result';
end
output_source = strcat(output_source, filesep,'HIGGS_TimeDomain');
if(~exist(output_source))
    mkdir(output_source);
end
%%
i = 2;
if(i == 1)
    sens_system = 'pseudo';
elseif(i == 2)
    sens_system = 'small';
else
    sens_system = 'large';
end
%%
[LeadFields,Nsensor,Nsource,Nsubj,Nsim,Nsegments,Nseed,Nnoise,Nfreqs,Ntpoints,db_source,db_sens,options,index_seed,index_full] = Initialize_simulation(sens_system);
%%
V_sim                           = cell(1,Nsim);
for sim=1:Nsim
    V_sim{:,sim}                = cell(1,Nsubj); % Sensor space with noise
end
for sim=1:Nsim
    for subj=1:Nsubj
        V_sim{:,sim}{:,subj}     = zeros(Nsensor,Ntpoints); % Sensor space with noise
    end
end
%%
for sim = 1:Nsim
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
    %%
    [~,theta0]     = gen_hggm02(Nseed,options);
    %%
    % creationg a covariance/precision tensor with specific spectral form and assumptions of the underlying spectral factorization,
    % the assumptions of the spectral factorization are encoded by the function "apply tensor phase"
    % Initializing
    F           = options.F;
    Iq          = eye(Nseed);
    theta       = zeros(Nseed,Nseed,length(F));
    sigma       = zeros(Nseed,Nseed,length(F));
    % crating tensors
    phi = tstudent_spectra(F, options.tstudent_a(2),options.tstudent_b(2),options.tstudent_dof(2),options.nu0(2));
    phi = phi/max(phi);%scaling to 0-1;
%     K = apply_tensor_phase(theta0, options.nu0, options.deltaf, options.Fmin, options.Fmax);
    [Nseed,~]         = size(theta0);
    [U0,D]        = svd(theta0);
    D             = sqrt(D);
    U0            = U0*D;
    Iq            = eye(Nseed);
    K0            = Iq - U0';
    pha0          = angle(K0);
    F             = options.Fmin:options.deltaf:options.Fmax;
    pha_rate      = F/options.nu0(2);
    pha           = bsxfun(@times,pha0,reshape(pha_rate,1,1,length(F)));
    K             = repmat(abs(K0),1,1,length(F)).*exp(1i *(pha));
    K = K./repmat(reshape(phi,1,1,length(F)),Nseed,Nseed,1);
    U = Iq - K;
    for count = 1:length(F)
        theta(:,:,count) = U(:,:,count)'*U(:,:,count);
        sigma(:,:,count) = inv(theta(:,:,count));
    end
    % set DC coef as zero
    theta(:,:,1)=0;
    sigma(:,:,1)=0;
    %
    process_energy = sigma_to_band_energy((options.Ntpoints)*sigma,...
        [options.band(2,1),options.band(2,2)],options.deltaf);
    [fourier_coef] = gen_tensor_hggm((options.Ntpoints)*sigma,Nsegments);% (options.Ntpoints)*
    process        = iFFT_time_series(fourier_coef);
    %%
    % Xi process
    Axi= gen_NnvarAdjacent(options.vertices,options.faces);
    Axi=Axi(index_full,index_full);
    % VAR
    NwarmTimepoints=100;
    for i=1:Nsubj
        rs0(:,:,:,i)=gen_var_time_series(full(Axi),eye(size(Axi,1)),Ntpoints,Nsegments,options.Fs,NwarmTimepoints);
    end
    %%
    band_range=(pi/(options.Nfreqs-1))*[options.band(2,1)/options.deltaf:1:options.band(2,2)/options.deltaf];
    H=calc_var_to_transfe(Axi,band_range);
    for i=1:size(H,3)
        Sxi(i)=norm(H(:,:,i),'fro').^2/size(Axi,1);% because the coefficient is a toeplitz matrix, so all nodes have same amplitude
    end
    noise_energy=sqrt(Sxi*options.Ntpoints);
    noise_energy= repmat(noise_energy,[Nnoise,1]);
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
            SourceNoiseRatio=db_source* sum(sum(abs(K(:,SeedsIdx(:,sim))*process_energy).^2))^(1/2) / sum(sum(abs(K(:,index_full)*noise_energy).^2))^(1/2);
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
        %%
        if isa(process,'gpuArray')
            for subj = 1:Nsubj
                V_sim{:,sim}{subj}(:,:,segment)            = gather(V{subj});    % Sensor space with noise
            end
        else
            for subj = 1:Nsubj
                V_sim{:,sim}{subj}(:,:,segment)            = V{subj};    % Sensor space with noise
            end
        end
    end %for segment
end



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
for k_sim = 1:Nsim
    disp(['higgs & two-step solution for simulation # ',num2str(k_sim),' '])
    %% Preprocessing
    alphaBand=[9.5,10.5];
    V_filt{k_sim}{1}   = calc_band_pass(V_sim{k_sim}{1},alphaBand,options.Fs);
    V_env{k_sim}{1}    = calc_envelope(V_filt{k_sim}{1},options.Fs);
    Svv_filt{k_sim}{1} = calc_covariance(V_filt{k_sim}{1});
    Svv_env{k_sim}{1}  = calc_covariance(V_env{k_sim}{1});
    %% connectivity leakage module
    Thetajj_est             = zeros(q,q,length(penalty) + 2);
    Sjj_est                 = zeros(q,q,length(penalty) + 2);
    %% h-hggm
    for k_penalty = 1:length(penalty)
        param.penalty  = penalty(k_penalty);
        [Thetajj_est(:,:,k_penalty),Sjj_est(:,:,k_penalty),llh{k_penalty}] = higgs(Svv_env{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    end
    %% eloreta + hggm
    param.gamma1        = 0.001;
    param.gamma2        = 0.05;
    param.delta_gamma   = 0.001;
    [Thetajj_est(:,:,4),Sjj_est(:,:,4),gamma_grid,gamma,gcv] = eloreta_hg_lasso(Svv_env{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    %% lcmv + hggm
    param.gamma         = sum(abs(diag(Svv_env{k_sim}{1})))/(length(Svv_env{k_sim}{1})*100);
    [Thetajj_est(:,:,5),Sjj_est(:,:,5)] = lcmv_hg_lasso(Svv_env{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    %% eloreta + Roi-nets
    [~,~,~,~,~,Tjv{1}] = eloreta_hg_lasso(Svv_filt{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    Thetajj_est(:,:,6)=roi_nets_network_analysis(V_filt{k_sim}{1},Tjv{1},options.Fs,aj);
    %% lcmv + Roi-nets
    [~,~,Tjv{2}] = lcmv_hg_lasso(Svv_filt{k_sim}{1},LeadFields{1}(:,SeedsIdx(:,k_sim)),param);
    Thetajj_est(:,:,7)=roi_nets_network_analysis(V_filt{k_sim}{1},Tjv{2},options.Fs,aj);
end
end










