function result = Main_realistic_em_penalty_test(output_sourse)


%%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: April 4, 2019


% Updates
% - Ariosky Areces Gonzalez

% Date: April 4, 2019


%% Test Penalties


output_sourse = strcat(output_sourse, filesep,'Realistic_em_penalty_test');
if(~isfolder(output_sourse))
mkdir(output_sourse);
end
 
for i = 1:1
    process_waitbar = waitbar(0,'Please wait...');
    if(i == 1)
        sens_system = 'small';
    else
        sens_system = 'large';
    end
    waitbar(1/4,process_waitbar,strcat('higgs & two-step solution - SS:',sens_system ));
    %%
    if strcmp(sens_system,'large') == 1
        load('HeadModel_large.mat'); % Load mesh vertices (all Lead Fields should be given on homemorph surfaces)
    elseif strcmp(sens_system,'small') == 1
        load('HeadModel_small.mat'); % Load mesh vertices (all Lead Fields should be given on homemorph surfaces)
    end
    %%
    vertices   = cortex.vertices;
    faces      = cortex.faces;
    %%
    waitbar(2/4,process_waitbar,strcat('higgs & two-step solution - SS:',sens_system ));
    %%  Definition of Cortical points
    if strcmp(sens_system,'large') == 1
        load('data_tips24_large.mat')
    elseif strcmp(sens_system,'small') == 1
        load('data_tips22_small.mat')
    end
    %%
    Nseed        = size(data_tips,2);
    Seeders      = zeros(Nseed,1);
    for cont = 1:Nseed
        coord_tmp = data_tips(cont).Position;
        vx       = coord_tmp(1);
        vy       = coord_tmp(2);
        vz       = coord_tmp(3);
        Seeders(cont) = pickpoint(vx,vy,vz,vertices,1E-3);
    end
    Seeders      = Seeders(randperm(Nseed));
    waitbar(3/4,process_waitbar,strcat('higgs & two-step solution - SS:',sens_system ));
    %%
    if strcmp(sens_system,'large') == 1
        load ('LeadFields_large.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
        subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
        Lvj        = LeadFields{1,subject};
    elseif strcmp(sens_system,'small') == 1
        load ('LeadFields_small.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
        subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
        Lvj        = LeadFields{1,subject};
    end
    %%
    if strcmp(sens_system,'large') == 1
        m               = 6000;                          % Sample number
    elseif strcmp(sens_system,'small') == 1
        m               = 600;                          % Sample number
    end
    q                   = size(Seeders,1);              % Number of generators
    p                   = size(Lvj,1);                  % Number of sensors
    nblocks             = 2;                            % Number of blocks in simulation
    options.config      = 2;                            % (2) overlapping blocks (1) nonoverlapping blocks
    options.var         = 2;                            % (2) complex variable (1) real variable
    options.extensions  = [ceil(Nseed/3); ceil(Nseed/3); Nseed - 2*ceil(Nseed/3)];  % patches extensions
    options.connections = [1 2; 2 3];        % patches connections
    % [Sjj_sim,j_sim,Thetajj_sim] = gen_hggm1(m,q,nblocks,options);
    % j_sim = transpose(j_sim);
    [Sjj_sim,j_sim,Thetajj_sim] = gen_hggm2(m,q,options);
    %%
    %% Generate data
    v0             = Lvj(:,Seeders)*j_sim; % data (observation)
    %% Biological noise
    d0  = 5E-3;
    index_full = [];
    for point = 1:Nseed
        Source          = Seeders(point);
        [index,findex]  = surfpatch(Source,vertices,faces,d0);
        index_full      = [index_full;index];
    end
    
    Nnoise              = length(index_full);
    bionoise            = randn(Nnoise,m) + 1i*randn(Nnoise,m);
    bionoise            = Lvj(:,index_full)*bionoise;
    bionoise            = sum(abs(v0(:)).^2)^(1/2)*bionoise/sum(abs(bionoise(:)).^2)^(1/2);
    %% Sensor noise
    sensnoise           = randn(p,m) + 1i*randn(p,m);
    sensnoise           = sum(abs(v0(:)).^2)^(1/2)*sensnoise/sum(abs(sensnoise(:)).^2)^(1/2);
    %% Corrupted data
    v                   = v0 + 0.1*bionoise + 0.1*sensnoise;
    %% Data empirical covariance
    Svv                 = cov(v');
    %% Likelihood test
    param.use_gpu       = 0;
    param.run_bash_mode = 0;
    param.str_band      = "none";
    param.maxiter_outer = 60;
    param.maxiter_inner = 30;
    p                   = length(Svv);
    param.p             = p;
    param.q             = q;
    param.Ip            = eye(p);
    param.Op            = ones(p,1);
    param.Iq            = eye(q);
    param.m             = m;
    aj                  = sqrt(log(q)/m);
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
    param.nu            = m;
    param.rth1          = 0.7;
    param.rth2          = 3.16;
    param.eigreg        = 1E-4;
    %%
    penalty             = [1 2 0]; % 1 (lasso) 2 (frobenious) 0 (naive)
    Thetajj_est         = zeros(q,q,length(penalty) + 2);
    llh_outer           = cell(1,length(penalty));
    llh_inner           = cell(1,length(penalty));
    waitbar(4/4,process_waitbar,strcat('higgs & two-step solution - SS:',sens_system ));
    for cont = 1:length(penalty)       
        param.penalty  = penalty(cont);
        [Thetajj_est(:,:,cont),Sjj_est,Tjv_est,llh_outer{cont}] = higgs(Svv,Lvj(:,Seeders),param);
    end
    
    %% eloreta + hggm
    param.gamma1        = 0.001;
    param.gamma2        = 0.05;
    param.delta_gamma   = 0.001;
    [Thetajj_est(:,:,4),Sjj_est,gamma_grid,gamma,gcv] = eloreta_hg_lasso(Svv,Lvj(:,Seeders),param);
    
    %% lcmv + hggm
    param.gamma         = sum(abs(diag(Svv)))/(length(Svv)*100);
    [Thetajj_est(:,:,5),Sjj_est] = lcmv_hg_lasso(Svv,Lvj(:,Seeders),param);
    delete(process_waitbar);
    
    %% Plot Results
    figure_partial_coherence_maps = figure('Position',[182,114,832,521]);
    load('colormap3')
    %%
    X = Thetajj_sim;
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,1); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('simulated PCoh')
    %%
    X = Thetajj_est(:,:,1);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,2); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('higss-lasso PCoh')
    %%
    X = Thetajj_est(:,:,2);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,3); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('higgs-ridge PCoh')
    %%
    X = Thetajj_est(:,:,3);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,4); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('higgs-naive PCoh')    
    %%
    X = Thetajj_est(:,:,4);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,5); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('eloreta-hglasso PCoh')    
    %%
    X = Thetajj_est(:,:,5);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(2,3,6); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('lcmv-hglasso PCoh')
    colormap(cmap);
    %%
    
    saveas( figure_partial_coherence_maps,strcat(output_sourse,filesep,'partial_coherence_maps_(',sens_system,')_.fig'));
    disp(strcat('Saving figure ---->  Partial Coherence Maps to  ---> ', output_sourse) );
    delete(figure_partial_coherence_maps);
    
end

end