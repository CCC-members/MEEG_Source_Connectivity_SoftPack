function result = Main_SSVEP_AnalyzeData(output_sourse)

%% General Parameters

output_sourse = strcat(output_sourse, filesep,'SSVEP_AnalyzeData');
mkdir(output_sourse);

process_waitbar = waitbar(0,'Please wait...');
ntask_blocks    = 3;
eeg_task_points = {7454:18579, 29811:40936, 52168:63294};
load('colormap2.mat')                % Colormap
fs              = 250;               % Sampling Frequency
deltaf          = 0.22;              % Frequency step
%%
%% Subjects and Conditions
name_list       = {'CAGDAS_KARSAN'}; % Subjects
freq_list       = {4};               % Experimental conditions 'frequency'
%%
%% Load head model

addpath(strcat('ssvep',filesep,'Data',filesep,'CK',filesep,'EEG'));
load('ck_leadfield.mat');
L(end,:)        = [];
load('ck_model.mat');
coor(end,:)     = [];
%%
%% Split cortical surface into Left and Right Hemispheres for plotting maps
labels_conn  = {'OL-L' 'TL-L' 'FL-L' 'OL-R' 'TL-R' 'FL-R'};
TL_slope     = -0.6992;
TL_intercept = -2.1651;
[qL,qR,qfull,indvL,indvR,indv,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos] = split_hemispheres(cortex,coor);
OL_L = find((verticesL(:,2) < -40.72).*(verticesL(:,3) < 10.07));
TL_L = find((verticesL(:,3) < TL_slope*verticesL(:,2) + TL_intercept).*(verticesL(:,2) > -40.72));
FL_L = find((verticesL(:,2) > 38).*(verticesL(:,3) < 10.07));
OL_R = find((verticesR(:,2) < -40.72).*(verticesR(:,3) < 10.07));
TL_R = find((verticesR(:,3) < TL_slope*verticesR(:,2) + TL_intercept).*(verticesR(:,2) > -40.72));
FL_R = find((verticesR(:,2) > 38).*(verticesR(:,3) < 10.07));
%%
%% Load Subject Data
name            = name_list{1};
load([name,'.mat'])
%%
%% Cycle by conditions
for cond = 1:length(freq_list)
    waitbar((cond)/(length(freq_list)),process_waitbar,strcat('Collecting Task data'));
    freq                              = freq_list{cond};
    var                               = [name,'_',num2str(freq),'Hz_ICAcorrected','.data'];
    %%
    %% Collecting Task data
    data_EEG                          = eval(var);
    data_EEG(31:end,:)                = [];
    data_EEG_task                     = [];
    for count_blocks = 1:ntask_blocks
        data_EEG_task                 = [data_EEG_task, data_EEG(:,eeg_task_points{count_blocks})];
    end
    %%
    %% Cross-spectra
    fmax                              = 2*freq_list{cond} +5;
    [Svv_full,F,Ns,PSD]               = xspectrum(data_EEG_task,fs,fmax,deltaf);
    %%
    %% Applying average reference
    Nf = length(F);
    for cont_freq = 1:Nf
        [Svv_full(:,:,cont_freq),K]   = applying_reference(Svv_full(:,:,cont_freq),full(L));    % applying average reference...
    end
    figure_occ_harmonics(F,PSD,freq,output_sourse);
    %%
    %% split Lead Field
    K_L  = K(:,indvL);
    K_R  = K(:,indvR);
    %%
    %% Pick frequency bins for analysis
    Svv        = zeros(size(Svv_full,1));
    ratio_list = [1/2 1 2];
    neigh_list = [-2 -1 0 1 2];
    for ratio = ratio_list
        [tmp,freqid]                  = min(abs(F - ratio*freq));
        for neigh = neigh_list
            Svv_tmp                   = Svv_full(:,:,(freqid + neigh));
            Svv_tmp                   = Svv_tmp/(sum(abs(diag(Svv_tmp)))/size(Svv_tmp,1));
            Svv                       = Svv + Svv_tmp;
        end
    end
    Svv = Svv/(length(neigh_list)*length(ratio_list));
    Ns  = Ns*length(neigh_list)*length(ratio_list);
    %%
    %% Activation Leakage Module
    disp('activation leakage module...');
    % Default Atlas (groups)
    nonovgroups = [];
    for ii = 1:length(K)
        nonovgroups{ii} = ii;
    end
    [miu,sigma_post,DSTF]             = cross_nonovgrouped_enet_ssbl({Svv},{K},Ns,nonovgroups);
    stat                              = sqrt(abs(miu))./abs(sigma_post);
    indms                             = find(stat > 0.5);
    J                                 = zeros(qfull,1);
    J(indms)                          = stat(indms); J = J/max(abs(J));
    %%
    ngen                              = 12;
    JL                                = J(indvL);
    indmsL                            = best_ranked(J(indvL),{OL_L TL_L FL_L},ngen);
    JR                                = J(indvR);
    indmsR                            = best_ranked(J(indvR),{OL_R TL_R FL_R},ngen);
    %%
    figure_hemisphere =  figure;
    subplot(1,2,1);
    Jtmp = zeros(length(indvL),1);
    Jtmp([indmsL{1};indmsL{2};indmsL{3}]) = JL([indmsL{1};indmsL{2};indmsL{3}]);
    patch('Faces',facesL,'Vertices',verticesL,'FaceVertexCData',JL,'FaceColor','interp',...
        'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    axis([-90 5 -110 110 -65 95]); view([-90 0]); axis off;
    colormap(cmap);
    caxis([0 1]);
    title(['left hemisphere',' ',num2str(freq),'Hz'],'color','k')
    %%
    subplot(1,2,2);
    Jtmp = zeros(length(indvR),1);
    Jtmp([indmsR{1};indmsR{2};indmsR{3}]) = JR([indmsR{1};indmsR{2};indmsR{3}]);
    patch('Faces',facesR,'Vertices',verticesR,'FaceVertexCData',JR,'FaceColor','interp',...
        'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    axis([-5 90 -110 110 -65 95]); view([90 0]); axis off;
    colormap(cmap);
    caxis([0 1]);
    title(['right hemisphere',' ',num2str(freq),'Hz'],'color','k')
    
    
    saveas( figure_hemisphere,strcat(output_sourse,filesep,'left_&_right_hemisphere.fig'));
    disp(strcat('Saving figure ---->  left_&_right_hemisphere to  ---> ', output_sourse) );
    delete(figure_hemisphere);
    
    %%
    %% Connectivity Leakage Module
    disp('connectivity leakage module...');
    % parameters
    param.maxiter_outer               = 60;
    param.maxiter_inner               = 30;
    param.m                           = Ns;
    param.rth                         = 3.55;
    param.axi                         = 1E-1;
    param.Axixi                       = eye(size(Svv,1));
    %%
    Kindms                            = cat(2,K_L(:,[indmsL{1};indmsL{2};indmsL{3}]),K_R(:,[indmsR{1};indmsR{2};indmsR{3}]));
    q                                 = size(Kindms,2);
    penalty                           = [1,2,0];
    Thetajj                           = zeros(q,q,length(penalty));
    Sjj                               = zeros(q,q,length(penalty));
    for k_penalty = 1:length(penalty)
        param.penalty  = penalty(k_penalty);
        [Thetajj(:,:,k_penalty),Sjj(:,:,k_penalty),llh{k_penalty},jj_on,xixi_on] = h_hggm(Svv,Kindms,param);
    end
    
    
    %%
    figure_partial_correlation_maps1 = figure;
    set(gcf,'Position',[300 300 1400 300]);
    X = Thetajj(:,:,1);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,1); imagesc(abs(X)); axis off
    assign_labels(labels_conn,ngen,5,3,2)
    title('h-hggm-lasso partial correlations')
    X = Thetajj(:,:,2);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,2); imagesc(abs(X)); axis off
    assign_labels(labels_conn,ngen,5,3,2)
    title('h-hggm-ridge partial correlations')
    X = Thetajj(:,:,3);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,3); imagesc(abs(X)); axis off
    assign_labels(labels_conn,ngen,5,3,2)
    title('h-hggm-naive partial correlations')
    colormap('hot');
    
    saveas( figure_partial_correlation_maps1,strcat(output_sourse,filesep,'partial_correlation_maps_freq_',num2str(freq),'.fig'));
    disp(strcat('Saving figure ---->  Partial Correlation Maps to  ---> ', output_sourse) );
    delete(figure_partial_correlation_maps1);
    
    %%
    conn = zeros(length(labels_conn),length(labels_conn),3);
    for area1 = 1:length(labels_conn)
        for area2 = 1:length(labels_conn)
            conn_tmp            = Thetajj(((area1-1)*ngen + 1):area1*ngen,((area2-1)*ngen + 1):area2*ngen,:);
            conn_tmp(:,:,1)     = conn_tmp(:,:,1) - diag(diag(conn_tmp(:,:,1)));
            conn_tmp(:,:,2)     = conn_tmp(:,:,2) - diag(diag(conn_tmp(:,:,2)));
            conn_tmp(:,:,3)     = conn_tmp(:,:,3) - diag(diag(conn_tmp(:,:,3)));
            conn_tmp            = reshape(conn_tmp,ngen^2,3);
            conn_tmp            = mean(abs(conn_tmp),1);
            conn(area1,area2,:) = conn_tmp;
        end
    end
    %%
    figure_partial_correlation_maps2 = figure;
    set(gcf,'Position',[300 300 1400 300]);
    X = conn(:,:,1);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,1); imagesc(abs(X)); axis off
    assign_labels(labels_conn,1,0.01,0.8,1)
    title('h-hggm-lasso partial correlations')
    X = conn(:,:,2);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,2); imagesc(abs(X)); axis off
    assign_labels(labels_conn,1,0.01,0.8,1)
    title('h-hggm-ridge partial correlations')
    X = conn(:,:,3);
    X = X - diag(diag(X));
    X = X/max(abs(X(:)));
    subplot(1,3,3); imagesc(abs(X)); axis off
    assign_labels(labels_conn,1,0.01,0.8,1)
    title('h-hggm-naive partial correlations')
    colormap('hot');
    
    saveas( figure_partial_correlation_maps2,strcat(output_sourse,filesep,'roi_partial_correlation_maps_freq_',num2str(freq),'.fig'));
    disp(strcat('Saving figure ---->  Partial Correlation Maps to  ---> ', output_sourse) );
    delete(figure_partial_correlation_maps2);
    
    figure_iterations_likelihood = figure;
    plot(llh{1}); hold on; plot(llh{2}); hold on; plot(llh{3});
    xlabel('iterations'); ylabel('likelihood');
    legend('h-hggm-lasso','h-hggm-ridge','h-hggm-naive')
    
    saveas( figure_iterations_likelihood,strcat(output_sourse,filesep,'iterations_likelihood_freq_',num2str(freq),'.fig'));
    disp(strcat('Saving figure ---->  iterations_likelihood  to  ---> ', output_sourse) );
    delete(figure_iterations_likelihood);
    
end
delete(process_waitbar);
end