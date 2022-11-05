function [] = Main_EEG_ECOG_AnalyzeData(output_source)

output_source = strcat(output_source, filesep,'EEG_ECOG');
if(~isfolder(output_source))
mkdir(output_source);
end
import +ROInets.*
                      process_waitbar = waitbar(0,'General Parameters...');
%%
%% General Parameters 
disp('General Parameters...');
load('colormap2.mat')                                   % Colormap
fs                    = 1000;                           % Sampling Frequency
Ntrials               = 30;                             % Number of transients
deltaf                = 0.25;                           % Frequency step
varf                  = 0.5;                            % filter width
peak_pos1             = 8;                              % starting frequency point of the band 
peak_pos2             = 14;                             % final frequency point of the band 
%  higgs parameters
param.maxiter_outer   = 60;
param.maxiter_inner   = 30;
param.rth1            = 0.7;
param.rth2            = 3.16;
%  Sensor array and Conditions
sens_list             = {'EEG','ECoG'};                       % sensor system
cond_list             = {'02','05_anesthesia'};               % Experimental conditions 'awake/anesthesia'
prep_list             = {'bandpass','notch','refremoval','cross_spectrum'};    % preprocessing
distance_label        = {'Kullback' 'Logeuclid' 'Riemann','Levenshtein','Opttransp','Diversity'};
methods_label         = {'onestep hgLASSO';'onestep hgRidge';'onestep hgNaive';'multistep eLORETA hgLASSO';'multistep LCMV hgLASSO'};
Dist_trials           = zeros(5,6,Ntrials);
Dist_stat             = cell(6,7,2);
Dist_global           = cell(6,7);
Dist_stat(2:end,1,1)  = methods_label;
Dist_stat(1,2:end,1)  = distance_label;
Dist_global(2:end,1)  = methods_label;
Dist_global(1,2:end)  = distance_label;
%  cortex
cortex                = load(strcat('EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Source_Space',filesep,'Cortex-mid_Su.mat'));
index                 = load(strcat('EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'001',filesep,'EEG_leadfield_reduced.mat'));
index                 = index.ind;
load(strcat('EEG_ECOG',filesep,'data',filesep,'5-Leadfields',filesep,'structural.mat'));
[nrm]                 = normals(vertices, faces, 'vertex');
                      
%% Collecting data
% ECoG
disp('Collecting data (ECoG)...');
sens                  = 2;                              % 1-EEG 2-ECoG
cond                  = 1;                              % 1-awake 2-anesthesia
prep                  = 2;                              % 1-bandpass 2-notch 3-refremoval
path                  = ['EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'002',filesep,sens_list{sens},cond_list{cond},prep_list{prep},'.mat'];
load(path);
ECoG                  = ECoGnotch;
ECoG                  = reshape(ECoG,128,size(ECoG,2)/Ntrials,Ntrials);
% EEG
disp('Collecting data (EEG)...');
sens                  = 1;                              % 1-EEG 2-ECoG
cond                  = 1;                              % 1-awake 2-anesthesia
prep                  = 2;                              % 1-bandpass 2-notch 3-refremoval
path                  = ['EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'002',filesep,sens_list{sens},cond_list{cond},prep_list{prep},'.mat'];
load(path);
EEG                   = EEGnotch;
EEG                   = reshape(EEG,18,size(EEG,2)/Ntrials,Ntrials);

%% Collecting Lead Fields
% ECoG
sens                  = 2;                              % 1-EEG 2-ECoG
Kecog                 = load(['EEG_ECOG',filesep,'data',filesep,'5-Leadfields',filesep,sens_list{sens},'-leadfield.mat']);
Kecog                 = Kecog.LFc;
Kecog_proj            = zeros(size(Kecog,1),size(Kecog,2)/3);
for chann = 1:size(ECoG,1)
    Kecog_proj(chann,:) = sum(reshape(Kecog(chann,:),3,size(Kecog,2)/3)'.*nrm,2);
end
% EEG
sens                  = 1;                              % 1-EEG 2-ECoG
Keeg                  = load(['EEG_ECOG',filesep,'data',filesep,'5-Leadfields',filesep,sens_list{sens},'-leadfield.mat']);
Keeg                  = Keeg.LF;
Keeg([end-1 end],:)   = [];
Keeg_proj             = zeros(size(Keeg,1),size(Keeg,2)/3);
for chann = 1:size(EEG,1)
    Keeg_proj(chann,:) = sum(reshape(Keeg(chann,:),3,size(Keeg,2)/3)'.*nrm,2);
end

%% Network screening and Lead Field reduction
% Cross-spectra (ECoG) all trials
data                  = reshape(ECoG,128,(size(ECoG,2))*Ntrials);       % all trials
Nw                    = 1;                              % Number of slepian windows
[Svv,Ns]              = xspectrum_band(data,fs,peak_pos1,peak_pos2,deltaf,varf,Nw);
% ENET-SSBL
parcels           = [];
counter               = 1;
for ii = 1:3:length(Kecog)
    parcels{counter} = [ii;ii+1;ii+1];
    counter = counter + 1;
end
[s2j,sigma2j_post]    = sSSBLpp_ultralarge(Svv,Kecog,parcels);
% stat
miu_stat              = sqrt(sum(reshape(s2j,3,length(Kecog)/3),1))';
sigma_post_stat       = sqrt(sum(reshape(sigma2j_post,3,length(Kecog)/3),1))';
stat                  = abs(miu_stat)./abs(sigma_post_stat);
indms                 = find(stat > 1);
J                     = zeros(qfull,1);
J(indms)              = stat(indms); J = J/max(abs(J));
JL                    = J(indvL);
% figure stat
ECoG_figure = figure;
patch('Faces',facesL,'Vertices',verticesL,'FaceVertexCData',JL,'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.99);
axis([-90 5 -35 35 -20 20]); view([-90 0]); axis off;
colormap(cmap);
caxis([0 1]);
title(['ECoG reconstructed sources',' ',num2str(peak_pos1),'-',num2str(peak_pos2),'Hz'],'color','k','FontSize',24)
saveas(ECoG_figure,strcat(output_source,filesep,['ECoG reconstructed sources',' ',num2str(peak_pos1),'-',num2str(peak_pos2),'Hz']));
delete(ECoG_figure);
% Reduction
indms_red             = intersect(indms,index);
q_rois                = cell(1,length(labels_conn));
indms_redL            = intersect(indms_red,indvL);
% Lead Field reduction ECoG
K_L_ecog              = [];
for roi = 1:length(labels_conn)
    indms_tmp         = intersect(indms_redL,rois{roi});
    q_rois{roi}       = length(indms_tmp);
    K_L_ecog               = [K_L_ecog Kecog_proj(:,indms_tmp)];
end
% Lead Field reduction EEG
K_L_eeg              = [];
for roi = 1:length(labels_conn)
    indms_tmp         = intersect(indms_redL,rois{roi});
    q_rois{roi}       = length(indms_tmp);
    K_L_eeg               = [K_L_eeg Keeg_proj(:,indms_tmp)];
end

%% ECoG
%  Cross-spectra (ECoG) trial
data                  = reshape(ECoG,128,(size(ECoG,2))*Ntrials);       % all trials
Nw                    = 1;                              % Number of slepian windows
[Svv,Ns]              = xspectrum_band(data,fs,peak_pos1,peak_pos2,deltaf,varf,Nw);
% Connectivity Leakage Module (ECoG)
disp('ECoG connectivity leakage module...');
%  Higgs
[p,q]                 = size(K_L_ecog);
param.use_gpu         = 1;
param.run_bash_mode   = 1;
param.str_band        = "none";
param.p               = p;
param.q               = q;
param.Ip              = eye(p);
param.Op              = ones(p,1);
param.Iq              = eye(q);
param.m               = Ns;
aj                    = sqrt(log(q)/Ns);                                           
Ajj_diag              = 0;                                                      
Ajj_ndiag             = 1;                                               
Ajj                   = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj              = aj;
param.Ajj             = Ajj;
param.axi             = 1E-1;
param.Axixi           = eye(p);
param.Axixi_inv       = eye(p);
param.ntry            = 0;
param.prew            = 1;
param.penalty         = 1;
param.nu              = Ns;
param.eigreg          = 1E-4;
[Thetajj(:,:,1),~,~] = higgs(Svv,K_L_ecog,param);
param.penalty                      = 2;
[Thetajj(:,:,2),~,~] = higgs(Svv,K_L_ecog,param);
param.penalty                      = 0;
[Thetajj(:,:,3),~,~] = higgs(Svv,K_L_ecog,param);
% eLORETA + HGGM
param.gamma1          = 0.01;
param.gamma2          = 0.5;
param.delta_gamma     = 0.01;
[Thetajj(:,:,4),~,~,~,~] = eloreta_hg_lasso(Svv,K_L_ecog,param);
% LCMV + HGGM
param.gamma           = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,5),~] = lcmv_hg_lasso(Svv,K_L_ecog,param);

%% EEG
%  Cross-spectra (EEG) trial
data                  = reshape(EEG,18,(size(EEG,2))*Ntrials);       % all trials
Nw                    = 1;                              % Number of slepian windows
[Svv,Ns]              = xspectrum_band(data,fs,peak_pos1,peak_pos2,deltaf,varf,Nw);
% Connectivity Leakage Module (EEG)
disp('EEG connectivity leakage module...');
% Higgs
[p,q]                 = size(K_L_eeg);
param.use_gpu         = 1;
param.run_bash_mode   = 1;
param.str_band        = "none";
param.p               = p;
param.q               = q;
param.Ip              = eye(p);
param.Op              = ones(p,1);
param.Iq              = eye(q);
param.m               = Ns;
aj                    = sqrt(log(q)/Ns);                                           
Ajj_diag              = 0;                                                      
Ajj_ndiag             = 1;                                               
Ajj                   = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj              = aj;
param.Ajj             = Ajj;
param.axi             = 1E-1;
param.Axixi           = eye(p);
param.Axixi_inv       = eye(p);
param.ntry            = 0;
param.prew            = 1;
param.penalty         = 1;
param.nu              = Ns;
param.eigreg          = 1E-4;
[Thetajj(:,:,6),~,~] = higgs(Svv,K_L_eeg,param);
param.penalty         = 2;
[Thetajj(:,:,7),~,~] = higgs(Svv,K_L_eeg,param);
param.penalty         = 0;
[Thetajj(:,:,8),~,~] = higgs(Svv,K_L_eeg,param);
% eLORETA + HGGM
param.gamma1          = 0.01;
param.gamma2          = 1;
param.delta_gamma     = 0.01;
[Thetajj(:,:,9),Sjj(:,:,9),gamma_grid{2},gamma{2},gcv{2}] = eloreta_hg_lasso(Svv,K_L_eeg,param);
% LCMV + HGGM
param.gamma           = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,10),Sjj(:,:,10)]     = lcmv_hg_lasso(Svv,K_L_eeg,param);

%% precision-matrix
node_wise_partial_correlations  = figure;
set(gcf,'Position',[300 300 1400 300]);
methods_label = {'onestep hgLASSO';'onestep hgRidge';'onestep hgNaive';'multistep eLORETA hgLASSO';'multistep LCMV hgLASSO';... % ECoG
    'onestep hgLASSO';'onestep hgRidge';'onestep hgNaive';'multistep eLORETA hgLASSO';'multistep LCMV hgLASSO'}; % EEG
Thetajj_norm = zeros(size(Thetajj,2),size(Thetajj,2),length(methods_label));
for meth = 1:length(methods_label)
    X                 = Thetajj(:,:,meth);
%     X                 = (X - diag(diag(X)))./sqrt((abs(diag(X))*abs(diag(X)')));
    X                 = (X - diag(diag(X)));
    X                 = (X + X')/2;
    Thetajj_norm(:,:,meth) = X;
    subplot(2,5,meth); imagesc(abs(X));
    title(methods_label{meth});
end
colormap('hot');
saveas(node_wise_partial_correlations,strcat(output_source,filesep,'node_wise_partial_coherences'));
delete(node_wise_partial_correlations);

%% divergence
node_wise_kullback  = figure;
set(gcf,'Position',[300 300 1400 300]);
methods_label = {'higgs-lasso';'higgs-ridge';'higgs-naive';'eloreta+hglasso';'lcmv+hglasso';... % ECoG
    'higgs-lasso';'higgs-ridge';'higgs-naive';'eloreta+hglasso';'lcmv+hglasso'}; % EEG
for meth = 1:(length(methods_label)/2)
    A            = Thetajj(:,:,meth);
%     normA        = sqrt(sum(abs(diag(A*A')))/q);
%     A            = A/normA;
    A            = (A + A')/2;
    B            = Thetajj(:,:,meth + length(methods_label)/2);
%     normB        = sqrt(sum(abs(diag(B*B')))/q);
%     B            = B/normB;
    B            = (B + B')/2;
    X            = map_kullback(A,B);
    X            = (X + X')/2;
    X            = X - diag(diag(X));
    subplot(1,5,meth); imagesc(abs(X));
    caxis([0 6000])
    title(methods_label{meth});
end
colormap('hot');
saveas(node_wise_kullback,strcat(output_source,filesep,'node_wise_kullback'));
delete(node_wise_kullback);

%% distances
dist              = zeros(5,6);
for cont = 1:size(Thetajj,3)/2
    A            = Thetajj(:,:,cont);
    A            = A/sqrt(sum(abs(diag(A*A')))/q);
    B            = Thetajj(:,:,cont+5);
    B            = B/sqrt(sum(abs(diag(B*B')))/q);
    dist(cont,1) = round(abs(distance_kullback(A,B)),4);
    dist(cont,2) = round(abs(distance_logeuclid(A,B)),4);
    dist(cont,3) = round(abs(distance_riemann(A,B)),4);
    dist(cont,4) = round(abs(distance_ld(A,B)),4);
    dist(cont,5) = round(abs(distance_opttransp(A,B)),4);
    dist(cont,6) = round(abs(distance_alphadiv(A,B,0.5)),4);
end
Dist_global(2:end,2:end)  = num2cell(dist);
%% Loop for trials
for trial = 1:Ntrials
%% ECoG
%  Cross-spectra (ECoG) trial
data                  = squeeze(ECoG(:,:,trial));       % trial
Nw                    = 1;                              % Number of slepian windows
[Svv,Ns]              = xspectrum_band(data,fs,peak_pos1,peak_pos2,deltaf,varf,Nw);
% Connectivity Leakage Module (ECoG)
disp('ECoG connectivity leakage module...');
%  Higgs
[p,q]                 = size(K_L_ecog);
param.use_gpu         = 1;
param.run_bash_mode   = 1;
param.str_band        = "none";
param.p               = p;
param.q               = q;
param.Ip              = eye(p);
param.Op              = ones(p,1);
param.Iq              = eye(q);
param.m               = Ns;
aj                    = sqrt(log(q)/Ns);                                           
Ajj_diag              = 0;                                                      
Ajj_ndiag             = 1;                                               
Ajj                   = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj              = aj;
param.Ajj             = Ajj;
param.axi             = 1E-1;
param.Axixi           = eye(p);
param.Axixi_inv       = eye(p);
param.ntry            = 0;
param.prew            = 1;
param.penalty         = 1;
param.nu              = Ns;
param.eigreg          = 1E-4;
[Thetajj(:,:,1),~,~] = higgs(Svv,K_L_ecog,param);
param.penalty                      = 2;
[Thetajj(:,:,2),~,~] = higgs(Svv,K_L_ecog,param);
param.penalty                      = 0;
[Thetajj(:,:,3),~,~] = higgs(Svv,K_L_ecog,param);
% eLORETA + HGGM
param.gamma1          = 0.01;
param.gamma2          = 0.5;
param.delta_gamma     = 0.01;
[Thetajj(:,:,4),~,~,~,~] = eloreta_hg_lasso(Svv,K_L_ecog,param);
% LCMV + HGGM
param.gamma           = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,5),~] = lcmv_hg_lasso(Svv,K_L_ecog,param);

%% EEG
%  Cross-spectra (EEG) trial
data                  = squeeze(EEG(:,:,trial));       % trial
Nw                    = 1;                              % Number of slepian windows
[Svv,Ns]              = xspectrum_band(data,fs,peak_pos1,peak_pos2,deltaf,varf,Nw);
% Connectivity Leakage Module (EEG)
disp('EEG connectivity leakage module...');
% Higgs
[p,q]                 = size(K_L_eeg);
param.use_gpu         = 1;
param.run_bash_mode   = 1;
param.str_band        = "none";
param.p               = p;
param.q               = q;
param.Ip              = eye(p);
param.Op              = ones(p,1);
param.Iq              = eye(q);
param.m               = Ns;
aj                    = sqrt(log(q)/Ns);                                           
Ajj_diag              = 0;                                                      
Ajj_ndiag             = 1;                                               
Ajj                   = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj              = aj;
param.Ajj             = Ajj;
param.axi             = 1E-1;
param.Axixi           = eye(p);
param.Axixi_inv       = eye(p);
param.ntry            = 0;
param.prew            = 1;
param.penalty         = 1;
param.nu              = Ns;
param.eigreg          = 1E-4;
[Thetajj(:,:,6),~,~] = higgs(Svv,K_L_eeg,param);
param.penalty         = 2;
[Thetajj(:,:,7),~,~] = higgs(Svv,K_L_eeg,param);
param.penalty         = 0;
[Thetajj(:,:,8),~,~] = higgs(Svv,K_L_eeg,param);
% eLORETA + HGGM
param.gamma1          = 0.01;
param.gamma2          = 1;
param.delta_gamma     = 0.01;
[Thetajj(:,:,9),Sjj(:,:,9),gamma_grid{2},gamma{2},gcv{2}] = eloreta_hg_lasso(Svv,K_L_eeg,param);
% LCMV + HGGM
param.gamma           = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,10),Sjj(:,:,10)]     = lcmv_hg_lasso(Svv,K_L_eeg,param);

%% distances
dist              = zeros(5,6);
for cont = 1:size(Thetajj,3)/2
    A            = Thetajj(:,:,cont);
    A            = A/sqrt(sum(abs(diag(A*A')))/q);
    B            = Thetajj(:,:,cont+5);
    B            = B/sqrt(sum(abs(diag(B*B')))/q);
    dist(cont,1) = round(abs(distance_kullback(A,B)),4);
    dist(cont,2) = round(abs(distance_logeuclid(A,B)),4);
    dist(cont,3) = round(abs(distance_riemann(A,B)),4);
    dist(cont,4) = round(abs(distance_ld(A,B)),4);
    dist(cont,5) = round(abs(distance_opttransp(A,B)),4);
    dist(cont,6) = round(abs(distance_alphadiv(A,B,0.5)),4);
end
Dist_trials(:,:,trial)  = dist;
end
Dist_stat(2:end,2:end,1)  = num2cell(squeeze(mean(Dist_trials,3)));
Dist_stat(2:end,2:end,2)  = num2cell(squeeze(std(Dist_trials,[],3)));
Table.Dist_trials         = Dist_trials;
Table.Dist_stat           = Dist_stat;
Table.Dist_global         = Dist_global;
save(strcat(output_source,filesep,'Table.mat'),'Table')

end