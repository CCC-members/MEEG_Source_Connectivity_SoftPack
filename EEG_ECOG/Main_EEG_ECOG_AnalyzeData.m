function [] = Main_EEG_ECOG_AnalyzeData(output_sourse)



output_sourse = strcat(output_sourse, filesep,'EEG_ECOG');
if(~isfolder(output_sourse))
mkdir(output_sourse);
end
  process_waitbar = waitbar(0,'General Parameters...');
%%
%% General Parameters 
disp('General Parameters...');
load('colormap2.mat')                                   % Colormap
fs                    = 1000;                           % Sampling Frequency
deltaf                = 0.25;                           % Frequency step
cond                  = 1;                              % 1-awake 2-anesthesia
prep                  = 2;                              % 1-bandpass 2-notch 3-refremoval
peak_pos1             = 9;
peak_pos2             = 14;
%  h-hggm parameters 
param.maxiter_outer   = 60;
param.maxiter_inner   = 30;
param.rth             = 3.55;
%  Sensor array and Conditions
sens_list             = {'EEG','ECoG'};                       % sensor system
cond_list             = {'02','05_anesthesia'};               % Experimental conditions 'awake/anesthesia'
prep_list             = {'bandpass','notch','refremoval'};    % preprocessing
%  cortex
cortex                = load(strcat('EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Source_Space',filesep,'Cortex-mid_Su.mat'));
index                 = load(strcat('EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'001',filesep,'EEG_leadfield_reduced.mat'));
index                 = index.ind;
 waitbar(1/4,process_waitbar,strcat('Split Left and Right Hemispheres' ));
%  Split cortical surface into Left and Right Hemispheres
[qL,qR,qfull,indvL,indvR,verticesL,verticesR,vertices,facesL,facesR,faces,elec] = split_hemispheres_monkey(cortex,1);
%  Get Rois
labels_conn           = {'I-OL' 'S-OL' 'P-TL' 'A-TL' 'I-PL' 'S-PL' 'I-FL' 'S-FL'};
IOL_SOL_slope         = -0.603724928366762;
IOL_SOL_intercept     = -12.968716332378225;
OL_TL_slope           = 1.756395839190329;
OL_TL_intercept       = 17.264430700028115;
PTL_ATL_slope         = 3.958790955237656;
PTL_ATL_intercept     = 11.396522288878636;
TL_PL_slope           = -0.819298245614035;
TL_PL_intercept       = 4.725245614035087;
PL_FL_slope           = -1.552157738095238;
PL_FL_intercept       = 33.539742559523810;
IPL_SPL_slope         = 0.146970304975923;
IPL_SPL_intercept     = 8.097142857142856;
IFL_SFL_slope         = -0.182281708094328;
IFL_SFL_intercept     = 12.785691523263225;
IOL                   = find((vertices(:,3) < IOL_SOL_slope*vertices(:,2) + IOL_SOL_intercept).*(vertices(:,3) > OL_TL_slope*vertices(:,2) + OL_TL_intercept));
SOL                   = find((vertices(:,3) > IOL_SOL_slope*vertices(:,2) + IOL_SOL_intercept).*(vertices(:,3) > OL_TL_slope*vertices(:,2) + OL_TL_intercept).*(vertices(:,3) < TL_PL_slope*vertices(:,2) + TL_PL_intercept));
PTL                   = find((vertices(:,3) > PTL_ATL_slope*vertices(:,2) + PTL_ATL_intercept).*(vertices(:,3) < TL_PL_slope*vertices(:,2) + TL_PL_intercept).*(vertices(:,3) < OL_TL_slope*vertices(:,2) + OL_TL_intercept));
ATL                   = find((vertices(:,3) < PTL_ATL_slope*vertices(:,2) + PTL_ATL_intercept).*(vertices(:,3) < TL_PL_slope*vertices(:,2) + TL_PL_intercept));
IPL                   = find((vertices(:,3) < PL_FL_slope*vertices(:,2) + PL_FL_intercept).*(vertices(:,3) > TL_PL_slope*vertices(:,2) + TL_PL_intercept).*(vertices(:,3) < IPL_SPL_slope*vertices(:,2) + IPL_SPL_intercept));
SPL                   = find((vertices(:,3) < PL_FL_slope*vertices(:,2) + PL_FL_intercept).*(vertices(:,3) > TL_PL_slope*vertices(:,2) + TL_PL_intercept).*(vertices(:,3) > IPL_SPL_slope*vertices(:,2) + IPL_SPL_intercept));
IFL                   = find((vertices(:,3) < IFL_SFL_slope*vertices(:,2) + IFL_SFL_intercept).*(vertices(:,3) > PL_FL_slope*vertices(:,2) + PL_FL_intercept));
SFL                   = find((vertices(:,3) < IFL_SFL_slope*vertices(:,2) + IFL_SFL_intercept).*(vertices(:,3) > PL_FL_slope*vertices(:,2) + PL_FL_intercept));
rois                  = {IOL, SOL, PTL, ATL, IPL , SPL, IFL, SFL};
%%
 waitbar(2/4,process_waitbar,strcat('Collecting data (ECoG)...' ));
%% Collecting data (ECoG)
disp('Collecting data (ECoG)...');
sens                  = 2;                                    % 1-EEG 2-ECoG
path                  = ['EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'002',filesep,sens_list{sens},cond_list{cond},prep_list{prep},'.mat'];
load(path)
data                  = eval([sens_list{sens},prep_list{prep}]);
elec                  = load(['EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Electrodes',filesep,sens_list{sens},'-elecs_Su.mat']);
elec                  = elec.electrodes;
%  Cross-spectra
waitbar(3/4,process_waitbar,strcat('Cross-spectra...' ));
[Svv_full,F,Ns,psd]   = xspectrum(data,fs,80,deltaf,6,[peak_pos1,peak_pos2]);
%  Pick frequency bins for analysis 
waitbar(3/4,process_waitbar,strcat('Pick frequency bins for analysis ' ));
peak_pos              = find(F == peak_pos1):find(F == peak_pos2);
Svv                   = mean(Svv_full(:,:,peak_pos),3);
Ns                    = Ns*length(peak_pos);
delete(process_waitbar);
% plot spectra
plot_spectra_figure = figure; 
set(gcf,'Position',[50 50 1400 700]);
min_psd = min(log(psd(:)));
max_psd = max(log(psd(:)));
plot_peak = min_psd*ones(length(F),1);
plot_peak(peak_pos) = max_psd;
plot(F, log(psd)); xlabel('frequency'); ylabel('log-PSD'); title([sens_list{sens},' power-spectral-density']);
hold on
plot(F,plot_peak,'--b');

saveas(plot_spectra_figure,strcat(output_sourse,filesep,[sens_list{sens},' power-spectral-density']));
delete(plot_spectra_figure);
%  Load head model
 process_waitbar = waitbar(0,'Load head model...');

Kecog                 = load(['EEG_ECOG',filesep,'data',filesep,'5-Leadfields',filesep,sens_list{sens},'-leadfield.mat']);
Kecog                 = Kecog.LFc;
%  Applying average reference
[output1,Kecog]       = applying_reference(eye(size(Kecog,1)),Kecog);    % applying average reference...
%% Activation Leakage Module (ECoG)
waitbar(1/4,process_waitbar,strcat('Activation Leakage Module (ECoG)...' ));
disp('ECoG activation leakage module...');
% Default Atlas (groups)
nonovgroups           = [];
counter               = 1;
for ii = 1:3:length(Kecog)
    nonovgroups{counter} = [ii;ii+1;ii+1];
    counter = counter + 1;
end
waitbar(3/4,process_waitbar,strcat('ENET-SSBL...' ));
%  ENET-SSBL
[miu,sigma_post,DSTF]             = cross_nonovgrouped_enet_ssbl({Svv},{Kecog},Ns,nonovgroups);
miu_stat                          = sqrt(sum(reshape(miu,3,length(Kecog)/3),1))';
sigma_post_stat                   = sqrt(sum(reshape(sigma_post,3,length(Kecog)/3),1))';
stat                              = abs(miu_stat)./abs(sigma_post_stat);
indms                             = find(stat > 1);
%  Plotting activity map
J                                 = zeros(qfull,1);
indms_red                         = intersect(indms,index);
J(indms)                          = stat(indms); J = J/max(abs(J));
JL                                = J(indvL);
delete(process_waitbar);

ECoG_figure = figure;
patch('Faces',facesL,'Vertices',verticesL,'FaceVertexCData',JL,'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.99);
axis([-90 5 -35 35 -20 20]); view([-90 0]); axis off;
colormap(cmap);
caxis([0 1]);
title(['ECoG reconstructed sources',' ',num2str(peak_pos1),'-',num2str(peak_pos2),'Hz'],'color','k','FontSize',24)

saveas(ECoG_figure,strcat(output_sourse,filesep,['ECoG reconstructed sources',' ',num2str(peak_pos1),'-',num2str(peak_pos2),'Hz']));
delete(ECoG_figure);
%%
%% Connectivity Leakage Module (ECoG)
process_waitbar = waitbar(0,strcat('ECoG connectivity leakage module...' ));
disp('ECoG connectivity leakage module...');
%  Computing normals
[nrm]                             = normals(vertices, faces, 'vertex');
%  Project and reduce Lead Field 
waitbar(1/5,process_waitbar,strcat('Project and reduce Lead Field' ));
Kproj                             = zeros(size(Kecog,1),size(Kecog,2)/3);
for chann = 1:length(Svv)
    Kproj(chann,:) = sum(reshape(Kecog(chann,:),3,size(Kecog,2)/3)'.*nrm,2);
end

K_L                               = [];
q_rois                            = cell(1,length(labels_conn));
indms_redL                        = intersect(indms_red,indvL);
for roi = 1:length(labels_conn)
    indms_tmp                     = intersect(indms_redL,rois{roi});
    q_rois{roi}                   = length(indms_tmp);
    K_L                           = [K_L Kproj(:,indms_tmp)];
end
%  H-HGGM
waitbar(2/5,process_waitbar,strcat('H-HGGM' ));
param.penalty                     = 1;
param.m                           = Ns;
param.axi                         = 1E-5;
param.sigma2xi                    = 1E0;
param.Axixi                       = eye(size(Svv,1));
q                                 = size(K_L,2);
Thetajj                           = zeros(q,q,10);
Sjj                               = zeros(q,q,10);
[Thetajj(:,:,1),Sjj(:,:,1),llh{1},jj_on,xixi_on] = h_hggm(Svv,K_L,param);
param.penalty                     = 2;
[Thetajj(:,:,2),Sjj(:,:,2),llh{2},jj_on,xixi_on] = h_hggm(Svv,K_L,param);
param.penalty                     = 0;
[Thetajj(:,:,3),Sjj(:,:,3),llh{3},jj_on,xixi_on] = h_hggm(Svv,K_L,param);

% eLORETA + HGGM
waitbar(3/5,process_waitbar,strcat('eLORETA + HGGM' ));
param.gamma1                      = 0.05;
param.gamma2                      = 0.1;
param.delta_gamma                 = 0.001;
[Thetajj(:,:,4),Sjj(:,:,4),gamma_grid{1},gamma{1},gcv{1}] = eloreta_hggm(Svv,K_L,param);
% LCMV + HGGM
waitbar(4/5,process_waitbar,strcat('LCMV + HGGM' ));
param.gamma                       = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,5),Sjj(:,:,5)]       = lcmv_hggm(Svv,K_L,param);

delete(process_waitbar);
%%
%% Collecting data (EEG)
process_waitbar = waitbar(1/3, strcat('Collecting data (EEG)' ));
disp('Collecting data (EEG)...');
sens            = 1;                                    % 1-EEG 2-ECoG
path            = ['EEG_ECOG',filesep,'data',filesep,'cleaned',filesep,'002',filesep,sens_list{sens},cond_list{cond},prep_list{prep},'.mat'];
load(path)
data            = eval([sens_list{sens},prep_list{prep}]);
elec            = load(['EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Electrodes',filesep,sens_list{sens},'-elecs_Su.mat']);
elec            = elec.electrodes;
%  Cross-spectra
waitbar(2/3,process_waitbar,strcat('Cross-spectra' ));
[Svv_full,F,Ns,psd] = xspectrum(data,fs,80,deltaf,6);
%  Pick frequency bins for analysis 
waitbar(2/3,process_waitbar,strcat('Pick frequency bins for analysis ' ));
peak_pos              = find(F == peak_pos1):find(F == peak_pos2);
Svv                   = mean(Svv_full(:,:,peak_pos),3);
Ns                    = Ns*length(peak_pos);
delete(process_waitbar);
% plot spectra
plot_spectra_figure_EEG = figure; 
set(gcf,'Position',[50 50 1400 700]);
min_psd = min(log(psd(:)));
max_psd = max(log(psd(:)));
plot_peak = min_psd*ones(length(F),1);
plot_peak(peak_pos) = max_psd;
plot(F, log(psd)); xlabel('frequency'); ylabel('log-PSD'); title([sens_list{sens},' power-spectral-density']);
hold on
plot(F,plot_peak,'--b');

saveas(plot_spectra_figure_EEG,strcat(output_sourse,filesep,[sens_list{sens},' power-spectral-density']));
delete(plot_spectra_figure_EEG);
%  Load head model
process_waitbar1 = waitbar(0,strcat('Load head model' ));
Keeg            = load(['EEG_ECOG',filesep,'data',filesep,'5-Leadfields',filesep,sens_list{sens},'-leadfield.mat']);
Keeg            = Keeg.LF;
%  Applying average reference
 waitbar(1/5,process_waitbar1,strcat(' Applying average reference' ));
[output1,Keeg] = applying_reference(eye(size(Keeg,1)),Keeg);    % applying average reference...
Keeg([end-1 end],:) = [];
%%
%% Connectivity Leakage Module (EEG)
 waitbar(2/5,process_waitbar1,strcat(' Connectivity Leakage Module (EEG)' ));
disp('EEG connectivity leakage module...');
% Project and reduce Lead Field 
Kproj                             = zeros(size(Keeg,1),size(Keeg,2)/3);
for chann = 1:length(Svv)
    Kproj(chann,:) = sum(reshape(Keeg(chann,:),3,size(Keeg,2)/3)'.*nrm,2);
end
K_L                               = [];
for roi = 1:length(labels_conn)
    indms_tmp                     = intersect(indms_redL,rois{roi});
    K_L                           = [K_L Kproj(:,indms_tmp)];
end
% H-HGGM
 waitbar(3/5,process_waitbar1,strcat(' H-HGGM' ));
param.penalty                     = 1;
param.m                           = Ns;
param.axi                         = 1E-1;
param.sigma2xi                    = 1E3;
param.Axixi                       = eye(size(Svv,1));
[Thetajj(:,:,6),Sjj(:,:,6),llh{4},jj_on,xixi_on] = h_hggm(Svv,K_L,param);
param.penalty                     = 2;
[Thetajj(:,:,7),Sjj(:,:,7),llh{5},jj_on,xixi_on] = h_hggm(Svv,K_L,param);
param.penalty                     = 0;
[Thetajj(:,:,8),Sjj(:,:,8),llh{6},jj_on,xixi_on] = h_hggm(Svv,K_L,param);
% eLORETA + HGGM
 waitbar(4/5,process_waitbar1,strcat('eLORETA + HGGM' ));
param.gamma1                      = 0.1;
param.gamma2                      = 5;
param.delta_gamma                 = 0.1;
[Thetajj(:,:,9),Sjj(:,:,9),gamma_grid{2},gamma{2},gcv{2}] = eloreta_hggm(Svv,K_L,param);
% LCMV + HGGM
waitbar(5/5,process_waitbar1,strcat('LCMV + HGGM' ));
param.gamma                       = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj(:,:,10),Sjj(:,:,10)]     = lcmv_hggm(Svv,K_L,param);
delete(process_waitbar1);
%%
%% Plotting results
% EEG/ECoG sensor distribution

EEG_ECoG_sensor_distribution_figure = figure;
elec                  = load(['EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Electrodes',filesep,sens_list{1},'-elecs_Su.mat']);
elec                  = elec.electrodes;
[qL,qR,qfull,indvL,indvR,verticesL,verticesR,vertices,facesL,facesR,faces,elec] = split_hemispheres_monkey(cortex,elec);
patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',zeros(length(vertices),1),'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.99);
axis off;
colormap(cmap);
hold on
scatter3(elec(:,1),elec(:,2),elec(:,3),'MarkerEdgeColor','b','MarkerFaceColor','b')
elec                  = load(['EEG_ECOG',filesep,'data',filesep,'4-Head_Model',filesep,'Electrodes',filesep,sens_list{2},'-elecs_Su.mat']);
elec                  = elec.electrodes;
[qL,qR,qfull,indvL,indvR,verticesL,verticesR,vertices,facesL,facesR,faces,elec] = split_hemispheres_monkey(cortex,elec);
hold on
scatter3(elec(:,1),elec(:,2),elec(:,3),'MarkerEdgeColor','r','MarkerFaceColor','r')

saveas(EEG_ECoG_sensor_distribution_figure,strcat(output_sourse,filesep,'EEG-ECoG sensor distribution'));
delete(EEG_ECoG_sensor_distribution_figure);

% node-wise partial correlations
node_wise_partial_correlations  = figure;
set(gcf,'Position',[300 300 1400 300]);
X = Thetajj(:,:,1);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(1,2,1); imagesc(abs(X));
title('ECoG partial correlations')
X = Thetajj(:,:,6);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(1,2,2); imagesc(abs(X));
title('EEG partial correlations')
colormap('hot');

saveas(node_wise_partial_correlations,strcat(output_sourse,filesep,'EEG partial correlations1'));
delete(node_wise_partial_correlations);

% Roi analysis
gen1 = 1:q_rois{1};
conn = zeros(length(labels_conn),length(labels_conn),2);
for roi1 = 1:length(labels_conn)
    gen2 = 1:q_rois{1};
    for roi2 = 1:length(labels_conn)
        conn_tmp            = Thetajj(gen1,gen2,[1 6]);
        if roi1 == roi2
            conn_tmp(:,:,1)     = conn_tmp(:,:,1) - diag(diag(conn_tmp(:,:,1)));
            conn_tmp(:,:,2)     = conn_tmp(:,:,2) - diag(diag(conn_tmp(:,:,2)));
        end
        conn_tmp            = reshape(conn_tmp,q_rois{roi1}*q_rois{roi2},2);
        conn_tmp            = mean(abs(conn_tmp),1);
        conn(roi1,roi2,:)   = conn_tmp;
        if roi2 < length(labels_conn)
            gen2                = (gen2(end) + 1):(gen2(end) + q_rois{roi2 + 1});
        end
    end
    if roi1 < length(labels_conn)
        gen1                = (gen1(end) + 1):(gen1(end) + q_rois{roi1 + 1});
    end
end

node_wise_partial_correlations2 = figure;
set(gcf,'Position',[300 300 1400 300]);
X = conn(:,:,1);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(1,2,1); imagesc(abs(X)); axis off
assign_labels(labels_conn,1,0.01,0.8,1)
title('ECoG partial correlations')
X = conn(:,:,2);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(1,2,2); imagesc(abs(X)); axis off
assign_labels(labels_conn,1,0.01,0.8,1)
title('EEG partial correlations')
colormap('hot');

saveas(node_wise_partial_correlations2,strcat(output_sourse,filesep,'EEG partial correlations2'));
delete(node_wise_partial_correlations2);

% likelihood

likelihood_figure = figure;
subplot(2,4,1)
plot(llh{1}); hold on; plot(llh{4});
xlabel('iterations'); ylabel('likelihood');
legend('ECoG','EEG')
title('h-hggm-lasso likelihood')
subplot(2,4,2)
plot(llh{2}); hold on; plot(llh{5});
xlabel('iterations'); ylabel('likelihood');
legend('ECoG','EEG')
title('h-hggm-ridge likelihood')
subplot(2,4,3)
plot(llh{3}); hold on; plot(llh{5});
xlabel('iterations'); ylabel('likelihood');
legend('ECoG','EEG')
title('h-hggm-naive likelihood')
% gcv
subplot(2,4,4)
[gcv_opt,idx_gamma] = min(gcv{1});
plot(gamma_grid{1},gcv{1},...
        '-',gamma_grid{1}(idx_gamma),...
        gcv_opt,'b*');
xlabel('regularization parameter'); ylabel('gcv value');
legend('ECoG')
title('eloreta-hggm gcv')
subplot(2,4,5)
[gcv_opt,idx_gamma] = min(gcv{2});
plot(gamma_grid{2},gcv{2},...
        '-',gamma_grid{2}(idx_gamma),...
        gcv_opt,'b*');
xlabel('regularization parameter'); ylabel('gcv value');
legend('EEG')
title('eloreta-hggm gcv')

saveas(likelihood_figure,strcat(output_sourse,filesep,'likelihood'));
delete(likelihood_figure);
%%
%% distances
distance_label    = {'alphadiv' 'kullback' 'ld' 'logeuclid' 'opttransp' 'riemann'};
methods_label     = {'h-hggm-lasso';'h-hggm-ridge';'h-hggm-naive';'eloreta+hggm';'lcmv+hggm'};
Table = cell(6,7);
Table(2:end,1)    = methods_label;
Table(1,2:end)    = distance_label;
dist              = cell(5,6);
for cont = 1:size(Thetajj,3)/2
    A            = Thetajj(:,:,cont)/sum(abs(diag(Thetajj(:,:,cont)))) + eye(q);
    B            = Thetajj(:,:,cont+5)/sum(abs(diag(Thetajj(:,:,cont+5)))) + eye(q);
    dist{cont,1} = round(abs(distance_alphadiv(A,B,0.01)),4);
    dist{cont,2} = round(abs(distance_kullback(A,B)),4);
    dist{cont,3} = round(abs(distance_ld(A,B)),4);
    dist{cont,4} = round(abs(distance_logeuclid(A,B)),4);
    dist{cont,5} = round(abs(distance_opttransp(A,B)),4);
    dist{cont,6} = round(abs(distance_riemann(A,B)),4);
end
Table(2:end,2:end)  = dist;
save(strcat(output_sourse,filesep,'Table.mat'))

end


