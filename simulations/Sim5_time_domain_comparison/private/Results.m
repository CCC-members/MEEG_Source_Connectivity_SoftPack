function [] = Results(output_source)

%% Run Results
%% Loading Simulation Substrate or Real Data
load('colormap2.mat')
% load('Pseudorand_Net.mat')
load([output_source,filesep,'Pseudorand_Net.mat'])
measures_label    = {'auc' 'sens' 'spec' 'prec' 'f1'};
methods_label     = {'higgs-lasso';'higgs-ridge';'higgs-naive';'eloreta-hglasso';'lcmv-hglasso'};
Nmeasures         = length(measures_label);
%%
%% h_hggm
% load('Solutions_higgs.mat')
load([output_source,filesep,'Solutions_higgs.mat'])
alphaPeakSliceNum=round(options.nu0(2)/options.deltaf)+1;
[measures]        = quality_measures(sol_higgs,theta_sim,alphaPeakSliceNum);
measures(isnan(measures)) = 0.5;
measures_mean     = round(100*mean(measures,3));
measures_std      = round(100*std(measures,0,3));
%%
%%
process_waitbar = waitbar(0,'Please wait...');
for meth = 1:5
    for meas = 1:5
        waitbar((meth*meas)/(25),process_waitbar,strcat('Outputing results....'));
        measures_higgs{meth,meas}  = [num2str(measures_mean(meth,meas)) '+/-' num2str(measures_std(meth,meas))];
    end
end
delete(process_waitbar);
%%
%% Table with quality measures
Table = cell(6,6);
Table(2:end,1)    = methods_label;
Table(1,2:end)    = measures_label;
Table(2:end,2:end)  = measures_higgs;


save(strcat(output_source,filesep,'Table_sens_system_',sens_system,'.mat'));
disp(strcat('Saving Table ---->  Table with quality measures to  ---> ', output_source) );

    
%%
%% Plot likelihood


Nsim = ceil(size(sol_higgs,2));
figure_likelihood = figure; 

subplot(3,2,1);
llh = zeros(length(sol_higgs{4,1}{1}{1}),Nsim);

for sim = 1:Nsim
  
    llh(:,sim) = sol_higgs{4,sim}{1}{1};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-lasso likelihood')


%%
%%

subplot(3,2,2); 
llh = zeros(length(sol_higgs{4,1}{1}{2}),Nsim);
for sim = 1:Nsim
   
    llh(:,sim) = sol_higgs{4,sim}{1}{2};  
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-ridge likelihood')


%%

subplot(3,2,3); 
llh = zeros(length(sol_higgs{4,1}{1}{3}),Nsim);
for sim = 1:Nsim
    
    llh(:,sim) = sol_higgs{4,sim}{1}{3};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-naive likelihood')


%%


subplot(3,2,4); 
for sim = 1:Nsim
    [gcv_opt,idx_gamma]       = min(sol_higgs{4,sim}{3});
    plot(sol_higgs{4,1}{2},sol_higgs{4,sim}{3},...
        '-',sol_higgs{4,1}{2}(idx_gamma),...
        gcv_opt,'b*');
    hold on;
end

ylabel('gcv value')
xlabel('regularization parameter')
title('eloreta-hglasso gcv function')

%% Plot corticaL map 
cortex.vertices = options.vertices;
cortex.faces    = options.faces;
[qL,qR,qfull,indvL,indvR,indv,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos_trans] = split_hemispheres(cortex,options.elec_pos);
J = zeros(qfull,1);
J(index_full) = 1;
%%
if strcmp(sens_system,'large') == 1 || strcmp(sens_system,'small') == 1
    subplot(3,2,5); 
    patch('Faces',facesL,'Vertices',verticesL,'FaceVertexCData',J(indvL),'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    if strcmp(sens_system,'large') == 1
        axis([-0.1 0.1 -0.005 0.1 -0.1 0.1]); axis off; view([-179.6 24.8]); % large
    elseif strcmp(sens_system,'small') == 1
        axis([-90 5 -110 110 -65 95]); axis off; view([-88.8 20]); % small
    end
    hold on 
    scatter3(elec_pos_trans(:,1),elec_pos_trans(:,2),elec_pos_trans(:,3),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('realistic-head model','color','k')
    axis equal 
    %%
    subplot(3,2,6);
    patch('Faces',facesR,'Vertices',verticesR,'FaceVertexCData',J(indvR),'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    if strcmp(sens_system,'large') == 1
        axis([-0.1 0.1 -0.1 0.005 -0.1 0.1]); axis off; view([7.60000000000002 27.2]); % large
    elseif strcmp(sens_system,'small') == 1
        axis([-5 90 -110 110 -65 95]); axis off; view([85.2 26.4]); % small
    end
    hold on 
    scatter3(elec_pos_trans(:,1),elec_pos_trans(:,2),elec_pos_trans(:,3),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('realistic-head model','color','k')
    axis equal
else
    subplot(3,2,5); 
    patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',J,'FaceColor','interp','Marker','o','MarkerFaceColor','y','EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    hold on 
    scatter(elec_pos(:,1),elec_pos(:,2),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    axis off;
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('pseudo-head model','color','k')
    axis equal
end

saveas( figure_likelihood,strcat(output_source,filesep,'higgs likelihood_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  higgs likelihood to  ---> ', output_source) );
delete(figure_likelihood);


%%
figure_partial_coherences = figure;
load('colormap3')
%% Plot partial correlations
X  = theta_sim{1}(:,:,alphaPeakSliceNum);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,1); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('simulated PCoh')
%%
X  = sol_higgs{3,1}(:,:,1);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,2); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('higgs-lasso PCoh')
%%
X  = sol_higgs{3,1}(:,:,2);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,3); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('higgs-ridge PCoh')
%%
X  = sol_higgs{3,1}(:,:,3);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,4); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('higgs-naive PCoh')
%%
X  = sol_higgs{3,1}(:,:,4);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,5); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('eloreta-hglasso PCoh')
%%
X  = sol_higgs{3,1}(:,:,5);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,3,6); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('lcmv-hglasso PCoh')
%%
colormap(cmap);

saveas( figure_partial_coherences,strcat(output_source,filesep,'partial_coherence_maps_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  Partial Coherence Maps to  ---> ', output_source) );
delete(figure_partial_coherences);
end