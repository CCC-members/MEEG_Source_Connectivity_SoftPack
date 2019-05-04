function [] = Results(output_sourse)

%% Run Results
%% Loading Simulation Substrate or Real Data
load('colormap2.mat')
load('Pseudorand_Net.mat')
measures_label    = {'auc' 'sens' 'spec' 'prec' 'f1'};
methods_label     = {'h-hggm-lasso';'h-hggm-ridge';'h-hggm-naive';'eloreta-hggm';'lcmv-h-hggm'};
Nmeasures         = length(measures_label);
%%
%% h_hggm
load('Solutions_h_hggm.mat')
[measures]        = quality_measures(sol_h_hggm,Theta_sim);
measures(isnan(measures)) = 0.5;
measures_mean     = round(100*mean(measures,3));
measures_std      = round(100*std(measures,0,3));
%%
%%
process_waitbar = waitbar(0,'Please wait...');
for meth = 1:5
    for meas = 1:5
        waitbar((meth*meas)/(25),process_waitbar,strcat('Outputing results....'));
        measures_h_hggm{meth,meas}  = [num2str(measures_mean(meth,meas)) '+/-' num2str(measures_std(meth,meas))];
    end
end
delete(process_waitbar);
%%
%% Table with quality measures
Table = cell(6,6);
Table(2:end,1)    = methods_label;
Table(1,2:end)    = measures_label;
Table(2:end,2:end)  = measures_h_hggm;


save(strcat(output_sourse,filesep,'Table_sens_system_',sens_system,'.mat'));
disp(strcat('Saving Table ---->  Table with quality measures to  ---> ', output_sourse) );

    
%%
%% Plot likelihood


Nsim = ceil(size(sol_h_hggm,2));
figure_likelihood = figure; 

subplot(3,2,1);
llh = zeros(length(sol_h_hggm{4,1}{1}{1}),Nsim);

for sim = 1:Nsim
  
    llh(:,sim) = sol_h_hggm{4,sim}{1}{1};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('h-hggm-lasso likelihood')


%%
%%

subplot(3,2,2); 
llh = zeros(length(sol_h_hggm{4,1}{1}{2}),Nsim);
for sim = 1:Nsim
   
    llh(:,sim) = sol_h_hggm{4,sim}{1}{2};  
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('h-hggm-ridge likelihood')


%%

subplot(3,2,3); 
llh = zeros(length(sol_h_hggm{4,1}{1}{3}),Nsim);
for sim = 1:Nsim
    
    llh(:,sim) = sol_h_hggm{4,sim}{1}{3};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('h-hggm-naive likelihood')


%%


subplot(3,2,4); 
for sim = 1:Nsim
    [gcv_opt,idx_gamma]       = min(sol_h_hggm{4,sim}{3});
    plot(sol_h_hggm{4,1}{2},sol_h_hggm{4,sim}{3},...
        '-',sol_h_hggm{4,1}{2}(idx_gamma),...
        gcv_opt,'b*');
    hold on;
end

ylabel('gcv value')
xlabel('regularization parameter')
title('eloreta-hggm gcv function')

%% Plot corticaL map 
cortex.vertices = vertices;
cortex.faces    = faces;
[qL,qR,qfull,indvL,indvR,indv,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos_trans] = split_hemispheres(cortex,elec_pos);
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
end

saveas( figure_likelihood,strcat(output_sourse,filesep,'h-hggm likelihood_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  h-hggm likelihood to  ---> ', output_sourse) );
delete(figure_likelihood);


%%
figure_partial_correlations = figure;
%% Plot partial correlations
X  = Theta_sim{1};
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,2,1); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('simulated partial correlations')
%%
X  = sol_h_hggm{3,1}(:,:,1);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,2,2); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('h-hggm-lasso partial correlations')
%%
X  = sol_h_hggm{3,1}(:,:,2);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,2,3); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('h-hggm-ridge partial correlations')
%%
X  = sol_h_hggm{3,1}(:,:,3);
X  = X - diag(diag(X));
X  = X/max(abs(X(:)));
subplot(2,2,4); imagesc(abs(X)); 
ylabel('generators')
xlabel('generators')
title('h-hggm-naive partial correlations')
%%
colormap('hot');

saveas( figure_partial_correlations,strcat(output_sourse,filesep,'partial_correlation_maps_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  Partial Correlation Maps to  ---> ', output_sourse) );
delete(figure_partial_correlations);
end