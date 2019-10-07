function result = Main_Jankova_test(output_sourse)
%%

%%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: April 4, 2019


% Updates
% - Ariosky Areces Gonzalez

% Date: April 4, 2019

%%
% Simulates a q-size Partial Correlation (PC) matrix made up a predefined
% number of blocks (nblocks) and uses m-samples of an Cyrcularly Symmetric
% Complex Gaussian (CSCG) generator to compute the Empirical Covariance
% (EC) matrix.
%%

%% Initialize generator parameters

output_sourse = strcat(output_sourse, filesep,'Jankova_test');
if(~isfolder(output_sourse))
mkdir(output_sourse);
end

m              = 600;           % Sample number
q              = 60;              % Number of variables
nblocks        = 3;               % Number of blocks in simulation
options.config = 2;               % (2) overlapping blocks (1) nonoverlapping blocks
options.var    = 2;               % (2) complex variable (1) real variable
ntrial         = 100;              % Number of trials
maxiter        = 30;              % Maximum number of iterations
rho            = sqrt(log(q)/m);  % Global Regularization parameter
rho_diag       = 0;               % Selection Mask diagonal
rho_ndiag      = 1;               % Selection Mask nondiagonal
lambda         = rho_diag*eye(q)+rho_ndiag*(ones(q)-eye(q));  % Selection Mask squared
z_stat         = [];              % z-statistic
llh            = zeros(ntrial,maxiter);
%%
process_waitbar = waitbar(0,'Please wait...');
for trial = 1:ntrial
    waitbar(trial/ntrial,process_waitbar,strcat('simulation # ',num2str(trial)));
    disp(['simulation # ',num2str(trial)]);
    %% Circularly Symmetric Complex Gaussian generator
    [EC,Data,PC_sim] = gen_hggm1(m,q,nblocks,options);
    %%
    scale_EC         = sqrt(mean(abs(diag(EC*EC'))));
    EC               = EC/scale_EC;
    %%
    % Estimates the PC by Gaussian Graphical Model Local Quadratic
    % Approximation (GGM_LQA).
    [PC,llh(trial,:)] = hggm_lasso_ssbl(EC,m,lambda,rho,maxiter);
    % Unbiased statistics
    PC_unb           = 2*PC - PC*(EC)*PC; % Unbiased PC
    PC_var           = sqrt(abs(diag(PC))*abs(diag(PC))' + abs(PC).^2); % Variances of PC_unb
    %% Plot trial
    if trial == 1
       hggm_lasso_PC_figure = figure;
        X            = PC_sim;
        X            = X - diag(diag(X));
        X            = X/max(abs(X(:)));
        subplot(2,2,1); imagesc(abs(X));
        ylabel('generators')
        xlabel('generators')
        title('simulated PC')
        X            = PC;
        X            = X - diag(diag(X));
        X            = X/max(abs(X(:)));
        subplot(2,2,2); imagesc(abs(X));
        ylabel('generators')
        xlabel('generators')
        title('hggm-lasso PC')
        X            = PC_unb;
        X_var        = PC_var;
        rth          = 3.16;
        rayleigh_mask = find(abs(X) < (rth/sqrt(m))*(X_var - diag(diag(X))));
        X(rayleigh_mask) = 0;
        X            = X - diag(diag(X));
        X            = X/max(abs(X(:)));
        subplot(2,2,3); imagesc(abs(X));
        ylabel('generators')
        xlabel('generators')
        title('Rayleigh corrected PC')
        X            = PC_unb;
        X            = X - diag(diag(X));
        X            = X/max(abs(X(:)));
        subplot(2,2,4); imagesc(abs(X));
        ylabel('generators')
        xlabel('generators')
        title('unbiased PC')
        colormap('hot')
        
        saveas( hggm_lasso_PC_figure,strcat(output_sourse,filesep,'simulated PC.fig'));
        disp(strcat('Saving figure ---->  simulated PC.fig  to  ---> ', output_sourse) );
        delete(hggm_lasso_PC_figure);
    end
    % Removing indices where Null-hypotheses does not hold
    ind              = find(abs(PC_sim) > 0);
    PC(ind)          = 0;
    PC_unb(ind)      = 0;
    PC_var(ind)      = 0;
    % Select variances over certain theshold to avoid numerical unstability
    ind_var          = find(PC_var/max(PC_var(:)) > 0.01);
    % PC_unb rated by its Variances (Z-transformed)
    z_stat_new       = sqrt(m)*abs(PC_unb(ind_var))./PC_var(ind_var);
    z_stat           = [z_stat;z_stat_new];
end
delete(process_waitbar);


 

h_hggm_lasso_likelihood_figure = figure;
plot(llh')
ylabel('likelihood')
xlabel('iterations')
title('h-hggm-lasso likelihood')

 saveas( h_hggm_lasso_likelihood_figure,strcat(output_sourse,filesep,'h-hggm-lasso likelihood.fig'));
 disp(strcat('Saving figure ---->  h-hggm-lasso likelihood.fig  to  ---> ', output_sourse) );
 delete(h_hggm_lasso_likelihood_figure);
%%
pdf_cdf_figure = figure;
subplot(1,2,1)
h        = histogram(z_stat(:),'Normalization','pdf');
bins     = h.BinEdges;
pdf      = 2*bins.*exp(-bins.^2);
hold on
plot(bins,pdf,'LineWidth',2);
title('Probability Density Function (pdf)')
xlabel('z-stat'); ylabel('pdf(z-stat)')
legend('z-stat histogram','Rayleigh distribution')


%%
subplot(1,2,2)
h        = histogram(z_stat(:),'Normalization','cdf');
bins     = h.BinEdges;
cpf      = 1- exp(-bins.^2);
hold on
plot(bins,cpf','LineWidth',2);
title('Commulative Distribution Function (cdf)')
xlabel('z-stat'); ylabel('cdf(z-stat)')
legend('z-stat histogram','Rayleigh distribution','location','northwest')


 saveas( pdf_cdf_figure,strcat(output_sourse,filesep,'PDF & CDF.fig'));
 disp(strcat('Saving figure ---->  Probability Density Function and Commulative Distribution Function to  ---> ', output_sourse) );
delete(pdf_cdf_figure);

%msgbox("The simulation HGGM Convergence & Jankova conditions has been suceffuly.. ","Result")
end

