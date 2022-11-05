function result = Main_Jankova_test(output_source)
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

output_source = strcat(output_source, filesep,'Jankova_test');
if(~isfolder(output_source))
    mkdir(output_source);
end

sizes  = [30,100,1000];
trial  = 1;
ntrial = length(sizes);
hgLASSO_figure = figure;
                             process_waitbar = waitbar(0,'Please wait...');
for q = sizes
                             waitbar(trial/ntrial,process_waitbar,strcat('size # ',num2str(q)));
    param.q              = q;
    param.m              = 100*q;                % Sample number
    param.nblocks        = 3*floor(log10(q));    % Number of blocks in simulation
    param.config         = 2;                    % (2) overlapping blocks (1) nonoverlapping blocks
    param.var            = 2;                    % (2) complex variable (1) real variable
    param.maxiter        = 30;                   % Maximum number of iterations
    param.a              = sqrt(log(q)/param.m); % Global Regularization parameter
    param.nu             = param.m;
    a_diag               = 0;                    % Selection Mask diagonal
    a_ndiag              = 1;                    % Selection Mask nondiagonal
    A                    = a_diag*eye(q)+a_ndiag*(ones(q)-eye(q));  % Selection Mask squared
    param.A              = A;
    param.run_bash_mode  = 1;
    param.use_gpu        = 1;
    param.eigreg         = 1E-4;
    %% Circularly Symmetric Complex Gaussian generator
    [Psi,Data,Theta_sim] = gen_hggm1(param);
    %%
    % Estimates the PC by Gaussian Graphical Model Local Quadratic
    % Approximation (GGM_LQA).
    Psi.X = Psi.X/sqrt(sum(Psi.d.^2)/q);
    Psi.d = Psi.d/sqrt(sum(Psi.d.^2)/q);
    [Theta,llh] = hg_lasso_lqa(Psi,param);
    % Unbiased statistics
    Theta_unb        = 2*Theta.X - Theta.X*(Psi.X)*Theta.X; % Unbiased PC
    Theta_var        = sqrt(abs(diag(Theta.X))*abs(diag(Theta.X))' + abs(Theta.X).^2); % Variances of PC_u
    %% Plot
    figure(hgLASSO_figure)
    set(gcf,'Position',[50 50 1400 700]);
    %
    X            = Theta_sim.X;
    X            = X - diag(diag(X));
    X            = X/max(abs(X(:)));
    subplot(ntrial,5,(trial-1)*5+1); imagesc(abs(X));
    ylabel('nodes')
    xlabel('nodes')
    title('Simulation')
    %
    X            = Theta.X;
    X            = X - diag(diag(X));
    X            = X/max(abs(X(:)));
    subplot(ntrial,5,(trial-1)*5+2); imagesc(abs(X));
    ylabel('nodes')
    xlabel('nodes')
    title('hgLASSO')
    %
    X            = Theta_unb;
    X            = X - diag(diag(X));
    X            = X/max(abs(X(:)));
    subplot(ntrial,5,(trial-1)*5+3); imagesc(abs(X));
    ylabel('nodes')
    xlabel('nodes')
    title('Debiasing')
    %
    X                = Theta_unb;
    X_var            = Theta_var;
    % Removing indices where Null-hypotheses does not hold
    ind              = find(abs(Theta_sim.X) > 0);
    Theta_unb(ind)   = 0;
    Theta_var(ind)   = 0;
    % Select variances over certain theshold to avoid numerical unstability
    ind_var          = find(Theta_var/max(Theta_var(:)) > 1E-2);
    % PC_unb rated by its Variances (Z-transformed)
    z_stat           = sqrt(param.m)*abs(Theta_unb(ind_var))./Theta_var(ind_var);
    subplot(ntrial,5,(trial-1)*5+4);
    h        = histogram(z_stat(:),'Normalization','pdf');
    bins     = h.BinEdges;
    pdf      = 2*bins.*exp(-bins.^2);
    hold on
    plot(bins,pdf,'LineWidth',2);
    xlim([0 5]);
    title('Statistics')
    xlabel('amplitudes'); ylabel('probabilities')
    legend('Empirical','Rayleigh')
    %
    rth              = 3.16;
    rayleigh_mask = find(abs(X) < (rth/sqrt(param.m))*(X_var - diag(diag(X))));
    X(rayleigh_mask) = 0;
    X                = X - diag(diag(X));
    X                = X/max(abs(X(:)));
    subplot(ntrial,5,(trial-1)*5+5); imagesc(abs(X));
    ylabel('nodes')
    xlabel('nodes')
    title('Correction')
    %
    load('colormap3')
    colormap(cmap)
    trial = trial + 1;
end

delete(process_waitbar);
saveas( hgLASSO_figure,strcat(output_source,filesep,'.fig'));
disp(strcat('Saving figure ---->  hgLASSO_figure.fig  to  ---> ', output_source) );
delete(hgLASSO_figure);

% msgbox("The simulation HGGM Convergence & Jankova conditions has been suceffuly.. ","Result")

end