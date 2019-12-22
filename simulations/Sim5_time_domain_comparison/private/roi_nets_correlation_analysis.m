function CorrMats=roi_nets_correlation_analysis(Data,envData,rho_ref)


    Settings.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
%     Settings.Regularize.path          = logspace(log10(rho_ref/10),log10(rho_ref/10),log10(rho_ref));              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
    Settings.Regularize.path          = logspace(-9,2,80);              % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization.
    Settings.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
    Settings.Regularize.adaptivePath  = true;                           % adapth the regularization path if necessary
%     Settings.Regularize.method        = 'Bayesian';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
%     Settings.Regularize.Prior         = struct('a', 1./3, 'b', 0);
    CorrMats = ROInets.run_correlation_analysis(Data,     ...
                                                envData,      ...
                                                Settings.Regularize);
end






