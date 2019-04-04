%% GGM Local Quadratic Approximation
function [PC,llh] = hggm_lasso_ssbl(EC,m,A,rho,maxiter)
%% sggm-lqa initial values 
q                 = length(EC);
Iq                = eye(q);
m2                = m^2;
rho2              = rho^2;
nu                = 2.5E2;
A2                = A.^2;
EC                = EC/sqrt(mean(abs(diag(EC*EC'))));
ECinv             = Iq/EC; 
ECinv             = (ECinv + ECinv')/2; 
idx               = (A > 0);
idx0              = (A == 0);
gamma2            = zeros(q);
llh               = zeros(maxiter,1);
%%
gamma2(idx0)      = rho*m*abs(ECinv(idx0)).^2;
PC                = ECinv;
%% Main cycle
for k_inner = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior
    DET              = 1 + 4*m2*rho2*A2(idx).*abs(PC(idx)).^2;
    gamma2(idx)      = (sqrt(DET) - 1)./(2*m*rho2*A2(idx));
    %%
    %% Standarization of the Empirical Covariance Matrix
    ninf             = max(gamma2(:));
    gamma2ninf       = gamma2/ninf;
    st_factor1       = nu*ninf^(-1/2)*gamma2ninf.^(1/2).*ECinv; % lqa corrected factor 1
    st_factor2       = (1 + gamma2ninf*nu); % lqa corrected factor 2
    %%
    ECst_inv         = st_factor1./st_factor2; % Inverted Standard Empirical Covariance Matrix
    ECst             = Iq/ECst_inv; % Standard Empirical Covariance Matrix
    ECst             = (ECst + ECst')/2; 
    %%
    %% Standard Precision Matrix estimator
    [U,D,V]          = svd(ECst);
    d                = abs(diag(D));
    PCst             = (1/2)*U*diag(sqrt(d.^2 + 4) - d)*V';
    PCst             = (PCst + PCst')/2;
    %%
    %% Unstandarized Precision Matrix estimator
    PC_tmp           = gamma2.^(1/2).*PCst;
    %%
    [U,D,V]          = svd(PC_tmp);
    d                = abs(diag(D));
    llh_tmp          = sum(log(d)) - sum(abs(diag(PC_tmp*EC))) -rho*sum(abs(A(:).*PC_tmp(:)));
    if k_inner == 1
        PC        = PC_tmp;
        llh(k_inner) = llh_tmp;
    elseif k_inner > 1 && llh_tmp >= llh(k_inner-1)
        PC        = PC_tmp;
        llh(k_inner) = llh_tmp;
    else
        llh(k_inner:end) = llh(k_inner-1);
        break
    end
end
end