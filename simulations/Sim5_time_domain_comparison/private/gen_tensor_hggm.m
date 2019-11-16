function [FtCoef,SigmaEpir,ThetaEpir]=gen_tensor_hggm(Sigma,Ntrial)
[n,~,t]=size(Sigma);
if nargout>2
    SigmaEpir= ndSparse.build([n,n,t+1]  );
    CEpir= ndSparse.build([n,n,t+1] );
end
FtCoef=ones(Ntrial,n,t)+1i*ones(Ntrial,n,t);
%%
tic
disp(['<gen_HggmToFtcoef.m> Generating ',num2str(Ntrial),' fourier coefficient samples from ', num2str(n), ' nodes ', num2str(size(Sigma,3)), ' slices cross spectrum;'])
%%
[userview systemview]=memory;
%%
if test_gpu
    FtCoef=gpuArray(FtCoef);
    for i=1:t
        Wisomph     = (1/2)*[real(Sigma(:,:,i)) -imag(Sigma(:,:,i)); imag(Sigma(:,:,i)) real(Sigma(:,:,i))];
        FtCoef_isomph = mvnrnd(zeros(1,2*n,'gpuArray'),gpuArray(full(Wisomph)),Ntrial);
        FtCoefRe      = FtCoef_isomph(:,1:n);
        FtCoefIm      = FtCoef_isomph(:,n+1:2*n);
        FtCoef(:,:,i) = FtCoefRe + 1i*FtCoefIm;% samples*nodes*frequencies
        %         SigmaEpir(:,:,i)           = (1/Ntrial)*sparse(FtCoef(:,:,i).'*conj(FtCoef(:,:,i)));%% sample covariance
        %         CEpir(:,:,i)           = (1/Ntrial)*(FtCoef(:,:,i).'*(FtCoef(:,:,i)));%% sample pseudo covariance
    end
    FtCoef=gather(FtCoef);
elseif t*n<100000 || userview.MemAvailableAllArrays<1E10
    %     if 1
    for i=1:t
        Wisomph     = (1/2)*[real(Sigma(:,:,i)) -imag(Sigma(:,:,i)); imag(Sigma(:,:,i)) real(Sigma(:,:,i))];
        FtCoef_isomph = mvnrnd(zeros(1,2*n),sparse(Wisomph),Ntrial);
        FtCoefRe      = FtCoef_isomph(:,1:n);
        FtCoefIm      = FtCoef_isomph(:,n+1:2*n);
        FtCoef(:,:,i)        = FtCoefRe + 1i*FtCoefIm;% samples*nodes*frequencies
        %         SigmaEpir(:,:,i)           = (1/Ntrial)*sparse(FtCoef(:,:,i).'*conj(FtCoef(:,:,i)));%% sample covariance
        %         CEpir(:,:,i)           = (1/Ntrial)*(FtCoef(:,:,i).'*(FtCoef(:,:,i)));%% sample pseudo covariance
    end
else
    startmatlabpool
    parfor i=1:t
        Wisomph     = (1/2)*[real(Sigma(:,:,i)) -imag(Sigma(:,:,i)); imag(Sigma(:,:,i)) real(Sigma(:,:,i))];
        FtCoef_isomph = mvnrnd(zeros(1,2*n),sparse(Wisomph),Ntrial);
        FtCoefRe      = FtCoef_isomph(:,1:n);
        FtCoefIm      = FtCoef_isomph(:,n+1:2*n);
        FtCoef(:,:,i) = FtCoefRe + 1i*FtCoefIm;% samples*nodes*frequencies
        %         SigmaEpir(:,:,i)           = (1/Ntrial)*sparse(FtCoef(:,:,i).'*conj(FtCoef(:,:,i)));%% sample covariance
        %         CEpir(:,:,i)           = (1/Ntrial)*(FtCoef(:,:,i).'*(FtCoef(:,:,i)));%% sample pseudo covariance
    end
    closematlabpool
end
toc
FtCoef=permute(FtCoef,[2,3,1]);% nodes*frequencies*samples
end
