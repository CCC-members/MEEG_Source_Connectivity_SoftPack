function process = iFFT_time_series(fourier_coef,single_side)
% Authors:
% - Ying Wang

% Date: Nov 24, 2019
if nargin<2
    single_side=true;
end
[numNodes,numFreqs,numTrials]=size(fourier_coef);
% Data1=ifft(FtCoef,[],2);
fourier_coef=permute(fourier_coef,[1,3,2]);%numNodes,numTrials,numFreqs
% FtCoef=(numFreqs/6)*reshape(FtCoef,[numNodes*numTrials,numFreqs]);
% FtCoef=(reshape(FtCoef,[numNodes*numTrials,numFreqs]).^2)/2;
fourier_coef=reshape(fourier_coef,[numNodes*numTrials,numFreqs]);
disp(['<iFFT_time_series.m> Doing inverse fourier transform for ',num2str(numNodes*numTrials),...
    ' trials(',num2str(numNodes),'*',num2str(numTrials),') ',...
    num2str(numFreqs),' frequency points;'])
if single_side
    if numNodes*numTrials*numFreqs>1e8
        fourier_coef=single(fourier_coef);
    end
    fourier_coef=cat(2,fourier_coef(:,1:end-1),zeros(numNodes*numTrials,1),conj(fourier_coef(:,end-1:-1:2)));
end

if test_gpu
    tic
    fourier_coef=gpuArray(fourier_coef);
    process = real(ifft(fourier_coef,[],2));
%     Data=gather(real(Data));
    toc
% if numNodes<1000
else
    tic
    process=ifft(fourier_coef,[],2);
    toc
end
if single_side
    numnumTimepoints=2*(numFreqs-1);
else
    numnumTimepoints=numFreqs;
end
process=reshape(process,[numNodes,numTrials,numnumTimepoints]);
process=permute(process,[1,3,2]);%numNodes,numnumTimepoints,numTrials


end
