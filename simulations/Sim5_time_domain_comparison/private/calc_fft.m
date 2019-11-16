function FtCoef=calc_fft(Data,singleSide,average)
if nargin<2 || isempty(singleSide)
    singleSide=true;
end
[numNodes,numnumTimepoints,numTrials]=size(Data);
% Data1=ifft(FtCoef,[],2);
Data=permute(Data,[1,3,2]);%numNodes,numTrials,numFreqs
% FtCoef=(numFreqs/6)*reshape(FtCoef,[numNodes*numTrials,numFreqs]);
% FtCoef=(reshape(FtCoef,[numNodes*numTrials,numFreqs]).^2)/2;
Data=reshape(Data,[numNodes*numTrials,numnumTimepoints]);


if numNodes*numTrials*numnumTimepoints>1e8
    Data=single(Data);
end

if test_gpu
    tic
    Data=gpuArray(Data);
    FtCoef=fft(Data,[],2);
    FtCoef=gather(real(FtCoef));
    toc
% if numNodes<1000
else
    tic
    FtCoef=fft(Data,[],2);
    toc
end
if singleSide
    numFreqs=numnumTimepoints/2+1;
    FtCoef=FtCoef(:,1:numFreqs);
else
    numFreqs=numnumTimepoints;
end
FtCoef=reshape(FtCoef,[numNodes,numTrials,numFreqs]);
FtCoef=permute(FtCoef,[1,3,2]);%numNodes,numFreqs,numTrials
if nargin>2&& average
    FtCoef=sum(abs(FtCoef),3);
end

end