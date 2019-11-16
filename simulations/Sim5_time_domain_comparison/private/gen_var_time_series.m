function [y,res]=gen_var_time_series(A,Sigma,numTimepoints,numTrials,sampleRate,numWarmTimepoints)

tic
% assert(isVarStable(A),'VAR Coef not stable');
if nargin<5
    sampleRate=100;
end
if nargin<6
    numWarmTimepoints=100;
end

[numNodes,~,maxLag]=size(A);

% disp(['<gen_var_time_seriese.m> Generating tir',num2str(numTrials),...
%     'trials signal with ', num2str(numNodes), ' nodes ', num2str(maxLag), ' slices cross spectrum'])
if isempty(Sigma)
    Sigma=eye(numNodes);
end

timeMax=numTimepoints/sampleRate;


% warmTime=1;
% numWarmTimepoints=warmTime*sampleRate;

warmTime=numWarmTimepoints/sampleRate;
wholeTime=timeMax+warmTime;
wholeNumTimepoints=numTimepoints+numWarmTimepoints;

t=0:1/sampleRate:(wholeTime-1/sampleRate);
c=zeros(numNodes,length(t));
% if isa(A,'gpuArray')
if test_gpu([],numNodes*wholeNumTimepoints*numTrials*4*2) 
    fprintf('<gen_var_time_series.m> Generating signals with variable(%d) order(%d) vector autoregressive model,trials=%d, numTimepoints=%d, using GPU;\n',...
        numNodes,maxLag,numTrials,numTimepoints);
    A=gpuArray(A);
    Sigma=gpuArray(Sigma);
    y=zeros(numNodes,wholeNumTimepoints,numTrials,'single','gpuArray');
    if isdiag(Sigma)
        res=randn(numNodes,wholeNumTimepoints,numTrials,'single','gpuArray')+c;
    else
        res=zeros(numNodes,wholeNumTimepoints,numTrials,'single','gpuArray');
        for i=1:numTrials
            res(:,:,i)=mvnrnd(zeros(numNodes,1),Sigma,wholeNumTimepoints,'single','gpuArray')'+c;
        end
    end
else
    fprintf('<gen_var_time_seriese.m> Generating signals with variable(%d) order(%d) vector autoregressive model,trials=%d, numTimepoints=%d, using CPU;\n',...
        numNodes,maxLag,numTrials,numTimepoints);
    y=zeros(numNodes,wholeNumTimepoints,numTrials);
    if isdiag(Sigma)
        res=randn(numNodes,wholeNumTimepoints,numTrials)+c;
    else
        res=zeros(numNodes,wholeNumTimepoints,numTrials);
        for i=1:numTrials
            res(:,:,i)=mvnrnd(zeros(numNodes,1),Sigma,wholeNumTimepoints)'+c;
        end
    end
end
y=permute(y,[1,3,2]);
res=permute(res,[1,3,2]);
y=var_time_series(A,y,res);
y=permute(y,[1,3,2]);
res=permute(res,[1,3,2]);
y=y(:,numWarmTimepoints+1:end,:);
res=res(:,numWarmTimepoints+1:end,:);
toc
end