function [y, t_ds, newFs]=calc_envelope(y,Fs)
if nargin<2
    Fs=200;
end
[m,n,s]=size(y);%m:variabel,n:length,s:segments or trials

T=double(n)/Fs;
time=0:1/Fs:T-(1/Fs);
timeRange=[];iSession=[];
[t, tI, tR] = time_range(time, timeRange, iSession);



Settings.EnvelopeParams.takeLogs    = false;                           % perform analysis on logarithm of envelope. This improves normality assumption
Settings.EnvelopeParams.absolute    = false;
Settings.EnvelopeParams.downsample  = false;
Settings.EnvelopeParams.windowLength = 1/20; % 1/20 s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012.
Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average
Settings.EnvelopeParams.saveMemory  = false;

% [y_temp, ~, ~]= ROInets.envelope_data(squeeze(y(:,:,1)),t,Settings.EnvelopeParams);
% y_env=zeros(m,size(y_temp,2),s);

% for i=1:s
%     [y_env(:,:,i), t_ds, newFs]= ROInets.envelope_data(squeeze(y(:,:,i)),t,Settings.EnvelopeParams);
% end
y=permute(y,[1,3,2]);
y=reshape(y,[m*s,n]);
% y= osl_filter(y,band,'fs',Fs);
if test_gpu([],m*s*n*4)
    [y, t_ds, newFs]= ROInets.envelope_data(gpuArray(y),t,Settings.EnvelopeParams);
else
    [y, t_ds, newFs]= ROInets.envelope_data(y,t,Settings.EnvelopeParams);
end
n1=size(y,2);
y=reshape(y,[m,s,n1]);
y=permute(y,[1,3,2]);





end
