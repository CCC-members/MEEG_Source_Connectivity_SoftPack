function [y, t_ds, newFs]=calc_envelope(y,Fs,F_deman)
if nargin<3
    F_deman=true;
end
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
Settings.EnvelopeParams.saveMemory  = false;
% set downsample to false then this two doesn't metter
    Settings.EnvelopeParams.windowLength = 1/20; % 1/20 s                       % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012.
    Settings.EnvelopeParams.useFilter    = true;                        % use a more sophisticated filter than a sliding window average


y=permute(y,[1,3,2]);
y=reshape(y,[m*s,n]);
if test_gpu([],m*s*n*4)
    [y, t_ds, newFs]= ROInets.envelope_data(gpuArray(y),t,Settings.EnvelopeParams);
    if F_deman
        y=ROInets.demean(y,2);
    end
else
    [y, t_ds, newFs]= ROInets.demean(ROInets.envelope_data(y,t,Settings.EnvelopeParams),2);
    if F_deman
         y=ROInets.demean(y,2);
    end
end
n1=size(y,2);
y=reshape(y,[m,s,n1]);
y=permute(y,[1,3,2]);





end
