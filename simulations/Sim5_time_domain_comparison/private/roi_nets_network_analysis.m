function Mat=roi_nets_network_analysis(V_filt,Tjv,Fs,rho_ref)
% Authors:
% - Ying Wang

% Date: Nov 24, 2019
%% get source time series using Inverse projector and filtered data
sourceSignal=zeros(size(Tjv,1),size(V_filt,2),size(V_filt,3));
for i= 1:size(V_filt,3)
    sourceSignal(:,:,i)=Tjv*V_filt(:,:,i);
end
%% apply orthogonalization to source singnal 
sourceSignal = calc_orthogonalize(sourceSignal,'householder');

%% apply hilbert envelope to narrow band source signal
[sourceSignalEnv, ~, ~]=calc_envelope(sourceSignal,Fs);
sourceSignalEnv=abs(sourceSignalEnv);

%% merge multi-segments data to one trial
sourceSignal=reshape(sourceSignal,[size(sourceSignal,1),size(sourceSignal,2)*size(sourceSignal,3)]);
sourceSignalEnv=reshape(sourceSignalEnv,[size(sourceSignalEnv,1),size(sourceSignalEnv,2)*size(sourceSignalEnv,3)]);

%% apply graphic lasso to estimate source covariance
CorrMatsCorrMats=roi_nets_correlation_analysis(sourceSignal,sourceSignalEnv,rho_ref);

%% if return a multi-trial precision 
Mat=mean(CorrMatsCorrMats.envPrecisionRegularized,3);


end
