function [theta0,sigma,theta,process,process_energe,SigmaEpir,ThetaEpir] = gen_processes(Nsegments,Nseed,options)
%%
% Authors:
% - Ying Wang
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: Nov 24, 2019



%% generates time series from a given spectral tensor
[~,theta0]     = gen_hggm02(Nseed,options);
[sigma,theta]  = spectral_precicion_tensor(theta0,options.nu0(2),...
options.tstudent_a(2),options.tstudent_b(2),options.tstudent_dof(2),...
options.deltaf,options.Fmin,options.Fmax);
process_energy = sigma_to_band_energy((options.Ntpoints)*sigma,...
    [options.band(2,1),options.band(2,2)],options.deltaf);
[fourier_coef] = gen_tensor_hggm((options.Ntpoints)*sigma,Nsegments);% (options.Ntpoints)*
process        = iFFT_time_series(fourier_coef);% sqrt(options.Ntpoints)*

%%
% process_energe2=fourier_coef_to_band_energy(fourier_coef,...
%     [options.band(2,1),options.band(2,2)],options.deltaf);
% plot(process_energe2(1,:)')
% hold on;plot(process_energe1(1,:)')

%% 1 single frequency
% p=squeeze(sqrt((options.Ntpoints)*sigma(1,1,:)));
% % p=squeeze(sqrt(sigma(1,1,:)));
% 
% X=sum(abs(fourier_coef(1,:,:)),3)/(Nsegments);
% plot(options.F,[p,X'])

%%
% SE=tsdata_to_cpsd(gather(process),options.Nfreqs-1);
% figure;plot_cpsd(sigma,{'model'},options.Fs,[],false);
% figure;plot_cpsd(SE,{'data'},options.Fs,[],false);
% figure;plot_cpsd(cat(4,sigma,SE),{'model','data'},options.Fs,[],false);
end
