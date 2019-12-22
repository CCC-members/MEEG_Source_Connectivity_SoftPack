function [DataXi,noise_eneger]=gen_process_and_theta_xi(Nnoise,Ntpoints,Nsegments,Nsubj,NoiseIdx,options)

%% Xi process
Axi= gen_NnvarAdjacent(options.vertices,options.faces);
Axi=Axi(NoiseIdx,NoiseIdx);
%% VAR process
NwarmTimepoints=100;
% if test_gpu([],Nnoise*(Ntpoints+NwarmTimepoints)*Nsegments*Nsubj*4*2)
%     DataXi=zeros(Nnoise,Ntpoints,Nsegments,Nsubj,'single','gpuArray');
% else
%     DataXi=zeros(Nnoise,Ntpoints,Nsegments,Nsubj);
% end
for i=1:Nsubj
    DataXi(:,:,:,i)=gen_var_time_series(full(Axi),eye(size(Axi,1)),Ntpoints,Nsegments,options.Fs,NwarmTimepoints);
end
%%
band_range=(pi/(options.Nfreqs-1))*[options.band(2,1)/options.deltaf:1:options.band(2,2)/options.deltaf];
H=calc_var_to_transfe(Axi,band_range);
for i=1:size(H,3)
    Sxi(i)=norm(H(:,:,i),'fro').^2/size(Axi,1);% because the coefficient is a toeplitz matrix, so all nodes have same amplitude
end
noise_eneger=sqrt(Sxi*options.Ntpoints);
noise_eneger= repmat(noise_eneger,[Nnoise,1]);
%%
% band_range=(pi/(options.Nfreqs-1))*options.F;
% [H,K]=calc_var_to_transfe(Axi,band_range);
% U = eye(size(K,1)) - K;
% for count = 1:length(options.F)
%     theta_xi(:,:,count) = U(:,:,count)'*U(:,:,count);
%     sigma_xi(:,:,count) = inv(theta_xi(:,:,count));
% end
% save('D:\code\MultiLagEegSimulation\LeakageCorrectionSimulation\result\DTF_xi.mat','-v7.3')
end
