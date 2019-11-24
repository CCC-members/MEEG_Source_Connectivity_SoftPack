function [DataXi,noise_eneger]=gen_process_xi(Nnoise,Ntpoints,Nsegments,Nsubj,NoiseIdx,options)
%%
% Authors:
% - Ying Wang
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: Nov 24, 2019


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
    Sxi(i)=norm(H(:,:,i),'fro').^2/size(Axi,1);% all nodes have same noise amplitude
end
noise_eneger=sqrt(Sxi*options.Ntpoints);
noise_eneger= repmat(noise_eneger,[Nnoise,1]);
end