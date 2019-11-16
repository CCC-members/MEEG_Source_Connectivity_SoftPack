function [DataXi,noise_eneger]=gen_process_xi(Nnoise,Ntpoints,Nsegments,Nsubj,NoiseIdx,options)

%% Xi process
% Axi= gen_Nnvar(options.vertices);
Axi= gen_NnvarAdjacent(options.vertices,options.faces);
% Axi=
Axi=Axi(NoiseIdx,NoiseIdx);
if test_gpu([],Nnoise*Ntpoints*Nsegments*Nsubj*4)
    DataXi=zeros(Nnoise,Ntpoints,Nsegments,Nsubj,'single','gpuArray');
else
    DataXi=zeros(Nnoise,Ntpoints,Nsegments,Nsubj);
    
end
for i=1:Nsubj
    DataXi(:,:,:,i)=gen_var_time_series(full(Axi),eye(size(Axi,1)),Ntpoints,Nsegments,options.Fs);
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