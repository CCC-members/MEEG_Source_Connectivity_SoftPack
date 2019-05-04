function [LeadFields,Nv,Nsubj,Nsamp,Nnoise,db_source,db_sens,Nsim,Nseed,options,index_seed,index_full,vertices,faces,elec_pos] = Initialize_simulation(sens_system)
%% Generates multiple subject Data from common spatio-temporal-trial simulated sources
%% Pick Lead Fields for analysis
if strcmp(sens_system,'large') == 1
    load ('LeadFields_large.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
    subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
    LeadFields = LeadFields(1,subject);
elseif strcmp(sens_system,'small') == 1
    load ('LeadFields_small.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
    subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
    LeadFields = LeadFields(1,subject);
elseif strcmp(sens_system,'pseudo') == 1
    load ('LeadFields_pseudo.mat'); % Load Lead Field cell array (1 X number of Lead Fields)
    subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
    LeadFields = LeadFields(1,subject);
end
%%
%% Surface
if strcmp(sens_system,'large') == 1
    load('HeadModel_large.mat'); % Load mesh vertices (all Lead Fields should be given on homemorph surfaces)
elseif strcmp(sens_system,'small') == 1
    load('HeadModel_small.mat'); % Load mesh vertices (all Lead Fields should be given on homemorph surfaces)
elseif strcmp(sens_system,'pseudo') == 1
    load('HeadModel_pseudo.mat'); % Load mesh vertices (all Lead Fields should be given on homemorph surfaces)
end
%%
vertices   = cortex.vertices;
faces      = cortex.faces;
elec_pos   = coor;
%%
%% Computing Simulation Dimensions
Nv         = size(vertices,1);    % Number of vertices
Nsubj      = size(LeadFields,2);  % Number of Subjects
if strcmp(sens_system,'pseudo') == 1 || strcmp(sens_system,'small') == 1
    Nsamp      = 600; % sample number
    d0         = 9E0; % patches geodesic radious
elseif strcmp(sens_system,'large') == 1
    Nsamp      = 6000; % sample number
    d0         = 1E-2; % patches geodesic radious
end
db_source  = 0.1;
db_sens    = 0.1;
%%
%%
% The Simulated Sources Configurations are built by marking the coordenates in the
% 3d visualization using function Plot_sources_Haufe:
%     Plot_sources_Haufe(Jtest,vertices,faces,'simple')
% Later the coordinates are find in the vertices vector:
% Ej. coordinate = (19.53,72.34,32.39), then:
%     ans = find(abs(vertices(:,2)-19.53)<0.01)
%% Simulated Sources Locations at Cortical Surface Mesh
%%  Definition of Cortical points
if strcmp(sens_system,'large') == 1
    load('data_tips24_large.mat')
elseif strcmp(sens_system,'small') == 1
    load('data_tips22_small.mat')
elseif strcmp(sens_system,'pseudo') == 1
    load('data_tips22_pseudo.mat')
end
%%
Nsim         = 100;
Nseed        = size(data_tips,2);
if strcmp(sens_system,'large') == 1 || strcmp(sens_system,'small') == 1
    Seeders      = zeros(1,Nseed);
    for cont = 1:Nseed
        coord_tmp = data_tips(cont).Position;
        vx            = coord_tmp(1);
        vy            = coord_tmp(2);
        vz            = coord_tmp(3);
        Seeders(cont) = pickpoint(vx,vy,vz,vertices,1E-3);
    end
else
    Seeders      = 1:Nseed;
end
Seeders      = sort(Seeders);
%%
%% gen_hgggm options
options.config      = 2;
options.var         = 2;
options.connections = [1 2; 2 3];
options.extensions  = [ceil(Nseed/3); ceil(Nseed/3); Nseed - 2*ceil(Nseed/3)];
%%
index_seed  = cell(1,Nseed);
index_full  = [];
for point = 1:Nseed
    if strcmp(sens_system,'large') == 1 || strcmp(sens_system,'small') == 1
        Source            = Seeders(point);
        [index,findex]    = surfpatch(Source,vertices,faces,d0);
        index_seed{point} = index;
        index_full        = [index_full; index];
    else
        index_seed{point} = point;
        index_full        = [index_full; point];
    end
end
Nnoise      = length(index_full);
end