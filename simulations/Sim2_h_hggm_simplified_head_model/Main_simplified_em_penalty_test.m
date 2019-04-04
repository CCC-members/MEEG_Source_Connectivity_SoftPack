function result = Main_simplified_em_penalty_test(output_sourse)


%% Test Penalties


output_sourse = strcat(output_sourse, filesep,'Simplified_em_penalty_test');
mkdir(output_sourse);


process_waitbar = waitbar(0,'Please wait...');
%%
%% Create partial correlations (ThetaJJ)
%  Generate source (state) empirical covariance (SJJ) and source (state) activity (J)
m                    = 600;             % Sample number
q                    = 22;              % Number of generators
p                    = 30;              % Number of sensors
nblocks              = 2;               % Number of blocks in simulation
options.config       = 2;               % (2) overlapping blocks (1) nonoverlapping blocks
options.var          = 2;               % (2) complex variable (1) real variable
options.extensions   = [ceil(q/3); ceil(q/3); q - 2*ceil(q/3)]; % patches extensions
options.connections  = [1 2; 2 3];      % patches connections
% [Sjj_sim,j_sim,Thetajj_sim] = gen_hggm1(m,q,nblocks,options);
% j_sim = transpose(j_sim);

[Sjj_sim,j_sim,Thetajj_sim] = gen_hggm2(m,q,options);
%% Creating pseudoLead Field (L)
Lvj             = zeros(p,q);
radj            = 60;
radv            = 85;
angj            = 2*pi/q;
angv            = 2*pi/p;
for contv = 1:p
    for contj = 1:q
        waitbar((contv*contj)/(p*q),process_waitbar,strcat('Creating pseudoLead Field (L)'));
        vectv            = [radv*cos((contv-1)*angv); radv*sin((contv-1)*angv)];
        vectj            = [radj*cos((contj-1)*angj); radj*sin((contj-1)*angj)];
        r                = vectv - vectj;
        r_unit           = r/sqrt(sum(abs(r).^2));
        miu              = vectj/sqrt(sum(abs(vectj).^2));
        Lvj(contv,contj) = (1/(4*pi))*miu'*r_unit/sqrt(sum(abs(r).^2))^2;
    end
end
delete(process_waitbar);


% LeadFields      = {Lvj};
% save('LeadFields_pseudo','LeadFields')
%% pseudo-cortex

process_waitbar = waitbar(0,'Please wait...');
vertices        = zeros(q,2);
for contj = 1:q
    waitbar((contj)/(q),process_waitbar,strcat('pseudo-cortex'));
    vertices(contj,:)    = [radj*cos((contj-1)*angj) radj*sin((contj-1)*angj)];
end
delete(process_waitbar);

process_waitbar = waitbar(0,'Please wait...');
faces           = [[1:q]' [2:q 1]'];
cortex.vertices = vertices;
cortex.faces    = faces;
coor            = zeros(p,2);
for contv = 1:p
    waitbar((contv)/(p),process_waitbar,strcat('HeadModel-pseudo ','cortex ','coor'));
    coor(contv,:)    = [radv*cos((contv-1)*angv) radv*sin((contv-1)*angv)];
end
save(strcat('simulations',filesep,'Sim2_h_hggm_simplified_head_model',filesep,'HeadModel_pseudo'),'cortex','coor')
delete(process_waitbar);
%%

process_waitbar = waitbar(0,'Please wait...');

[U,D,V]         = svd(Lvj,'econ');
dmin            = min(diag(D));
%% Generate data
v0              = Lvj*j_sim; % data (observation)
%% Biological noise
bionoise        = randn(q,m) + 1i*randn(q,m);
bionoise        = Lvj*bionoise;
bionoise        = sum(abs(v0(:)).^2)^(1/2)*bionoise/sum(abs(bionoise(:)).^2)^(1/2);
%% Sensor noise
sensnoise       = randn(p,m) + 1i*randn(p,m);
sensnoise       = sum(abs(v0(:)).^2)^(1/2)*sensnoise/sum(abs(sensnoise(:)).^2)^(1/2);
%% Corrupted data
v               = v0 + 0.1*bionoise + 0.1*sensnoise;
%% Data empirical covariance
Svv             = cov(v');
%% Likelihood test
param.maxiter_outer = 60;
param.maxiter_inner = 30;
param.m             = m;
param.rth           = 3.16;
param.axi           = 1E-5;
param.Axixi         = eye(length(Svv));
%%
penalty             = [1 2 0]; % 1 (lasso) 2 (frobenious) 0 (naive)
Thetajj_est         = zeros(q,q,length(penalty));
llh_outer           = cell(1,length(penalty));
llh_inner           = cell(1,length(penalty));
for cont = 1:length(penalty)
    waitbar((cont)/(length(penalty)),process_waitbar,strcat('Likelihood test'));
    param.penalty  = penalty(cont);
    [Thetajj_est(:,:,cont),Sjj_est,llh_outer{cont},xixi_on,jj_on] = h_hggm(Svv,Lvj,param);
end
delete(process_waitbar);
%% Plot Results
figure_partial_correlation_maps = figure;
%%
X = Thetajj_sim;
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(2,2,1); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('simulated partial correlations')
%%
X = Thetajj_est(:,:,1);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(2,2,2); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('h-hggm-lasso partial correlations')
%%
X = Thetajj_est(:,:,2);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(2,2,3); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('h-hggm-ridge partial correlations')
%%
X = Thetajj_est(:,:,3);
X = X - diag(diag(X));
X = X/max(abs(X(:)));
subplot(2,2,4); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('h-hggm-naive partial correlations')
colormap('hot');

saveas( figure_partial_correlation_maps,strcat(output_sourse,filesep,'partial_correlation_maps.fig'));
disp(strcat('Saving figure ---->  Partial Correlation Maps to  ---> ', output_sourse) );
delete(figure_partial_correlation_maps);

end