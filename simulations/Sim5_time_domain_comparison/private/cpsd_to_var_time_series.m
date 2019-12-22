function process= cpsd_to_var_time_series(sigma,Ntpoints,Nsegments,Fs)


[G,q] = cpsd_to_autocov(sigma);
[A,SIG] = autocov_to_var(G);
%     [Data,E,mtrunc] = var_to_tsdata(A,SIG,options.numTimepoints,options.numTrials,100);
[process,E]=gen_var_time_series(A,SIG,Ntpoints,Nsegments,Fs);
