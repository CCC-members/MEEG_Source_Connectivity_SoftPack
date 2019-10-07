function [llh] = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_pst,param)
axi               = param.axi;
aj                = param.aj;
Ajj               = param.Ajj;
penalty           = param.penalty;
%%
[U,D]             = eig(Sigmajj_pst);
d_pst             = abs(diag(D));
[U,D]             = eig(Thetajj);
d                 = abs(diag(D));
llh_pst           = sum(log(d_pst));
% llh_pst           = 0;
llh_xixi          = sum(log(1./sigma2xi)) - sum(abs(diag(diag(1./sigma2xi)*Sxixi))) - axi*sum(1./sigma2xi);
if penalty == 0 % naive
    llh_jj                = sum(log(d)) - sum(abs(diag(Thetajj*Sjj)));
elseif penalty == 1 % lasso
    llh_jj                = sum(log(d)) - sum(abs(diag(Thetajj*Sjj))) - aj*sum(abs(Ajj(:).*Thetajj(:)));
elseif penalty == 2 % ridge
    llh_jj                = sum(log(d)) - sum(abs(diag(Thetajj*Sjj))) - (aj^2/2)*sum(abs(diag(Thetajj*Thetajj')));
end
llh               = llh_pst + llh_xixi + llh_jj;
end