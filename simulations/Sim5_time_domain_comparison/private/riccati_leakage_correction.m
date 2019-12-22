function [Sigma]=riccati_leakage_correction(Lvj,Tjv,Sjj,param)

R            = Tjv*Lvj;
aj           = param.aj;
% [V1,dr,Uconj] = svd(R);
% [V2,ds,Vconj] = svd(Sjj);
[U,V,X,C,S] = gsvd(R',Sjj');
% MATLAB
% A   = U*C*X'
% B   = V*S*X'
% C'*C + S'*S = I 

% Our case
% R   = V*diag(dr)*U';
% Sjj = V*diag(dr)*V';

% R'  = U*diag(dr)'*V' = U*C*X';
% Sjj'= V*diag(ds)'*V' = V*S*X';

% R   = V*diag(dr)*U'  = X*C*U';
% Sjj = V*diag(dr)*V'  = X*S*V';

% V=X; % V should equal to X
% U=U;
% dr=diag(C);
% ds=diag(S);

%% check
% R1=	X*C'*U';
% Sjj1=X*S'*V';
% figure(1); 
% imagesc(abs(R1-R));colorbar
% figure(2); 
% imagesc(abs(Sjj1-Sjj));colorbar
%%

% for i=1:length(dr)
%     dsigma(i,:)=roots([dr(i)^2,-ds(i),-0^2]);
% end
% sigma1=U*diag(dsigma(:,1))*U';
% Sjj=U*diag(dsigma(:,2))*U';

d = abs(diag(X'*V)).*diag(S);
c = abs(diag(C));
dsigma=2*aj^2./((-d.*c.^2)+sqrt((d.^2).*(c.^4)+4*(aj^2)*(c.^4)));
Sigma=U*diag(dsigma)*U';


end

