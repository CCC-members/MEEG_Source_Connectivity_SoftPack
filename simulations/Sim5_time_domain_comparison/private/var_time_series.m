function Y=var_time_series(A,Y,RES)
% Authors:
% - Ying Wang

% Date: Nov 24, 2019
%% VAR (vector autoregressive model) 
% model: 
% Y_t=c+A_1*Y_t-1+A_2*Y_t-2+...A_p*Y_t-p+��t
% A: n*n*p
% Y: n*m*t
% RES: n*m*t

% n: num of viriable
% m: num of trial
% t: num of sample points
[n,~,p]=size(A);
[~,m,t]=size(Y);
% A=sparse(reshape(A,[n,n*p]));
% A=sparse(reshape(flip(A,3),[n,n*p]));
A=reshape(flip(A,3),[n,n*p]);
% Y=sparse(reshape(permute(Y,[1,3,2]),n*t,m));
Y=reshape(permute(Y,[1,3,2]),n*t,m);

for i=1:t
    if i-p<=0
        Y(n*(i-1)+1:n*i,:)=A(:,1:n*(i-1))*Y(n*(i-1)-n*(i-1)+1:n*(i-1),:)+RES(:,:,i);
        continue
    end
    Y(n*(i-1)+1:n*i,:)=A*Y(n*(i-1)-n*p+1:n*(i-1),:)+RES(:,:,i);
end
Y=permute(reshape(Y,n,t,m),[1,3,2]);
end
