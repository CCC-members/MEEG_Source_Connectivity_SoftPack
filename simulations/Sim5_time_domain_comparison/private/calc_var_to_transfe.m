% Authors:
% - Ying Wang

% Date: Nov 24, 2019

function [H,K]=calc_var_to_transfe(A,w,approxi)
[n,~,p]=size(A);
fprintf('<calc_var_to_transfe.m> Calculating transfer function from coefficient of variable(%d) order(%d) vector autoregressive model; \n',...
    n,p);
if nargin<3
    approxi=false;
end
gpuInput=isa(A,'gpuArray');
% if gpuInput || test_gpu
if test_gpu
    tic
    A0=gpuArray(full(A));
    A=zeros(size(A0,1),size(A0,2),size(A0,3),'gpuArray');
    H=zeros(size(A0,1),size(A0,2),length(w),'gpuArray');
    K=zeros(size(A0,1),size(A0,2),length(w),'gpuArray');
    k=1:p;
    for j=1:length(w)
        e=exp(-1i*k*w(j));
        if p==1
            A=A0*e;
        else
            for i=1:p
                A(:,:,i)=A0(:,:,i)*e(i);
            end
            A=sum(A,3);
        end
        
        if approxi
            L = chol(A,'lower');
            U = L\eye(n,'gpuArray');
            H(:,:,j) = L'\U;
        else
            K(:,:,j)=A;
            H(:,:,j)=inv(eye(n,'gpuArray')-K(:,:,j));
        end
    end
    if ~gpuInput
        K=gather(K);
        H=gather(H);
    end
    toc
else
    tic
    if issparse(A)
        A0=A;
    else
        A0=full(A);%ndSparse
    end
    A=zeros(size(A0,1),size(A0,2),size(A0,3));
    H=zeros(size(A0,1),size(A0,2),length(w));
    K=zeros(size(A0,1),size(A0,2),length(w));
    k=1:p;
    for j=1:length(w)
        e=exp(-1i*k*w(j));
        if p==1
            A=A0*e;
        else
            for i=1:p
                A(:,:,i)=A0(:,:,i)*e(i);
            end
            A=sum(A,3);
        end
        if approxi
            L = chol(A,'lower');
            U = L\eye(n);
            H(:,:,j) = L'\U;
        else
            K(:,:,j)=A;
            H(:,:,j)=inv(eye(n)-K(:,:,j));
        end
    end
    toc
end



end
