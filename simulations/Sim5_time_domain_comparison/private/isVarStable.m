function varargout=isVarStable(varargin)
%% test all modulus of eigenvalus of companion matrix is less then one
A=varargin{1};
[~,n,p]=size(A);
% tensor to matrix
F=reshape(A,[n,n*p]);
I=eye((p-1)*n);
I=[I,zeros((p-1)*n,n)];
F=[F;I];
% tic
modu=eigs(F,1,'LM');
if isnan(modu)
modu=max(eig(F));
end
modu=abs(modu);
% toc
if modu(1)<=1%<1
    varStablity=true;
else
    varStablity=false;
end

switch nargout
    case 0
        varargout{1}=varStablity;
    case 1
        varargout{1}=varStablity;
    case 2
        varargout{1}=varStablity;
        varargout{2}=modu;
    case 3
        varargout{1}=varStablity;
        varargout{2}=D;
        varargout{3}=V;
end


% the companion matrix have modulus less than one
% the stability may related with Gerschgorin theory, where we can selecte
% our A matrices
