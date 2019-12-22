function [A] = gen_NnvarAdjacent(vertices,faces,SeedsIdx)

% Rigel, Oct 2019

% D = squareform(pdist(double(vertices)));
% D = squareform(pdist(vertices));
if size(faces,2)==3
TR = triangulation(double(faces),double(vertices));
E = edges(TR);
elseif size(faces,2)==2
    E=faces;
end
E=cat(1,E,flip(E,2));

ampDiag=0.9;%0.04 0.9
A=ampDiag*eye(size(vertices,1));

ampOffdiag=-0.04;%-0.45 -0.04

% D(logical(eye(size(D))))=max(max(D));
% [~,I]=min(D,[],2);
A(sub2ind(size(A),E(:,1),E(:,2)))=ampOffdiag;

if nargin>2 && ~isempty(SeedsIdx)
    A=A(SeedsIdx,SeedsIdx);
end

A=sparse(A);



%% STABLIZE
n=size(A,1);
% dmin                     = min(eig(S));
dmax=eigs(A,1,'LM');
if abs(dmax)>1
    if dmax < 0
        A  = A + (abs(dmax)-0.9)*speye(n);
    else
        A  = A - (abs(dmax)-0.9)*speye(n);
    end
end
dmax=eigs(A,1,'LM');

assert(abs(dmax)<1,'Coefficient matrix for VAR(1) not stable,adjust the diagonal or off-diagonal')

end









