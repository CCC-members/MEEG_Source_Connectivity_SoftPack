function y=calc_orthogonalize(y,protocol)
%   orthogonalisation methods for source-spread correction. It can be:
%     'none'          - No correction applied. 
%     'symmetric'     - Apply orthogonalisation on the parcel time-courses.
%                       This produces orthonormal parcel time-courses
%                       which are as close as possible to the original
%                       time-courses.
%     'closest'       - Apply orthogonalisation on the parcel time-courses.
%                       Start as for the symmetric method, then converge to
%                       a (not orthonormal) orthogonal matrix which is as
%                       close as possible to the original time-courses. 
%     'householder'   - Orthogonalise using a more numerically stable
%                       alternative to the Gram-Schmidt process, dealing
%                       with ROI time-courses in a random order. 
%
if nargin<2
    protocol='closest';
end


[m,n,s]=size(y);%m:variabel,n:length,s:segments or trials
% y=permute(y,[1,3,2]);
% y=reshape(y,[m*s,n]);
% y= osl_filter(y,band,'fs',Fs);
% y=reshape(y,[m,s,n]);
% y=permute(y,[1,3,2]);
% y = remove_source_leakage(y, protocol);
% y = arrayfun(@(i) ROInets.remove_source_leakage(y(:,:,i),protocol), 1:s,'UniformOutput',false);
% y =reshape(cell2mat(arrayfun(@(i) ROInets.remove_source_leakage(y(:,:,i),protocol), 1:s,'UniformOutput',false)),m,n,s);
for i=1:s
y(:,:,i)=ROInets.remove_source_leakage(y(:,:,i),protocol);
end
end

%% 
% closed two warning in [L, d, rho, W] = closest_orthogonal_matrix(A)


