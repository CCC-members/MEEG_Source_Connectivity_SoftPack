function K = apply_tensor_phase(theta0, nu0, deltaf, Fmin, Fmax)
%% emulates an spectral factorization with unitary covariance for the 
%  noise proces Iq, and spectral factors of the central slice (nu0 
% frequency) equal to U0 (this comes fro a hermitic singular value 
% decomposition), then changes the phase of the spectyral factors across 
% frequencies assuming that the underlying MVAR has unique time delays 
% encoded in the phase 
%%
[q,~]         = size(theta0);
[U0,D]        = svd(theta0);
D             = sqrt(D);
U0            = U0*D;
Iq            = eye(q);
K0            = Iq - U0';
pha0          = angle(K0);
F             = Fmin:deltaf:Fmax;
pha_rate      = F/nu0;
pha           = bsxfun(@times,pha0,reshape(pha_rate,1,1,length(F)));
K             = repmat(abs(K0),1,1,length(F)).*exp(1i *(pha));
% K             = repmat(K0,1,1,length(F));

end
