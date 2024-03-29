% KLS = Kullback_Leibler_sym(P,Q)
% Kullback_leibler distance.

function KL = distance_kullback(P,Q)
KL  = 0.5*trace(P/Q + P\Q - 2*eye(length(P)));