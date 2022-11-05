function KL = map_kullback(P,Q)
KL  = 0.5*(P/Q + P\Q);
KL  = abs(KL - diag(diag(KL)));