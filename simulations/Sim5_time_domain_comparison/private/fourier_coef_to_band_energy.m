function energy=fourier_coef_to_band_energy(fourier_coef,band,deltaf)
f=[round(band(1)/deltaf)+1:1:round(band(2)/deltaf)+1];
[n,N,m]=size(fourier_coef);
energy=zeros(n,length(f));
% energy=zeros(m,1);
for i=1:n
energy(i,:)=sum(abs(fourier_coef(i,f,:)),3)/m;
end
end


