function y=calc_band_pass(y,band,Fs)
if nargin<2
    band = [9,11];
end
if nargin<3
    Fs=200;
end

[m,n,s]=size(y);%m:variabel,n:length,s:segments or trials
y=permute(y,[1,3,2]);
y=reshape(y,[m*s,n]);
y= osl_filter(y,band,'fs',Fs);
y=reshape(y,[m,s,n]);
y=permute(y,[1,3,2]);
end
