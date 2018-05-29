function [PQ,PZ]=cfrftV3_precomp(m,alpha)
%
% [PQ,PZ]=cfrftV3_precomp(m,alpha)
%
% Generate tables to be used by cfrftV3. 
% The tables contain precomputed factors that are for computing the
% fractional Fourier transform. Avoiding recomputing these factors speeds
% up the computation.
%
% Input
%   alpha    A vector with the spacings with which the fractional Fourier
%            transform will be computed in subsequent calls to cfrftV3.
%            After calling the current function, instead of calling
%            w=cfrftV3(x,alpha), where alpha is the required spacing, call
%            w=cfrftV3(x,0,PQ,PZ,k), where PQ and PZ are the tables
%            computed by the current function. The fractional Fourier 
%            transform will be then computed with spacing alpha(k).
%   m        The length of the signal that will be tramsformed by
%            subsequent calls to cfrftV3.
%
% Returns the precomputed tables PQ and PZ.
%
% Revised:
% Yoel Shkolnisky  May 17, 2010.

n=length(alpha);

lm= -fix(m/2);
hm= fix((m-0.5)/2);
ofs=floor(3*m/2)+1;

j=lm:hm;
j2= -m:m;

PQ=zeros(m,n);
PZ=zeros(3*m,n);

for k=1:n    
    E=1i*pi*alpha(k);
    PQ(:,k)=exp(-E*j.^2/m);

    z=zeros(1,3*m);
    z(-m+ofs:-m+ofs+length(j2)-1)=exp(E*j2.^2/m);
    PZ(:,k)=fft(z);    
end