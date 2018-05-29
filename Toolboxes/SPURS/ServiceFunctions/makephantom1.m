function [ sqrmat ] = makephantom1(N,percent)
if percent == 0
    factor = N/256;
else
    factor = 1/256;
end

offset = 0;
scale = 4;
    
sqrmat = zeros(6,11);
%  Return square phantom data.
sqrmat(1,1) = (2 + offset)/scale	;% amplitude
sqrmat(2,1) = 80*factor	;% x-length
sqrmat(3,1) = 80*factor	;% y-length
sqrmat(4,1) = 0*factor		;% x-shift
sqrmat(5,1) = 0*factor		;% y-shift
sqrmat(6,1) = 0		;% rotation angle

sqrmat(1,2) = (2 + offset)/scale	;% amplitude
sqrmat(2,2) = 2*factor		;% x-length
sqrmat(3,2) = 20*factor	;% y-length
sqrmat(4,2) = -30*2*factor	;% x-shift
sqrmat(5,2) = 20*2*factor	;% y-shift
sqrmat(6,2) = 0		;% rotation angle

sqrmat(1,3) = (2 + offset)/scale	;% amplitude
sqrmat(2,3) = 2*factor		;% x-length
sqrmat(3,3) = 20*factor	;% y-length
sqrmat(4,3) = -25*2*factor	;% x-shift
sqrmat(5,3) = 20*2*factor	;% y-shift
sqrmat(6,3) = 0		;% rotation angle

sqrmat(1,4) = (-1 + offset)/scale	;% amplitude
sqrmat(2,4) = 8*factor		;% x-length
sqrmat(3,4) = 20*factor	;% y-length
sqrmat(4,4) = -15*2*factor	;% x-shift
sqrmat(5,4) = 20*2*factor	;% y-shift
sqrmat(6,4) = 0		;% rotation angle

sqrmat(1,5) = (1 + offset)/scale	;% amplitude
sqrmat(2,5) = 14*factor	;% x-length
sqrmat(3,5) = 20*factor	;% y-length
sqrmat(4,5) = 2*2*factor		;% x-shift
sqrmat(5,5) = 20*2*factor	;% y-shift
sqrmat(6,5) = 0		;% rotation angle

sqrmat(1,6) = (2 + offset)/scale	;% amplitude
sqrmat(2,6) = 20*factor	;% x-length
sqrmat(3,6) = 20*factor	;% y-length
sqrmat(4,6) = 30*2*factor	;% x-shift
sqrmat(5,6) = 20*2*factor	;% y-shift
sqrmat(6,6) = 0		;% rotation angle

sqrmat(1,7) = (2 + offset)/scale	;% amplitude
sqrmat(2,7) = 2*factor		;% x-length
sqrmat(3,7) = 20*factor	;% y-length
sqrmat(4,7) = -35*2*factor	;% x-shift
sqrmat(5,7) = 20*2*factor	;% y-shift
sqrmat(6,7) = 0		;% rotation angle

sqrmat(1,8) = (2 + offset)/scale	;% amplitude
sqrmat(2,8) = 30*factor	;% x-length
sqrmat(3,8) = 2*factor		;% y-length
sqrmat(4,8) = -20*2*factor	;% x-shift
sqrmat(5,8) = -18*2*factor	;% y-shift
sqrmat(6,8) = 0		;% rotation angle

sqrmat(1,9) = (2 + offset)/scale	;% amplitude
sqrmat(2,9) = 18*factor	;% x-length
sqrmat(3,9) = 12*factor	;% y-length
sqrmat(4,9) = 5*2*factor		;% x-shift
sqrmat(5,9) = -18*2*factor	;% y-shift
sqrmat(6,9) = 0		;% rotation angle

sqrmat(1,10) = (1.5 + offset)/scale	;% amplitude
sqrmat(2,10) = 12*factor	;% x-length
sqrmat(3,10) = 20*factor	;% y-length
sqrmat(4,10) = 21*2*factor	;% x-shift
sqrmat(5,10) = -18*2*factor	;% y-shift
sqrmat(6,10) = 0		;% rotation angle

sqrmat(1,11) = (-1 + offset)/scale	;% amplitude
sqrmat(2,11) = 8*factor	;% x-length
sqrmat(3,11) = 32*factor	;% y-length
sqrmat(4,11) = 32*2*factor	;% x-shift
sqrmat(5,11) = -18*2*factor	;% y-shift
sqrmat(6,11) = 0	;% rotation angle

end

