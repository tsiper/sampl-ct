function [ I ] = invF2( Ihat )
%F2 Receives an image I and makes a 2D dft on it using fft

I = invF1( invF1( Ihat.' ).' );

end