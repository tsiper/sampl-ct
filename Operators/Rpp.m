function [ y ] = Rpp( x ,varargin )
%RPP Applies the pseudo-polar radon transform on the input x, to generate a
% pseudo-polar sinogram y

y = real(invF1(App(x,varargin{:})));

end

