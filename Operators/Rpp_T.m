function [ x ] = Rpp_T( y )
%RPP_T Applies the adjoint  pseudo-polar radon transform on y to get x
% This operates directly on the sinogram, and applies a preconditioner too

x = real(App_T(F1(y)));

end

