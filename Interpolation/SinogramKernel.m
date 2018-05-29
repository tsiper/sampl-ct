function [ F ] = SinogramKernel( theta,t,B,R,W )
% Calculates the kernel in over the values of theta and t, if one of the
% variables is a vector than the output is a vector. If both t and theta are
% vectors than a 2D function will be the output

% Making sure that the inputs are vectorized columns
theta_vec = theta(:);
t_vec     = t(:);

% Initializing F
F = zeros(length(t_vec),length(theta_vec));

% Iterating over all detector t positions 
for i=1:length(t_vec) 
    % iterating over all the theta angles spap2
    for j=1:length(theta_vec)
        theta = theta_vec(j);
        t     = t_vec(i);

        % Bulidlng the kernel according to my derivations, taking the limits into account
        if (theta==0)&&(t~=0)
            F(i,j) = (2*t*sin(W*t)*(B+R*W)+2*R*(cos(W*t)-1)) / (pi*t^2);
        elseif (theta==0)&&(t==0)
            F(i,j) = 2/pi * W * (B+R*W/2);
%         elseif (abs(t)+1e-3 > abs(theta*R)) && (abs(t)-1e-3 < abs(theta*R))
        elseif (abs(t) == abs(theta*R))
            F(i,j) = (2*theta*R*W*sin(B*theta)+cos(B*theta)-cos(theta*(B+2*R*W))) / (2*pi*theta^2*R);
        else
            F(i,j) = 1/pi./theta .* ( (cos(W*t-theta*(B+R*W))-cos(B*theta))./(t-theta*R) - ...
                           (cos(W*t+theta*(B+R*W))-cos(B*theta))./(t+theta*R) );
        end
    end
end

end

