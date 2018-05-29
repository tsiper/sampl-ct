% function rr=slowradon3(im);
%
% Compute the 3D pseudo-Radon transform directly, according to its definition.
% The computation requires O(n^6) operations for a nxnxn image.
%
% Input:
%   im - 3D array of size nxnxn.
%        first  index - x direction
%        second index - y direction
%        third  index - z direction
%     In each of the directions, the range 1...n is treated as -n/2...n/2-1.
%     The input image should be of size nxnxn with even n.
%
% Output:
%   xr - 4D array 3x(n+1)x(n+1)x(3n+1) containing the 3-D discrete Radon transform of the input image.
%     The output array has the following format:
%        first  index  - 1 for x-planes, 2 for y-planes, 3 for z-planes
%        second index  - plane intercept
%        third  index  - slope p
%        fourth index  - slope q
%
% Yoel Shkolnisky 29/01/03

function rr=slowradon3(im)

% verify that the input is a 3D image of size nxnxn
verifyImage(im);

% Initialize output data structure
s=size(im);
n=s(1); % at this point n is even
m = 3*n+1;
rr = zeros(3,3*n+1,n+1,n+1);

% Compute the 3-D discrete Radon transform for x-planes
for p=-n/2:n/2  %first slope
   s1=2*p/n;
   
    for q=-n/2:n/2 %second slope
      s2=2*q/n;
      
      for t=-3*n/2:3*n/2    %intercept
        acc=0;
        for v=-n/2:n/2-1
            for w=-n/2:n/2-1
               acc = acc+I1(im,n,s1*v+s2*w+t,v,w);
            end
        end
        coord = toUnaliasedCoord([t,p,q],[m,n+1,n+1]);
        rr(1,coord{:}) = acc;
      end
    end
end

% Compute the 3-D discrete Radon transform for y-planes
for p=-n/2:n/2  %first slope
   s1=2*p/n;
   
    for q=-n/2:n/2 %second slope
      s2=2*q/n;
      
      for t=-3*n/2:3*n/2    %intercept
        acc=0;
        for u=-n/2:n/2-1
            for w=-n/2:n/2-1
               acc = acc+I2(im,n,u,s1*u+s2*w+t,w);
            end
        end
        coord = toUnaliasedCoord([t,p,q],[m,n+1,n+1]);
        rr(2,coord{:}) = acc;
      end
    end
end

% Compute the 3-D discrete Radon transform for z-planes
for p=-n/2:n/2  %first slope
   s1=2*p/n;
   
    for q=-n/2:n/2 %second slope
      s2=2*q/n;
      
      for t=-3*n/2:3*n/2    %intercept
        acc=0;
        for u=-n/2:n/2-1
            for v=-n/2:n/2-1
               acc = acc+I3(im,n,u,v,s1*u+s2*v+t);
            end
        end
        coord = toUnaliasedCoord([t,p,q],[m,n+1,n+1]);
        rr(3,coord{:}) = acc;
      end
    end
end


%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%

% computation of I1 - continuous extension of I along the x direction
function s = I1(im,n,x,v,w)
m=3*n+1;
acc = 0;

for u=-n/2:n/2-1
   coord = toUnaliasedCoord([u,v,w],[n,n,n]);
   acc = acc + im(coord{:})*dirichlet(x-u,m);
end
s=acc; 


% computation of I2 - continuous extension of I along the y direction
function s = I2(im,n,u,y,w)
m=3*n+1;
acc = 0;

for v=-n/2:n/2-1
   coord = toUnaliasedCoord([u,v,w],[n,n,n]);
   acc = acc + im(coord{:})*dirichlet(y-v,m);
end
s=acc; 


% computation of I3 - continuous extension of I along the z direction
function s=I3(im,n,u,v,z)
m=3*n+1;
acc = 0;

for w=-n/2:n/2-1
   coord = toUnaliasedCoord([u,v,w],[n,n,n]);
   acc = acc + im(coord{:})*dirichlet(z-w,m);
end
s=acc; 
