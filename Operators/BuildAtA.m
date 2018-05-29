function [ AtA ] = BuildAtA( n ,DecFactor,DecMethod,OverrideFlag )
%ATA Performs the 2D convolution required for the image filteration

% Building for double size
n = 2*n;

m = 2*n+1;
ThisFilePath = fileparts(mfilename('fullpath'));
FileName = [ThisFilePath '/../ConvMatrices/AtA_',num2str(n/2),'_',DecMethod,num2str(DecFactor),'.mat'];

ThisFilePath = fileparts(mfilename('fullpath')); 
FileName = [ThisFilePath '/../ConvMatrices/AtA_',num2str(n/2),'_',DecMethod,num2str(DecFactor),'.mat']; 

%% Calcaulting the convolution kernel
if nargin < 4
    OverrideFlag = 0;
end

% Check if matrix already exists
if (exist(FileName, 'file') == 2) && (~OverrideFlag)
    hStruct = load(FileName);
    h = hStruct.h;
    % Returning the convlution operator
    AtA =@(x) ApplyAtA(x,h);
    return;
end
    
h = zeros(n);
Dop = DecOperator(DecFactor,DecMethod);
Dmat = Dop(ones(m,2*n+2));

Mvec = M(complex(ones(m,1)));
% Mvec = M_old(complex(ones(m,1)));
% Mvec = ones(2*n+1,1);
% wb = MyWaitbar(0,'Building the AtA matrix');
cpb = ConsoleProgressBar();
cpb.setText('Calculating AtA');  % update user text
cpb.start();
% for u=-n/2:n/2-1
%     for v=-n/2:n/2-1
for u=-n/2:n/2-1
    for v=-n/2:n/2-1
        for k=-n:n
%             for l=-n/2:n/2
            for l=-n/2:n/2
                if Dmat(k+n+1,l+n/2+1)
%                     h(u+n/2+1,v+n/2+1) = h(u+n/2+1,v+n/2+1) ...
%                         + Mvec(k+n+1)*( exp(2j*pi*k/m * (u-2*l*v/n)) ...% );...
%                         + exp(2j*pi*k/m*(v-2*l*u/n)));
%                 end
                    h(u+n/2+1,v+n/2+1) = h(u+n/2+1,v+n/2+1) ...
                        + Mvec(k+n+1)*(exp(2j*pi*k/m * (u-2*l*v/n)) ...% );...
                        + exp(2j*pi*k/m*(v-2*l*u/n)));
                        
                        
                end

            end
        end
    end
    cpb.setValue((u+n/2+1)/n);
%     (u+n/2+1)/n
%     MyWaitbar((u+n/2+1)/n,wb);
end
% close(wb);
% cpb.stop();

% We want only the real part
h = real(h);
save(FileName,'h');

% Returning the convlution operator
AtA =@(x) ApplyAtA(x,h);

end

function [ x_hat ] = ApplyAtA(x,h)
    N = size(x,1);
%     M = size(h,1)
    x_pad  = padarray(x,[N/2 N/2]); 
    x_conv = real(FreqConv2(x_pad,h));
    x_hat  = x_conv(N/2+1:3*N/2 , N/2+1:3*N/2) / (2*N+1)^2;

end


