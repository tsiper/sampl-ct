function [ Dop ] = DecOperator( DecFactor,DecMethod, varargin)
%DECOPERATOR Returns a decimation operator handle, D, that decimates the
%sinogram measurements. 
%
%  Input:
%    DecFactor -    The decimation factor (integer)
%    DecMethod -    The method for decimating:
%                   {'uniform','random','limited','detector'}
%    Mod_Override - (optional) if 1 then we don't mod over the limited
%                   detector
%
%  Output:
%    D - The decimation operator, that can be accessed by D(y)
%
%  -------------------------------------------------------------

Dop = @(y) D(y,DecFactor,DecMethod, varargin{:});


end


function [ y_hat ] = D( y , factor , method , down_freq , mod_override)
%D Decimates the input y
global RandSeed;
global detector_pattern;
[m,n] = size(y);
y_hat = zeros(m,n);

if nargin < 5
    mod_override = 0;
end

% Checking if this is a sinogram in space domain or frequncy domain
RealFlag = isreal(y);

switch method
    case 'uniform'
        % In this case we split for two cases
        y1 = fliplr(y(:,1:end/2));
        y2 = y(:,end/2+1:end);
        y_hat1 = zeros(m,n/2);
        y_hat2 = zeros(m,n/2);
        for j=1:factor:n/2
            y_hat1(:,j) = y1(:,j);
            y_hat2(:,j) = y2(:,j);
        end
        y_hat1 = fliplr(y_hat1);
        y_hat = [y_hat1,y_hat2];
    case 'random'        
        num_of_measurements = round(n/factor);
        % Make sure we have the same random permuatation each time
        rng(RandSeed);
        perm_inds = randperm(n);
        for j=perm_inds(1:num_of_measurements)
            y_hat(:,j) = y(:,j);
        end
    case 'limited'
        y_hat(:,1:n-floor(n/factor)) = y(:,1:n-floor(n/factor));
    case 'detector'
        if ~RealFlag
            y_space = y;
        else
            y_space = y;
        end
        y_filt  = zeros(m,n);
        d = 0;
        tmp=1;
        for j=1:n
            if ~mod_override, d = mod(tmp-1,factor); end
            for i=1:factor:m-factor-d+1
              %  y_filt(i+d:i+d+factor-1,j) = sum(y_space(i+d:i+d+factor-1,j))/factor;
              %  y_filt(i+d:i+d+factor-1,j) = y_space(i+d,j);
                 y_filt(i+d:i+d+factor-1,j) = median(y_space(i+d:i+d+factor-1,j));
            end
%             y_filt = [y_filt, resample(downsample(y_space(:,j),factor),factor,1)]; %#ok<AGROW>
            if mod(j,down_freq)==0
                tmp=tmp+1;
            end
        end
        if m>size(y_filt,1)
            y_filt = [y_filt; zeros(m-size(y_filt,1),size(y_filt,2))];
        end
        y_hat = y_filt;

    case 'pattern algorithm'
            y_space = y;
            y_filt  = zeros(m,n);
            if (size(detector_pattern) == 0)
              [ detector_pattern ] = pattern_algorithm(factor, y_space );
            end
            for i=1:max(detector_pattern(:))
                relevant_find=find(detector_pattern==i)';
                if ~isempty(relevant_find)
                    y_filt(relevant_find)=sum(y_space(relevant_find))/factor;
                end
            end
            for i=1:length(detector_pattern(:))
                if (detector_pattern(i) == 0)
                    y_filt(i) = y_space(i);
                end
            end
            y_hat=y_filt;
end
   
% Assigning the remaining values to our decimated measurements y_hat

end
