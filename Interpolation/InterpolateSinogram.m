function [ ppSinogram ] = InterpolateSinogram( Sinogram, N, theta,InterpMethod)
%INTERPOLATESPLINE 
% Interpolates a sinogram to the PP paradigm
% Input:
%     Sinogram     - The input regular radon parallel-beam sinogram 
%     N            - The resolution of our objective
%     theta        - The vector describing the acquisition angle for each column
%                    in the sinogram
%     InterpMethod - can be 'spline','cubic','pchip','linear','nearest' or 'filt'
%                

% InterpMethod 
if strcmpi(InterpMethod,'filt')
    % Running my subspace based interpolation
    ppSinogram = InterpolateFast( Sinogram, N, theta);
    
else % Running on of the regular approaches
    
    % Padding to columns so it complies with the PP standards
    Sinogram = padcols(Sinogram,2*N+1);
    [n,m] = size(Sinogram);
    
    % The Pseudo-Polar scanning angles
    l = -N/2:N/2;
    theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;
    
    y_rd_angle  = zeros(n,length(theta_pp));
    % Interpolating over the angles
    for i=1:size(Sinogram,1)
        % Now using interp1 to span over the theta_pp angles
        y_rd_angle(i,:) = interp1(theta,Sinogram(i,:),theta_pp,InterpMethod,'extrap');
    end
    y_rd_angle = padcols(y_rd_angle,2*N+1);
    ppSinogram = zeros(size(y_rd_angle));
    % Interpolating over the detectors, namely the t's
    for j=1:2*N+2
        t_vec = -(n-1)/2:(n-1)/2;

        d = y_rd_angle(:,j);
        
        
        if (theta_pp(j) >= -45) && (theta_pp(j) <= 45)
            t_pp_vec = (t_vec + 1/2)*cosd((theta_pp(j)));
            %         ppSinogram(:,j) = spline(t_vec,d,t_pp_vec)*abs(cosd(theta_pp(j)));
            ppSinogram(:,j) = interp1(t_vec,d,t_pp_vec,InterpMethod,'extrap')*abs(cosd(theta_pp(j)));
            
        else
            t_pp_vec = (t_vec + 1/2 -2*(theta_pp(j)-45)/90);
            
            t_pp_vec = t_pp_vec*sind((theta_pp(j)));
            %         t_vec = t_vec - abs( mod(theta_pp(j)+45,90) - 45 )/180; %          t_pp_vec = t_pp_vec - abs( mod(theta_pp(j)+45,90) - 45 )/180;
            %         ppSinogram(:,j) = spline(t_vec,(y_rd_angle(:,j)),t_pp_vec)*abs(sind(theta_pp(j)));
            ppSinogram(:,j) = interp1(t_vec,d,t_pp_vec,InterpMethod,'extrap')*abs(sind(theta_pp(j)));
        end
    end
    % y_pp_spline = flipud(y_pp_spline);
end

end

