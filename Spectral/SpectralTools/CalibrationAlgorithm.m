function [ Q_MLE,  a_Error, as_init_fit,M_calib] = CalibrationAlgorithm( calibration_data )
%CALIBRATIONALGORITHM Summary of this function goes here
%   Detailed explanation goes here

% covariance from calibrator data in middle of calibration region
% [dim1,dim2,dim3] = size(cdat.As);
% distance_from_center = cdat.As-repmat(shiftdim(cdat.a_calib_region/2,-1),dim1,dim2,1);
% [~, idx_noise_dot] = min(sqrt(sum(distance_from_center.^2, 3)));
% C_L = diag(1./cdat.ns(idx_noise_dot,:)); %FOR US: check if optimal 
% M_calib = (cdat.As\cdat.p_k)'; % note transpose
% C_Li = inv(C_L);
% Q_MLE = inv(M_calib'*C_Li*M_calib)*M_calib'*C_Li;
% as_init_fit = cdat.p_k*(Q_MLE');
% a_Error = cdat.As - as_init_fit;
middle = calibration_data.a_calib_region/2;
[adim1,adim2,adim3] = size(calibration_data.As);
as_matrix = reshape(calibration_data.As,adim1*adim2,adim3);
[pdim1,pdim2,pdim3] = size(calibration_data.p_k);
ps_matrix = reshape(calibration_data.p_k,pdim1*pdim2,pdim3);
ns_matrix = reshape(calibration_data.ns,pdim1*pdim2,pdim3);

sub_as = as_matrix-repmat(middle,adim1*adim2,1);
distance = sqrt(sum(abs(sub_as).^2,2));
[~,ind_min] = min(distance);

C_L = diag(1./ns_matrix(ind_min,:));
C_Li = inv(C_L);
M_calib = (as_matrix\ps_matrix)';
Q_MLE = pinv(M_calib'*C_Li*M_calib)*M_calib'*C_Li;
as_init = ps_matrix*(Q_MLE');
a_Error_matrix = as_matrix-as_init;
a_Error = reshape(shiftdim(a_Error_matrix,-1),adim1,adim2,adim3);
as_init_fit = reshape(shiftdim(as_init,-1),adim1,adim2,adim3);

