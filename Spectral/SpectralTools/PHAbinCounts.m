function [p_k_picture, ns4bins, N0] = PHAbinCounts(As_picture,specdat,idx_threshold) 
% function ns4bins = PHAbinCounts(As,specdat,idx_threshold) 
% PHA counts for each bin for each value of As
% inputs:
%  As:  an nthickness by nbasis column matrix 
%  specdat: a spectrum data struct for incident spectrum with fields
%        energies: x-ray energies
%        specnum: the photon number spectrum
%        mus: the basis materials attenuation coefficients at the energies. 
%        idx_threshold: indices into the spectrum array for inter-bin energies
% outputs:
%   ns4bins: the total photon counts for each bin along columns for A values along rows
%   N0 = initial photons at each energy bin, before transmission
%   p_k = data at detector
% REA 4/20/10 - 8/Sep/2013 10:04
% see also PoissonLambda. This is a simplified and faster version

nreqargs = 3;
assert(nargin >= nreqargs);
if nargin ==2
    idx_threshold = [];
end
% reshaping radon materials to a vector
[dim1,dim2,dim3] = size(As_picture);
As = shiftdim(reshape(As_picture,1,dim1*dim2,dim3),1);
nthick = dim1*dim2;
nbins = numel(idx_threshold)+1;       

bin_edges = [0; idx_threshold(:);numel(specdat.energies)];
bin_idxs = cell(nbins,1);
for k = 1:nbins
    bin_idxs{k}= (bin_edges(k)+1):bin_edges(k+1);
end

transmission = exp(-As*(specdat.mus')); % size = (nthick,nenergy)

% ns_spectrum is the number of photons for each calibration point at each
% energy
N_detector = bsxfun(@times,transmission,specdat.specnum_N0(:)'); 
N_detector = round(N_detector);
ns4bins = zeros(nthick,nbins);
N0 = zeros(1, nbins);
for kbin = 1:nbins
    ns4bins(:,kbin) = sum(N_detector(:,bin_idxs{kbin}),2);
    N0(kbin) = sum(specdat.specnum_N0(bin_idxs{kbin}));
end
p_k = -log(bsxfun(@rdivide,ns4bins,N0)+eps);%to prevent log(0)
ns4bins = reshape(shiftdim(ns4bins,-1),dim1,dim2,nbins);
p_k_picture = reshape(shiftdim(p_k,-1),dim1,dim2,nbins);