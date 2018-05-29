function [Ethreshold] = OptimNbinsThresholds(specdat, nbins, varargin)
%  function [idx,Ethreshold,specnum] = OptimNbinsThresholds(specdat,nbins,varargin)
%  Determine optimal thresholds for given spectrum, according to the
%  principal of equal number of photons in each energy bin, at the
%  detector (after passing through a body).
%  
%  inputs:
%     specdat: an array giving specnum or a spectrum structure with fields: specnum and egys
%     nbins: the number of bins
%     zA: (optional = 0) additional attenuator. given as ('a_space', [array of
%     Aspace values])
%  outputs:
%     idx: an (nbins-1) integer array of the indexes into the energies array for the optimal thresholds
%     Ethreshold: the corresponding energies
%     
%  REA 2/13-8/23/10

  % process the optional arguments
nreqargs = 2;
assert(nargin>=nreqargs);

do_attenuate = 0;
attenuated_specnum = 0;

if(nargin>nreqargs)
  i=1;
  while(i<=size(varargin,2))
     switch lower(varargin{i})
%     case 'max'
%         minmax=1;
     case 'a_space'
         a_coor=varargin{i+1};
         do_attenuate = 1;
         i=i+1;
     otherwise
        error('Unknown argument %s given',varargin{i});
     end
     i=i+1;
  end
end

assert( isstruct(specdat) );
    % attenuate the spectrum if zmus present
if do_attenuate    
    Transmission_E = exp(-specdat.mus*a_coor(:));
    attenuated_specnum = specdat.specnum_N0(:).*Transmission_E(:);
    specnum = attenuated_specnum;
else
    specnum = specdat.specnum_N0;
end

if sum(specnum) < 1000 % make sure sum(specnum) big enough so can do integer calcs below
    specnum = 1000*specnum/sum(specnum);
end
integral = cumsum(specnum);
idx = zeros(nbins-1,1);
for k = 1:(nbins-1)
   idx(k) = find(integral< floor(k*integral(end)/nbins),1,'last');
end
Ethreshold = specdat.energies(idx);