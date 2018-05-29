function Psi = scalesShearsAndSpectra ( l, varargin)
%SCALESSHEARSANDSPECTRA compute shearlet spectra
% Compute the shearlet spectra of a given size l. The number of scales
% and a boolean indicating real or complex shearlets are optional
% parameters.
% Using a parameter value list further details can be provided.	
% The output Psi is a 3-d-matrix with the shearlets in Fourier domain
% ordered with ascending scale and within each scale ordered by the
% direction of the respective shears (see comments below for further details).
%
% INPUT:
%  l				(vector) dimensions of image
%  numOfScales		(int) number of scales OR
%	 				(3-d-matrix) precomputed Psi (optional)
%                   Note that for internal use the first element of
%                   varargin can also be a struct passed from
%                   shearletTransformSpect.m.
%  realCoefficients (bool) real/complex shearlets  (optional)
%
% OUTPUT:
%  Psi				(3-d-matrix) spectrum of shearlets
%
% PARAMETERS: (as optional parameter value list, arbitrary order)
%  'shearletSpect'  (string or function handle) shearlet spectrum
%  'shearletArg'	(arbitrary) further parameters for shearlet
%  'realReal'       (bool) guarantees really real shearlets
%  'maxScale'       ('max','min') maximal or minimal finest scale
%
% EXAMPLES:
%  see shearletTransformSpect.m
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

%% parse input
	if( nargin == 1 || ~isstruct(varargin{1}) )
		p = inputParser;
		
		addRequired(p,'l',@checkLength);
		addOptional(p,'numOfScales',defaultNumberOfScales(size(A)),@checkNumOfScales);
		addOptional(p,'realCoefficients',1,@isnumeric);
		p = parseShearletParameterInputs( p );
		
		parse(p,l,varargin{:});
		pp = p.Results;
	elseif(isstruct(varargin{1}))
		pp = varargin{1};
	else
		error('Wrong input arguments!');
    end
    
    if(isempty(pp.numOfScales))
        pp.numOfScales = defaultNumberOfScales(pp.l);
    end
		
%% rectangular images
    if( pp.l(2) ~= pp.l(1) )
        rectangular = 1;
    else
        rectangular = 0;
    end

    %% computation of the shearlets
    
    %for better symmetry each l should be odd
    l_orig = pp.l;
	lm = logical(1-mod(pp.l,2)); %1 for even, 0 for odd
	pp.l(lm) =  pp.l(lm) + 1;
	
	% create meshgrid
	%largest value where psi_1 is equal to 1
	if(strcmp(pp.maxScale,'max'))
		X = 2^(2*(pp.numOfScales-1)+1); % = 2^(2*numOfScales - 1)
	elseif(strcmp(pp.maxScale,'min'))
		X = 2^(2*(pp.numOfScales-1)); % = 2^(2*numOfScales - 2)
	else
		error('Wrong option for maxScale, must be "max" or "min"');
	end
	%xi_x = linspace(-X,X-1/l(2)*2*X,l(2)); %not exactly symmetric
    xi_x_init = linspace(0,X,(pp.l(2)+1)/2);
	xi_x_init = [-fliplr(xi_x_init(2:end)) xi_x_init];
	if(rectangular)
        xi_y_init = linspace(0,X,(pp.l(1)+1)/2);
        xi_y_init = [-fliplr(xi_y_init(2:end)) xi_y_init];
    else
        xi_y_init = xi_x_init;
	end
	%create grid, from left to right, bottom to top
	[xi_x, xi_y] = meshgrid(xi_x_init,fliplr(xi_y_init));

    %cones
	C_hor = abs(xi_x) >= abs(xi_y); %with diag
	C_ver = abs(xi_x) < abs(xi_y);
    
	% number of shears: |-2^j,...,0,...,2^j| = 2 * 2^j + 1 
	% now: inner shears for both cones: 
	% |-(2^j-1),...,0,...,2^j-1| 
	% = 2 * (2^j - 1) + 1
	% = 2^(j+1) - 2 + 1 = 2^(j+1) - 1
	% outer scales: 2 ("one" for each cone)
	% shears for each scale: hor: 2^(j+1) - 1, ver: 2^(j+1) - 1, diag: 2
	%  -> hor + ver + diag = 2*(2^(j+1) - 1) +2 = 2^(j + 2)
	%  + 1 for low-pass
	shearsPerScale = 2.^((0:pp.numOfScales-1)+2);
	numOfAllShears = 1 + sum(shearsPerScale);

	%init
	Psi = zeros([pp.l, numOfAllShears]);
	% frequency domain:
	% k  2^j 0 -2^j
	%
	%     4  3  2  -2^j
	%      \ | /
	%   (5)- x -1  0
	%      / | \
	%              2^j
	%
	%        [0:-1:-2^j][-2^j:1:2^j][2^j:-1:1] (not 0)
	%           hor          ver        hor
	%
	% start with shear -2^j (insert in index 2^j+1 (with transposed
	% added)) then continue with increasing scale. Save to index 2^j+1 +- k,
	% if + k save transposed. If shear 0 is reached save -k starting from
	% the end (thus modulo). For + k just continue.
	% 
	% then in time domain:
	%
	%  2  1  8
	%   \ | /
	%  3- x -7
	%   / | \
	%  4  5  6
	%

	%lowpass
    Psi(:,:,1) = pp.shearletSpect(xi_x, xi_y, NaN, NaN, pp.realCoefficients, pp.shearletArg, 'scaling');

	%loop for each scale
	for j = 0:pp.numOfScales-1
		%starting index
		idx = 2^j + 1;
		start_index = 1 + sum(shearsPerScale(1:j));
		shift = 1;
		for k = -2^j:1:2^j
			%shearlet spectrum
            P_hor = pp.shearletSpect(xi_x, xi_y, 2^(-2*j), k*2^(-j), pp.realCoefficients, pp.shearletArg);
			if(rectangular)
                P_ver = pp.shearletSpect(xi_y, xi_x, 2^(-2*j), k*2^(-j), pp.realCoefficients, pp.shearletArg);
			else
				%the following three terms are equivalent
				%the matrix is supposed to be mirrored at the counter
				%diagonal
                %P_ver = fliplr(flipud(P_hor'));
				%P_ver = rot90(rot90(P_hor)');
				P_ver = rot90(P_hor,2)';
			end
			if(~pp.realCoefficients)
				P_ver = rot90(P_ver,2); %workaround to cover left-upper part
			end
			if( k==-2^j )
				if(pp.realCoefficients)
					Psi(:,:,start_index + idx) = P_hor .* C_hor + P_ver .* C_ver;
				else
					Psi(:,:,start_index + idx) = P_hor .* C_hor + P_ver .* C_ver;
				end
            elseif( k==2^j )
				Psi(:,:,start_index + idx + shift) = P_hor .* C_hor + P_ver .* C_ver;
            else
				new_pos = mod(idx-shift,shearsPerScale(j+1));
				if(new_pos == 0)
					new_pos = shearsPerScale(j+1);
				end
				Psi(:,:,start_index + new_pos) = P_hor;
				Psi(:,:,start_index + idx + shift) = P_ver;

				%update shift
				shift = shift + 1;
			end
		end
	end
    
    %generate output with size l
	Psi = Psi(1:l_orig(1), 1:l_orig(2), :);
	
	%modify spectra at finest scales to obtain really real shearlets
	%the modification has only to be done for dimensions with even length
	if(pp.realCoefficients && pp.realReal)
		idx_finest_scale = (1 + sum(shearsPerScale(1:end-1))) + 1;
		scale_idx = idx_finest_scale + [1:idx_finest_scale/2 idx_finest_scale/2+2:(shearsPerScale(end)-1)];
		if(lm(1)) %even number of rows -> modify first row
			idx = 2:l_orig(2);
 			Psi(1,idx,scale_idx) = 1/sqrt(2)*Psi(1,idx,scale_idx) + 1/sqrt(2)*Psi(1,l_orig(2):-1:2,scale_idx);
		end
		if(lm(2)) %even number of columns -> modify first column
			idx = 2:l_orig(1);
 			Psi(idx,1,scale_idx) = 1/sqrt(2)*Psi(idx,1,scale_idx) + 1/sqrt(2)*Psi(l_orig(1):-1:2,1,scale_idx);
		end
	end
	
end