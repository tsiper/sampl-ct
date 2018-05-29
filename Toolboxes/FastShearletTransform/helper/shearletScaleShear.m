function varargout = shearletScaleShear(a,b,c,d)
% compute index from scale j, shear k and cone and vice versa. Optionally
% return values of shearlets (or coefficients) for given index or given
% scale j, shear k and cone.
%
% OPTIONS
%% scale, shear and cone from index
% [j,k,cone] = shearletScaleShear(a)
% INPUT:
%  a	(int) index
%
% OUTPUT:
%  j		(int) scale j (>= 0)
%  k		(int) shear k, -2^j <= k <= 2^j
%  cone		(char) cone {h,v,x,0}
%% return data for index
% ST = shearletScaleShear(a,b)
% INPUT:
%  a    (3-d-matrix) shearlets or shearlet coefficients
%  b	(int) index
%
% OUTPUT:
%  ST		(matrix) layer of ST for index [ST(:,:,index)]
%% index from scale, shear and cone
% index = shearletScaleShear(a,b,c)
% INPUT:
%  a		(int) scale j (>= 0)
%  b		(int) shear k, -2^j <= k <= 2^j
%  c		(char) cone {h,v,x,0}
%
% OUTPUT:
%  index	(int) respective index
%% return data for j,k and cone
% ST = shearletScaleShear(a,b,c,d)
% INPUT:
%  a		(3-d-matrix) shearlets or shearlet coefficients
%  b		(int) scale j (>= 0)
%  c		(int) shear k, -2^j <= k <= 2^j
%  d		(char) cone {h,v,x,0}
%
% OUTPUT:
%  index	(int) respective index
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

%% display informations
disp = 0;

%% different cases
if(nargin == 1)
	%compute j and k from index
	index = a;
 	[j,k,cone] = index2jk(index);
 	varargout = {j,k,cone};
	if( disp == 1)
		fprintf('index %d represents:\n',index)
		fprintf('scale j: %d (a = %.4f)\n',j,4^(-j));
		fprintf('shear k: %d (s = %.4f)\n',k,2^(-j)*k);
		fprintf('cone   : %s\n',cone);
	end
elseif(nargin == 2)
	%return data for index
	ST = a;
	index = b;
	varargout = {ST(:,:,index)};
	[j,k,cone] = index2jk(index);
	
	if( disp == 1)
		fprintf('index %d represents:\n',index)
		fprintf('scale j: %d (a = %.4f)\n',j,4^(-j));
		fprintf('shear k: %d (s = %.4f)\n',k,2^(-j)*k);
		fprintf('cone   : %s\n',cone);
	end
elseif(nargin == 3) 
		%compute index from j and k and cone
		j = a;
		k = b;
		cone = c;
		index = jk2index(j,k,cone);
		varargout = {index};
		
		if( disp == 1)
			fprintf('index %d represents:\n',index)
			fprintf('scale j: %d (a = %.4f)\n',j,4^(-j));
			fprintf('shear k: %d (s = %.4f)\n',k,2^(-j)*k);
			fprintf('cone   : %s\n',cone);
		end
elseif(nargin == 4)
	%return data for j and k and cone
	ST = a;
	j = b;
	k = c;
	cone = d;
	
	index = jk2index(j,k,cone);
	varargout = {ST(:,:,index)};
	
	if( disp == 1)
		fprintf('index %d represents:\n',index)
		fprintf('scale j: %d (a = %.4f)\n',j,4^(-j));
		fprintf('shear k: %d (s = %.4f)\n',k,2^(-j)*k);
		fprintf('cone   : %s\n',cone);
	end
end

end


function index = jk2index(j,k,cone)
% helper function, compute index from j, k and cone
	%lowpass
	index = 1;
	if(isnan(j) && isnan(k) && strcmp(cone,'0'))
		return;
	else
		% sum of lower scales
		index = index + sum(2.^(2+(0:j-1)));
				
		%get detail index from shear (and cone!)
		switch cone
			case 'h'
				if( k <= 0 )
					index = index - k;
				else
					index = index + 4*2^j - k;
				end
			case 'v'
				index = index + 2^j + (k + 2^j);
			case 'x'
				index = index + (2 + sign(k)) * 2^j;
		end
			
		% sligth adjustment ( k=0 <=> index = 1)
		index = index + 1;
	end
end

function [j,k,cone] = index2jk(index)
% helper function, compute j, k and cone from index
	if(index <= 1) %lowpass, j and k not needed
		j = NaN;
		k = NaN;
		cone = '0';
	else
		%substract 1 for the lowpass
		index = index - 1;
		
		%determine scale j
		%substract number of shears in each scale:
		% 2^(j+0), 2^(j+1), 2^(j+2), ...
		j = 0;
		while ( index > 2^(2 + j) )
			index = index - 2^(j+2);
			j = j + 1;
		end
		
		% shift to zero (first index <=> k=0)
		index = index - 1;
		
		% determine cone
		% index | 0 1 ... 2^j ... 2*2^j ... 3*2^j ... 4*2^j -1 
		% k     | 0 -1   -2^j       0        2^j        1
		% cone  | h ... h  x  v ... v ... v   x h ...   h
		index2 = index / 2^j;
		
		if ( index2 < 1)
			k = -index;
			cone = 'h';
		elseif ( index2 == 1 )
			k = -2^j;
			cone = 'x';
		elseif ( index2 < 3 )
			k = index -2*2^j;
			cone = 'v';
		elseif (index2 == 3)
			k = 2^j;
			cone = 'x';
		else
			k = -(index - 4*2^j);
			cone = 'h';
		end
	end
end



	
	