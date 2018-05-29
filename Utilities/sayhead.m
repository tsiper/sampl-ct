function [  ] = sayhead( varargin )
%SAY Output a headline

outStr = sprintf(varargin{:});

say(outStr,'head');
say(repmat('=',1,length(outStr)),'head');

end
