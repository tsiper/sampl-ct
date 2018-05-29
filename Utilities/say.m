function [  ] = say( varargin )
%SAY Outputs everything the function gets with a timestamp, and also writes to log.

% 
persistent SayCounter;
if isempty(SayCounter);
    SayCounter = 1;
end

% Adding compatibility for sayhead
if strcmpi(varargin{end},'head')
    Affix = '\t';
    varargin{end} = [];
else
    Affix = '\t\t';
end
    

% Formatting the out string
outString = [...
    num2str(SayCounter),'\t',      ...
    datestr(clock,'HH:MM:SS'),Affix,...
    sprintf(varargin{:})            ...
    ];
% Promoting the serial number
SayCounter = SayCounter+1;


% Writing to screen
fprintf([outString,'\n']);

% Openining the Log file for writing
% fileID = fopen(['Log/',datestr(date),'-SuperResLog.txt'],'a');

% Writing to file
% fprintf(fileID,[outString,'\r\n']);

% Closing the file
% fclose(fileID);

end