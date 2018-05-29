p = mfilename('fullpath');
[pathstr,~,~]=fileparts(p);
addpath(genpath(pathstr));
fsepidx=strfind(pathstr,filesep);
pathstr=pathstr(1:fsepidx(end));
pathstr=[pathstr 'common'];
addpath(pathstr);