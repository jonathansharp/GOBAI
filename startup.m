% startup.m

% add chinook utilities path
if exist('/raid/sharp/matlab','dir') == 7
    addpath(genpath('/raid/sharp/matlab/utilities/'));
end

% add current path
addpath(genpath(pwd));

% suppress OpenGL warning
warning('off','MATLAB:hg:AutoSoftwareOpenGL');
warning('off','MATLAB:prnRenderer:opengl');
warning('off','MATLAB:opengl:deprecatedRendererQuery');
warning('off','stats:gmdistribution:cluster:MissingData');
