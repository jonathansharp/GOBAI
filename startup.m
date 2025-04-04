% startup.m

% add chinook utilities path
if exist('/raid/sharp/matlab','dir') == 7
    addpath(genpath('/raid/sharp/matlab/utilities/'));
end

% add current path
addpath(genpath(pwd));

% suppress OpenGL warning
warning('off','MATLAB:hg:AutoSoftwareOpenGL');
