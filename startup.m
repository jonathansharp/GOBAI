% startup.m

% add chinook matlab path
if exist('/raid/sharp/matlab','dir') == 7
    addpath('/raid/sharp/matlab');
end

% add current path
addpath(genpath(pwd));