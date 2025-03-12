% convert longitude from either 0:360 to -180:180 or vice-versa
% can input 'format','0:360' or 'format','-180:180'

function lon = convert_lon(lon,varargin)

frmt = 'none'; % set default
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'format')
        frmt = varargin{i+1};
    end
end

if strcmp(frmt,'none')
    if any(lon > 180)
        % convert from 0:360 to -180:180
        lon(lon>180) = lon(lon>180)-360;
        lon(lon<-180) = lon(lon<-180)+360;
    elseif any(lon < 0)
        % convert from -180:180 to 0:360
        lon(lon<0) = lon(lon<0)+360;
        lon(lon>360) = lon(lon>360)-360;
    else
        % do nothing
    end
elseif strcmp(frmt,'0:360') || strcmp(frmt,'0-360') || strcmp(frmt,'0 to 360')
    % convert from -180:180 to 0:360
    lon(lon<0) = lon(lon<0)+360;
    lon(lon>360) = lon(lon>360)-360;
elseif strcmp(frmt,'-180:180') || strcmp(frmt,'-180-180') || strcmp(frmt,'-180 to 180')
    % convert from 0:360 to -180:180
    lon(lon>180) = lon(lon>180)-360;
    lon(lon<-180) = lon(lon<-180)+360;
end
