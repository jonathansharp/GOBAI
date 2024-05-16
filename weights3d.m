function [vol,area,h] = weights3d(lon,lat,depth)

% replicate dimensions
lon3d = repmat(lon,1,length(lat),length(depth));
lat3d = repmat(lat',length(lon),1,length(depth));
% volume = length x width x depth
L = (((lon3d + 0.5) - (lon3d- 0.5)) .* 111.320 .* cosd(lat3d)); % longitude distance (km)
W = (((lat3d + 0.5) - (lat3d - 0.5)) .* 110.574); % latitude distance (km)
% detemine depth bounds and heights of grid cells
D = single(nan(size(depth)));
d1 = 0;
for z = 1:length(depth)
    d0 = d1;
    f = depth(z) - d0;
    d1 = depth(z) + f;
    D(z) = d1 - d0;
end
clear d0 d1 f 
D = repmat(permute(D,[3 2 1]),length(lon),length(lat),1); % meters
% calculate volumes and areas
vol = (L*1e3).*(W*1e3).*D; % meters cubed
area = (L*1e3).*(W*1e3); % meters squared
h = D; % meters
