% replicate_dims
%
% DESCRIPTION:
% This function is used to replicate the 
% spatiotemporal dimensions of gridded products.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 11/1/2024

function TS = replicate_dims(base_grid,TS,tdim)

% replicate dmensions
TS.lon_cos_1 = repmat(cosd(TS.Longitude-20),1,TS.ydim,TS.zdim,tdim);
TS.lon_cos_2 = repmat(cosd(TS.Longitude-110),1,TS.ydim,TS.zdim,tdim);
TS.longitude = repmat(TS.Longitude,1,TS.ydim,TS.zdim,tdim);
TS.latitude = repmat(TS.Latitude',TS.xdim,1,TS.zdim,tdim);
if strcmp(base_grid,'RG') || strcmp(base_grid,'RFROM')
    TS.pressure = repmat(permute(TS.Pressure,[3 2 1]),TS.xdim,TS.ydim,1,tdim);
else
    TS.depth = repmat(permute(TS.Depth,[3 2 1]),TS.xdim,TS.ydim,1,tdim);
    TS.pressure = -gsw_p_from_z(TS.depth,TS.latitude);
end
if strcmp(base_grid,'RG')
    TS.time = repmat(permute(TS.Time,[4 3 2 1]),TS.xdim,TS.ydim,TS.zdim,1);
end

