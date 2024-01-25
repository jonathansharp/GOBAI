% replicate_RG_dim
%
% DESCRIPTION:
% This function is used to replicate the 
% spatiotemporal dimensions of the Roemmich 
% and Gilson temperature and salinity data product.
%
% AUTHOR: J. Sharp, UW CICOES / NOAA PMEL
%
% DATE: 1/24/2024

function TS = replicate_RG_dim(TS,tdim)

% replicate dmensions
TS.lon_cos_1 = repmat(cosd(TS.Longitude-20),1,TS.ydim,TS.zdim,tdim);
TS.lon_cos_2 = repmat(cosd(TS.Longitude-110),1,TS.ydim,TS.zdim,tdim);
TS.longitude = repmat(TS.Longitude,1,TS.ydim,TS.zdim,tdim);
TS.latitude = repmat(TS.Latitude',TS.xdim,1,TS.zdim,tdim);
TS.pressure = repmat(permute(TS.Pressure,[3 2 1]),TS.xdim,TS.ydim,1,tdim);
if tdim > 1
    TS.time = repmat(permute(TS.Time,[4 3 2 1]),TS.xdim,TS.ydim,TS.zdim,1);
end
