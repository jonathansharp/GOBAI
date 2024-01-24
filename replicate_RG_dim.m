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

function TS = replicate_RG_dim(TS)

% replicate dmensions
TS.lon_cos_1_3D = repmat(cosd(TS.Longitude-20),1,TS.ydim,TS.zdim,12);
TS.lon_cos_2_3D = repmat(cosd(TS.Longitude-110),1,TS.ydim,TS.zdim,12);
TS.latitude_3D = repmat(TS.Latitude',TS.xdim,1,TS.zdim,12);
TS.pressure_3D = repmat(permute(TS.Pressure,[3 2 1]),TS.xdim,TS.ydim,1,12);
TS.time_3D = repmat(permute(TS.Time,[4 3 2 1]),TS.xdim,TS.ydim,TS.zdim,1);
