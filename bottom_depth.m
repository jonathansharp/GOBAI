function depth = bottom_depth(latitude,longitude)

% uses ETOPOv2022 bathymetry file to calculate bottom depth from
% input latitude and longitude

% convert longitude if needed
longitude(longitude>180) = longitude(longitude>180)-360;

if size(latitude) == size(longitude)
    lat_mat = latitude;
    lon_mat = longitude;
    lat_vec = lat_mat(:);
    lon_vec = lon_mat(:);
    depth_vec = nan(size(lat_vec));
else
    [lat_mat,lon_mat] = meshgrid(latitude,longitude);
    lat_vec = lat_mat(:);
    lon_vec = lon_mat(:);
    depth_vec = nan(size(lat_vec));
end

% load ETOPOv2022
ETOPO.lon = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lon');
ETOPO.lat = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','lat');
ETOPO.z = ncread('ETOPO_2022_v1_60s_N90W180_bed.nc','z');

depth_vec = ...
    interp2(ETOPO.lat,ETOPO.lon,ETOPO.z,...
    double(lat_vec),double(lon_vec),'nearest');

depth = reshape(depth_vec,size(lat_mat));

end