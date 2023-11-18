function depth = bottom_depth(latitude,longitude)

% uses ETOPO2 bathymetry file to calculate bottom depth from
% input latitude and longitude

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
    
load('ETOPO2/ETOPO2.mat');

depth_vec = ...
    interp2(ETOPO2.lat,ETOPO2.lon,ETOPO2.bottomdepth,...
    double(lat_vec),double(lon_vec),'nearest');

depth = reshape(depth_vec,size(lat_mat));

end