% calculate volume from latitude, longitude, and depth
function vol = calculate_volume(lat,lon,depth,depth_bnds)

rE = 6.371e6; % radius of Earth (m)
dlat = lat(3) - lat(2); % spacing between latitudes
area_3d = nan(length(lon),length(lat),length(depth));
depth_3d = repmat(permute(diff(depth_bnds),[3 2 1]),length(lon),length(lat),1);
for i = 1:length(lat)
    lat_area = 2*pi*rE^2*abs(sind(lat(i)-dlat/2)-sind(lat(i)+dlat/2));
    area_3d(:,i,:) = lat_area/length(lon); % m^2
end

vol = area_3d.*depth_3d; % m^3