% calculate volume from latitude, longitude, and depth
function area = calculate_area(lat,lon)

rE = 6.371e6; % radius of Earth (m)
dlat = lat(3) - lat(2); % spacing between latitudes
area = nan(length(lon),length(lat));
for i = 1:length(lat)
    lat_area = 2*pi*rE^2*abs(sind(lat(i)-dlat/2)-sind(lat(i)+dlat/2));
    area(:,i) = lat_area/length(lon); % m^2
end