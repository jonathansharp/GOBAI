gobai = load_gobai('O2','v2.2','/raid/Data/GOBAI-O2/v2.2/',{'oxy' 'temp'});

idx_lat = find(gobai.lat == 20.5);
idx_lon = find(gobai.lon == -150.5);

figure; plot(1:240,squeeze(gobai.oxy(idx_lon,idx_lat,1,:)));
figure; plot(1:240,squeeze(gobai.temp(idx_lon,idx_lat,1,:)));



figure; plot(1:240,squeeze(gobai.oxy(idx_lon,idx_lat,30,:)));


figure; plot(1:240,squeeze(gobai.oxy(idx_lon,idx_lat,35,:)));

worldmap('world'); scatterm(30.5,-150.5); plot_land('map')