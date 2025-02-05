% load gobai
gobai = load_gobai('O2','ann-v2.2','/raid/Data/GOBAI-O2/v2.2/',{'oxy' 'temp'});

% calculate weights
[vol,area,h] = weights3d(gobai.lon,gobai.lat,gobai.pres);

% extract 2023 and baseline oxygen
gobai.oxy_2023 = gobai.oxy(:,:,:,end);
gobai.oxy_base = mean(gobai.oxy(:,:,:,1:end-1),4,'omitnan');

% define index
idx = isnan(mean(gobai.oxy_2023,4,'omitnan'));
h(idx) = NaN;

% calculate weighted means 
gobai.oxy_2023_0_200 = sum(gobai.oxy_2023.*h,3,'omitnan')./sum(h,3,'omitnan');
gobai.oxy_base_0_200 = sum(gobai.oxy_base.*h,3,'omitnan')./sum(h,3,'omitnan');


figure;
worldmap('world');
pcolorm(double(gobai.lat),double(gobai.lon),double(gobai.oxy_2023_0_200-gobai.oxy_base_0_200)');
colorbar;
caxis([-5 5]);
colormap(cmocean('balance','pivot',0));