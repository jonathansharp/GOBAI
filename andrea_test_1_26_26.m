gobai = '/raid/Data/GOBAI-O2/v2.3/GOBAI-O2-v2.3.nc';
gobai_dic = '/raid/Data/GOBAI-DIC/v1.0/GOBAI-DIC-v1.0.nc';
time  = double(ncread(gobai,'time')) + datenum('1950-01-01');
gpres  = ncread(gobai,'pres');
glat   = ncread(gobai,'lat');
glon   = ncread(gobai,'lon');
xlon   = find(glon == 215.5);
xlat   = find(glat == 50.5)
xp     = find(gpres == 440);
dic    = squeeze(double(ncread(gobai_dic,'dic',[xlon xlat xp 1],[1 1 1 Inf])));
plot(time,dic)