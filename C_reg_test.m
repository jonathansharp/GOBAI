fpath = '/raid/Data/GOBAI-O2/v2.3/';
fname = 'GOBAI-O2-v2.3.nc';

lon = ncread([fpath fname],'lon');
lat = ncread([fpath fname],'lat');
pres = ncread([fpath fname],'pres');
time = ncread([fpath fname],'time');

for n = 1:1000
    for t = 1:length(time)
        o2 = ncread([fpath fname],'oxy',[1 1 1 t],[Inf Inf Inf 1]);
        temp = ncread([fpath fname],'temp',[1 1 1 t],[Inf Inf Inf 1]);
        sal = ncread([fpath fname],'sal',[1 1 1 t],[Inf Inf Inf 1]);
        uncer = ncread([fpath fname],'uncer',[1 1 1 t],[Inf Inf Inf 1]);
        o2_sat = o2satv2b(sal,temp);
        aou = o2_sat - o2;
        creg = aou.*(117/160);
        keyboard
    end
    
end