function GOBAI = load_gobai(ver,var)

path = ['/raid/Data/GOBAI-' var '/'];
GOBAI = netcdfreader([path 'GOBAI-' var '-' ver '.nc']);
% temporal index
idx_t = find(datenum(2004,1,1+double(GOBAI.time)) > datenum(2004,0,0) & ...
    datenum(2004,1,1+double(GOBAI.time)) < datenum(2018,0,0));
% climatology and amplitude
GOBAI.oxy_clim = nan(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres),12);
for m = 1:12
    GOBAI.oxy_clim(:,:,:,m) = mean(GOBAI.oxy(:,:,:,m:12:max(idx_t)),4);
end
GOBAI.oxy_amp = max(GOBAI.oxy_clim,[],4)-min(GOBAI.oxy_clim,[],4);
% annual average and IAV
GOBAI.oxy_ann = nan(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres),length(GOBAI.time(idx_t))/12);
for y = 1:length(GOBAI.time(idx_t))/12
    GOBAI.oxy_ann(:,:,:,y) = mean(GOBAI.oxy(:,:,:,(y-1)*12+1:(y-1)*12+12),4,'omitnan');
end
GOBAI.oxy_IAV = std(GOBAI.oxy_ann,[],4,'omitnan');
% trends
GOBAI.oxy_tr = nan(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres));
for a = 1:length(GOBAI.lon)
    for b = 1:length(GOBAI.lat)
        for c = 1:length(GOBAI.pres)
            if any(isnan(squeeze(GOBAI.oxy_ann(a,b,c,:))))
                GOBAI.oxy_tr(a,b,c) = NaN;
            else
                [~,~,x] = leastsq2(2004:2017,squeeze(GOBAI.oxy_ann(a,b,c,:)),2003,0,0);
                GOBAI.oxy_tr(a,b,c) = x(2);
            end
        end
    end
end
% means
GOBAI.temp = mean(GOBAI.temp(:,:,:,idx_t),4,'omitnan');
GOBAI.sal = mean(GOBAI.sal(:,:,:,idx_t),4,'omitnan');
GOBAI.oxy = mean(GOBAI.oxy(:,:,:,idx_t),4,'omitnan');
GOBAI.uncer = mean(GOBAI.uncer(:,:,:,idx_t),4,'omitnan');
% clean up
GOBAI = rmfield(GOBAI,{'oxy_clim' 'oxy_ann'});