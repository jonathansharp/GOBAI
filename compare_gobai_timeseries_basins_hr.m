%% file information
ver = 'HR-v1.1'; % version
var = 'DIC'; % variable
var2 = 'dic';
fpath = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
name1 = ['GOBAI-' var '-HR-v202602-1993-01.nc'];
% download dimensions
GOBAI.lon = ncread([fpath 'monthly/' name1],'longitude');
GOBAI.lat = ncread([fpath 'monthly/' name1],'latitude');
GOBAI.pres = ncread([fpath 'monthly/' name1],'mean_pressure');

%% load and process basin masks
% define RECCAP2 mask file
filen = 'RECCAP2_region_masks_all_v20221025.nc';
% load and RECCAP2 dimensions
mask.lat = ncread(filen,'lat');
mask.lon = ncread(filen,'lon');
% load and interpolate RECCAP2 basin masks
basins = {'atlantic' 'pacific' 'indian' 'arctic' 'southern'};
basin_names = {'Atlantic' 'Pacific' 'Indian' 'Arctic' 'Southern'};
for b = 1:length(basins)
    mask.(basins{b}) = ncread(filen,basins{b});
    GOBAI.mask.(basins{b}) = interp2(mask.lon',mask.lat,...
        single(mask.(basins{b}))',GOBAI.lon',GOBAI.lat,'nearest')';
end

%% define coverage mask
GOBAI.mask = true(length(GOBAI.lon),length(GOBAI.lat),length(GOBAI.pres));
files = dir([fpath 'monthly/']); idx_files = false(length(files),1);
for f = 1:length(files)
    idx_files(f) = contains(files(f).name,'GOBAI');
end
files(~idx_files) = [];
for t = 1:length(files)
    inf = ncinfo([fpath 'monthly/' files(t).name]);
    for w = 1:inf.Dimensions(find(strcmp({inf.Dimensions.Name},'time'))).Length
        gobai_tmp = ncread([fpath 'monthly/' files(t).name],...
            var2,[1 1 1 w],[Inf Inf Inf 1]);
        GOBAI.mask(isnan(gobai_tmp)) = false;
    end
end
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(~GOBAI.mask) = NaN;
% save basin mask
nccreate([fpath 'GOBAI-' var '-HR-v202602-mask.nc'],'mask','Datatype','int8',...
    'Dimensions',{'lon' length(GOBAI.lon) 'lat' length(GOBAI.lat) 'pres' length(GOBAI.pres)});
ncwrite([fpath 'GOBAI-' var '-HR-v202602-mask.nc'],'mask',int8(GOBAI.mask));
% save basin volume
nccreate([fpath 'GOBAI-' var '-HR-v202602-vol.nc'],'vol','Datatype','double',...
    'Dimensions',{'lon' length(GOBAI.lon) 'lat' length(GOBAI.lat) 'pres' length(GOBAI.pres)});
ncwrite([fpath 'GOBAI-' var '-HR-v202602-vol.nc'],'vol',GOBAI.vol);

%% plot inventory for each basin
for b = 1:length(basins)
    GOBAI.time = []; cnt = 1;
    GOBAI.([basins{b} '_mean']) = [];
    for t = 1:length(files)
        time_temp = ncread([fpath 'monthly/' files(t).name],'time');
        GOBAI.time = [GOBAI.time;datenum(1950,1,1)+time_temp];
        inf = ncinfo([fpath 'monthly/' files(t).name]);
        for w = 1:inf.Dimensions(find(strcmp({inf.Dimensions.Name},'time'))).Length
            gobai_tmp = ncread([fpath 'monthly/' files(t).name],...
                var2,[1 1 1 w],[Inf Inf Inf 1]);
            gobai_tmp(~GOBAI.mask) = NaN;
            GOBAI.([basins{b} '_mean'])(cnt) = sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./...
                sum(GOBAI.vol(:),'omitnan');
            GOBAI.([basins{b} '_mean'])(cnt) = sum(gobai_tmp(:).*GOBAI.vol(:),'omitnan')./...
                sum(GOBAI.vol(:),'omitnan');
            cnt = cnt+1;
        end
    end
    plot(double(GOBAI.time),GOBAI.([basins{b} '_mean']),'LineWidth',2);

    %% figure information
    f = gcf;
    f.Position(3) = f.Position(3)*2;
    datetick('x','yyyy');
    title([basin_names{b} ' Mean GOBAI-' var ' Versions']);
    ylabel(['Weighted Average ' var]);
    %legend({'v1.0' 'v1.1-HR'},'Location','northeast');
    legend({'v1.1-HR'},'Location','southeast');

    %% export figure
    exportgraphics(gcf,['/raid/Data/GOBAI-' var '/' basins{b} '_gobai_comparison.png']);
    close
    
    %% clean up
    clear f GOBAI mask_3d mask_4d var ver fpath

end