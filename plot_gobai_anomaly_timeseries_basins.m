%% file information
ver = 'v2.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path

%% load and process basin masks
% define RECCAP2 mask file
filen = 'RECCAP2_region_masks_all_v20221025.nc';
% load and process RECCAP2 mask latitude
mask.lat = ncread(filen,'lat');
mask.lat = mask.lat(26:170);
% load and process RECCAP2 mask longitude
mask.lon = ncread(filen,'lon');
mask.lon = convert_lon(mask.lon);
mask.lon = [mask.lon(21:end);mask.lon(1:20)];
% load RECCAP2 basin masks
basins = {'atlantic' 'pacific' 'indian' 'arctic' 'southern'};
basin_names = {'Atlantic' 'Pacific' 'Indian' 'Arctic' 'Southern'};
% reorder for 20-380
for b = 1:length(basins)
    mask.(basins{b}) = ncread(filen,basins{b});
    mask.(basins{b}) = [mask.(basins{b})(21:end,26:170);mask.(basins{b})(1:20,26:170)];
end
clear filen

% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
if strcmp(ver,'v2.1')
GOBAI.time = GOBAI.time + datenum(2004,0,0);
else
GOBAI.time = GOBAI.time + datenum(1950,0,0);
end

for b = 1:length(basins)

    %% figure information
    f = gcf;
    f.Position(3) = f.Position(3)*2;
    f.Position(4) = f.Position(4)*2;
    tiledlayout(f,2,1);

    % download oxygen
    GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
    GOBAI.uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');
    % calculate weights
    GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
    GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;

    % process basin
    mask_3d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres));
    mask_4d = repmat(mask.(basins{b})>0,1,1,length(GOBAI.pres),length(GOBAI.time));
    GOBAI.oxy(~mask_4d) = NaN;
    GOBAI.vol(~mask_3d) = NaN;
    
    % calculate global mean
    global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
        length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
        GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
    [yf,yr,x] = leastsq2(double(GOBAI.time),global_mean,double(GOBAI.time(1)),...
        3,[365.2425/2 365.2425 365.2425*12]);

    %% plot timeseries with fit
    nexttile;
    plot(double(GOBAI.time),yf,double(GOBAI.time),global_mean,'LineWidth',2);
    datetick('x','yyyy');
    title([basin_names{b} ' Mean GOBAI-O_{2}']);
    ylabel('Weighted Average [O_{2}]');
    legend({'Least-squares fit' 'Weighted Global Mean [O_{2}]'})
    
    %% plot anomaly
    nexttile;
    plot(double(GOBAI.time),yr,'LineWidth',2);
    datetick('x','yyyy');
    title([basin_names{b} ' Mean GOBAI-O_{2} Anomaly']);
    ylabel('Weighted Average [O_{2}] Anomaly');
    legend({'Residuals from fit'})
    
    %% export figure
    exportgraphics(gcf,['/raid/Data/GOBAI-O2/global_gobai_anomaly_' ...
        basin_names{b} '_' ver '.png']);
    close

end

clear