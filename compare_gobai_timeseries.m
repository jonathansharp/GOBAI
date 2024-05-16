%% v1.0
% file information
ver = 'v1.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
figure; plot(double(datenum(0,0,0)+GOBAI.time),global_mean,'LineWidth',2); hold on

%% v2.0
% file information
ver = 'v2.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(0,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% v2.1
% file information
ver = 'v2.1'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(2004,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% v 2.2
% file information
ver = 'v2.2'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
% calculate global mean
global_mean = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(sum(GOBAI.vol(:),'omitnan'));
plot(double(datenum(1950,0,0)+GOBAI.time),global_mean,'LineWidth',2);

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
datetick('x','yyyy');
title('Global Mean GOBAI-O_{2} Versions');
ylabel('Weighted Average [O_{2}]');
legend({'v1.0' 'v2.0' 'v2.1' 'v2.2'},'Location','northeast');

%% export figure
exportgraphics(gcf,'global_gobai_comparison.png');
close
clear


%% v2.1 AOU
% file information
ver = 'v2.0'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path
% download dimensions
GOBAI.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
GOBAI.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
GOBAI.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
GOBAI.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time');
% download oxygen, temperature, salinity
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
GOBAI.temp = ncread([path 'GOBAI-' var '-' ver '.nc'],'temp');
GOBAI.sal = ncread([path 'GOBAI-' var '-' ver '.nc'],'sal');
% save (or load) NetCDFs
fname = ['/raid/Data/GOBAI-' var '/' ver '/GOBAI-' var '-' ver '_anc.nc'];
vars = {'sa','ct','dens','aou'};
if ~exist(fname,'file')
    % calculate density, aou
    pres3d = repmat(permute(GOBAI.pres,[3 2 1]),length(GOBAI.lon),length(GOBAI.lat),1);
    lon3d = repmat(GOBAI.lon,1,length(GOBAI.lat),length(GOBAI.pres));
    lat3d = repmat(GOBAI.lat',length(GOBAI.lon),1,length(GOBAI.pres));
    GOBAI.sa = single(nan(size(GOBAI.oxy)));
    GOBAI.ct = single(nan(size(GOBAI.oxy)));
    GOBAI.dens = single(nan(size(GOBAI.oxy)));
    GOBAI.aou = single(nan(size(GOBAI.oxy)));
    for t = 1:length(GOBAI.time)
        GOBAI.sa(:,:,:,t) = gsw_SA_from_SP(GOBAI.sal(:,:,:,t),pres3d,lon3d,lat3d);
        GOBAI.ct(:,:,:,t) = gsw_CT_from_t(GOBAI.sa(:,:,:,t),GOBAI.temp(:,:,:,t),pres3d);
        GOBAI.dens(:,:,:,t) = gsw_rho(GOBAI.sa(:,:,:,t),GOBAI.ct(:,:,:,t),pres3d);
        GOBAI.aou(:,:,:,t) = gsw_O2sol(GOBAI.sa(:,:,:,t),GOBAI.ct(:,:,:,t),pres3d,lon3d,lat3d)-GOBAI.oxy(:,:,:,t);
    end
    % save file
    for v = 1:length(vars)
        nccreate(fname,vars{v},'Dimensions',{'lon',length(GOBAI.lon),'lat',length(GOBAI.lat),...
            'pres',length(GOBAI.pres),'time',length(GOBAI.time)},'Datatype','single');
        ncwrite(fname,vars{v},GOBAI.(vars{v}));
    end
else
    % load variables
    for v = 1:length(vars)
        GOBAI.(vars{v}) = ncread(fname,vars{v});
    end
end
% calculate weights
[GOBAI.vol,GOBAI.area,GOBAI.heights] = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
GOBAI.area(isnan(mean(GOBAI.area,4,'omitnan'))) = NaN;
% restrict to 200-1000 dbar
GOBAI.pres = GOBAI.pres(20:44);
GOBAI.oxy = GOBAI.oxy(:,:,20:44,:);
GOBAI.aou = GOBAI.aou(:,:,20:44,:);
GOBAI.dens = GOBAI.dens(:,:,20:44,:);
GOBAI.vol = GOBAI.vol(:,:,20:44,:);
GOBAI.heights = GOBAI.heights(:,:,20:44,:);
% calculate global mean
global_inv_oxy = sum(reshape(GOBAI.oxy.*GOBAI.dens,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(1e21); % umol/kg * kg/m^3 * m^3 * 1 Pmol/10^21 umol
global_inv_aou = sum(reshape(GOBAI.aou.*GOBAI.dens,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.vol(:),'omitnan')./(1e21);
global_inv_creg = global_inv_aou.*(117/170);
% calculate global mean using trapezoidal integration
GOBAI.oxy_col_int = squeeze(trapz(GOBAI.pres,GOBAI.oxy.*GOBAI.dens,3));
GOBAI.aou_col_int = squeeze(trapz(GOBAI.pres,GOBAI.aou.*GOBAI.dens,3));
area_weights = GOBAI.area(:,:,1);
global_inv_oxy_2 = sum(reshape(GOBAI.oxy_col_int,[length(GOBAI.lon)*...
    length(GOBAI.lat) length(GOBAI.time)]).*...
    area_weights(:),'omitnan')./(1e21); % umol/kg * kg/m^3 * m^3 * 1 Pmol/10^21 umol
global_inv_aou_2 = sum(reshape(GOBAI.aou_col_int,[length(GOBAI.lon)*...
    length(GOBAI.lat) length(GOBAI.time)]).*...
    area_weights(:),'omitnan')./(1e21); % umol/kg * kg/m^3 * m^3 * 1 Pmol/10^21 umol
global_inv_creg_2 = global_inv_aou_2.*(117/170);

% calculate column inventories (mol/m2)
global_column_inv_aou = ...
    sum(mean(GOBAI.aou,4).*mean(GOBAI.dens,4).*GOBAI.heights,3)./(10^6);

% plot AOU column inventory
figure; worldmap([-90 90],[20 380]);
title('AOU in mol m^{-2} (200 - 1000 dbar)')
pcolorm(double(GOBAI.lat),double(GOBAI.lon),double(global_column_inv_aou)');
colorbar; plot_land('map');
mlabel off; plabel off;
exportgraphics(gcf,['Figures/aou_col_inv_' ver '.png']);
close

% aou and o2
figure;
yyaxis left; 
plot(double(datenum(1950,0,0)+GOBAI.time),global_inv_oxy_2,'LineWidth',2); hold on
ylabel('[O_{2}] from 200 to 1000 dbar (Pmol)');
datetick('x');
yyaxis right;
plot(double(datenum(1950,0,0)+GOBAI.time),-global_inv_aou_2,'LineWidth',2);
ylabel('[-AOU from 200 to 1000 dbar (Pmol)');
exportgraphics(gcf,['Figures/oxy_aou_inv_' ver '.png']);
close

% C_reg
figure;
plot(double(datenum(2004,0,0)+GOBAI.time),global_inv_creg_2.*12.011,'LineWidth',2);
ylabel('C_{reg.} from 200 to 1000 dbar (Pg)');
exportgraphics(gcf,['Figures/creg_inv_' ver '.png']);
close

% trend fit
num_pers = 2;
[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(double(datenum(2004,0,0)+GOBAI.time),global_inv_oxy_2,...
        double(datenum(2004,0,0)+GOBAI.time(1)),num_pers,[365.25/2 365.25]);
% Determined lagged autocorrelation of residuals
[acov,acor,lag,dof] = autocov((1:228)',yr,228);
% get decadal trend in Pmol
dec_O2 = x(2)*365.25*10*10^3;
% Scale error estimated from covariance matrix by effective degrees of
% freedom, 12 months, and 95% confidence
eDOF_scale = 2+num_pers*2;
edof = dof-eDOF_scale;
if edof < 1; edof = 1; end
err95_inv = err(2)*(sqrt(228)/sqrt(edof))*365.25*10*10^3*2;
disp(['O_{2} Inv.: ' num2str(dec_O2) ' +/- ' num2str(err95_inv) ' Tmol/dec']);

[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(double(datenum(2004,0,0)+GOBAI.time),global_inv_aou_2,...
        double(datenum(2004,0,0)+GOBAI.time(1)),2,[365/2 365]);
dec_AOU = x(2)*365*10*10^3;
disp(['AOU Inv.: ' num2str(dec_AOU) ' Tmol/dec']);
disp(['AOU Per.: ' num2str(100*(-dec_AOU/dec_O2)) ' %']);

[yf,yr,x,err,corrmat,r2,n2] = ...
    leastsq2(double(datenum(2004,0,0)+GOBAI.time),global_inv_creg_2,...
        double(datenum(2004,0,0)+GOBAI.time(1)),2,[365/2 365]);
dec_Creg = x(2)*365*10*10^3;
disp(['C_reg Inv.: ' num2str(dec_Creg) ' Tmol/dec']);

writematrix([GOBAI.time,global_inv_oxy',...
    global_inv_aou',global_inv_creg'],'gobai_timeseries.csv')
