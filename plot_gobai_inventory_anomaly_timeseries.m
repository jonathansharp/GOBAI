%% file information
ver = 'v2.3'; % version
var = 'O2'; % variable
path = ['/raid/Data/GOBAI-' var '/' ver '/']; % file path

%% figure information
f = gcf;
f.Position(3) = f.Position(3)*2;
f.Position(4) = f.Position(4)*2;
tiledlayout(f,2,1);

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
% download oxygen
GOBAI.oxy = ncread([path 'GOBAI-' var '-' ver '.nc'],'oxy');
GOBAI.uncer = ncread([path 'GOBAI-' var '-' ver '.nc'],'uncer');
% calculate weights
GOBAI.vol = weights3d(GOBAI.lon,GOBAI.lat,GOBAI.pres);
GOBAI.vol(isnan(mean(GOBAI.oxy,4,'omitnan'))) = NaN;
GOBAI.dens = 1025;
GOBAI.mass = GOBAI.vol .* GOBAI.dens;
% calculate global inventroy (umol/kg * kg/m3 * m3)
global_inv = sum(reshape(GOBAI.oxy,[length(GOBAI.lon)*...
    length(GOBAI.lat)*length(GOBAI.pres) length(GOBAI.time)]).*...
    GOBAI.mass(:),'omitnan')./(10^18);
[yf,yr,x] = leastsq2(double(GOBAI.time),global_inv,double(GOBAI.time(1)),...
    3,[365.2425/2 365.2425 365.2425*12]);

%% plot timeseries with fit
nexttile;
plot(double(GOBAI.time),yf./1000,double(GOBAI.time),global_inv./1000,'LineWidth',2);
datetick('x','yyyy');
%ylim([149 153]);
title('Global GOBAI-O_{2} Inventory');
ylabel('O_{2} Inventory (Pmol)');
legend({'Least-squares fit' 'Weighted Global Mean [O_{2}]'})

%% obtain enso index
enso_table = readtable('meiv2.data.txt');
enso_index = table2array(enso_table(26:46,2:13))';
enso_index = enso_index(:);
enso_time = datenum(repelem(2004:2024,12)',repmat(1:12,1,21)',15);
enso_index_m = movmean(enso_index,12);

%% plot anomaly
nexttile; hold on; box on;
title('Global GOBAI-O_{2} Inventory Anomaly');
yyaxis left;
plot(double(GOBAI.time),yr,'LineWidth',2);
datetick('x','yyyy');
%ylim([-0.3 0.3]);
ylabel('O_{2} Inventory Anomaly (Tmol)');
yyaxis right;
corr_coeff = corr(enso_index_m,yr);
plot(enso_time,enso_index_m,'LineWidth',2);
ylabel('Multivariate ENSO Index');
text(datenum(2020,1,1),1.5,['\itr\rm = ' num2str(corr_coeff)])

%% export figure
exportgraphics(gcf,['/raid/Data/GOBAI-O2/global_gobai_inventory_anomaly_' ver '.png']);
close
clear