function plot_anom_pres(lat,lon,var,min_depth,max_depth,GOBAI)

% get variable label
label = getlabel(var);

% get colorbar limits
color_lims = getcolorlims(var);

% convert input lon to 0:360
lon(lon<0) = lon(lon<0) + 360;

% match input lat
if length(lat) > 1
    lat_idx = find(GOBAI.lat(1,:,1,1)' > min(lat) & GOBAI.lat(1,:,1,1)' < max(lat));
else
    lat_idx = find(abs(GOBAI.lat(1,:,1,1)'-lat) == min(abs(GOBAI.lat(1,:,1,1)'-lat)));
end

% match input lon
if length(lon) > 1
    lon_idx = find(GOBAI.lon(:,1,1,1) > min(lon) & GOBAI.lon(:,1,1,1) < max(lon));
else
    lon_idx = find(abs(GOBAI.lon(:,1,1,1)-lon) == min(abs(GOBAI.lon(:,1,1,1)-lon)));
end

% determine weights by lat and lon
area_weights = ...
(((GOBAI.lat + 0.25) - ...
    (GOBAI.lat - 0.25)) .* 110.574) .* ... % lat distance
(((GOBAI.lon + 0.25) - ...
    (GOBAI.lon - 0.25)) .* ...
    111.320.*cosd(GOBAI.lat(:,:,1,1))); % lon distance
area_weights = area_weights(lon_idx,lat_idx);
area_weights = repmat(permute(area_weights,[4 3 1 2]),...
    size(GOBAI.time,4),size(GOBAI.pressure,3),1,1);

% extract variables to plot
time = squeeze(GOBAI.time(1,1,1,:)); % time array
pressure = squeeze(GOBAI.pressure(1,1,:,1)); % pressure array
variable = GOBAI.([var '_resid_notrend']); % variable to plot
variable = cell2mat(permute(variable(lon_idx,lat_idx,:),[4 3 1 2]));

% define area weights over land as NaN
area_weights(isnan(variable)) = NaN;

% calculate weighted means
variable = variable.*(area_weights./sum(sum(area_weights,3,'omitnan'),4,'omitnan'));
variable = sum(sum(variable,3,'omitnan'),4,'omitnan')'; % weighted mean variable to plot
variable(variable==0) = NaN;
sigma = permute(GOBAI.sigma0(lon_idx,lat_idx,:,:),[4 3 1 2]);
sigma = sigma.*(area_weights./sum(sum(area_weights,3,'omitnan'),4,'omitnan'));
sigma = sum(sum(sigma,3,'omitnan'),4,'omitnan')'; % weighted mean sigma
sigma(sigma==0) = NaN;

% create contour plot
figure; hold on;
set(gcf,'units','normalized','position',[0 0.5 1 0.5]);
set(gca,'ydir','reverse','fontsize',28);
idx_pres = ~isnan(pressure);
contourf(time,pressure(idx_pres),variable(idx_pres,:),'linestyle','none');
contour(time,pressure(idx_pres),sigma(idx_pres,:),[20:26 26.5:0.25:28],...
    'k','linewidth',3,'linestyle',':','ShowText','on');
caxis(color_lims);
c=colorbar;
colormap(cmocean('balance','pivot',0));
ylim([min_depth max_depth]);
ylabel('Depth (dbar)');
datetick('x','yyyy','keeplimits');
c.Label.String = label;
%c.Label.FontSize = fntsz;

function label = getlabel(var)

if contains(var,'temp')
    label = '\DeltaT (degC)';
elseif contains(var,'sal')
    label = '\DeltaS';
elseif contains(var,'oxy') && contains(var,'sat')
    label = '\Delta%O_{2(sat)}';
elseif contains(var,'oxy')
    label = '\DeltaO_{2} (\mumol kg^{-1})';
elseif contains(var,'no3')
    label = '\Delta[NO_{3}] (\mumol kg^{-1})';
elseif contains(var,'spice')
    label = 'Spice';
elseif contains(var,'sigma0')
    label = 'Sigma';
end

function color_lims = getcolorlims(var)

if contains(var,'temp')
    color_lims = [-2 2];
elseif contains(var,'sal')
    color_lims = [-1 1];
elseif contains(var,'oxy') && contains(var,'sat')
    color_lims = [-10 10];
elseif contains(var,'oxy')
    color_lims = [-20 20];
elseif contains(var,'no3')
    color_lims = [-50 50];
elseif contains(var,'spice')
    color_lims = [-1 1];
elseif contains(var,'sigma0')
    color_lims = [-1 1];
else
    color_lims = [-100 100];
end