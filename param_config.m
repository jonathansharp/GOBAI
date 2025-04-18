% process input parameter name to GOBAI functions

function param_props = param_config(p)

if strcmp(p,'o2')
    param_props.dir_name = 'O2'; % directory name
    param_props.file_name = 'o2'; % file name
    param_props.fig_name = 'oxygen'; % for saving figures
    param_props.label = '[O_{2}]'; % axis labels
    param_props.argo_name = 'DOXY'; % Argo name
    param_props.temp_name = 'OXY'; % my temporary name
    param_props.glodap_name = 'G2oxygen'; % glodap name
    param_props.woa_name = 'OXYGEN'; % woa folder
    param_props.edges = 0:5:500;
    param_props.units = '(umol kg^{-1})';
    param_props.long_param_name = 'Dissolved Oxygen Amount Content';
    param_props.dec_points = 2;
    param_props.cmap = cmocean('ice');
elseif strcmp(p,'no3')
    param_props.dir_name = 'NO3'; % directory name
    param_props.file_name = 'no3'; % file name
    param_props.fig_name = 'nitrate'; % for saving figures
    param_props.label = '[NO_{3}]'; % axis labels
    param_props.argo_name = 'NITRATE'; % Argo name
    param_props.temp_name = 'NIT'; % my temporary name
    param_props.glodap_name = 'G2nitrate'; % glodap name
    param_props.woa_name = 'NITRATE'; % woa folder
    param_props.edges = 0:0.5:50;
    param_props.units = '(umol kg^{-1})';
    param_props.long_param_name = 'Dissolved Nitrate Amount Content';
    param_props.dec_points = 2;
    param_props.cmap = flipud(cmocean('speed'));
elseif strcmp(p,'dic')
    param_props.dir_name = 'DIC'; % directory name
    param_props.file_name = 'dic'; % file name
    param_props.fig_name = 'dic'; % for saving figures
    param_props.label = '{/itC}_{T}'; % axis labels
    param_props.argo_name = 'PH_IN_SITU_TOTAL'; % Argo name
    param_props.temp_name = 'PH'; % my temporary name
    param_props.glodap_name = 'G2tcarbon'; % glodap name
    param_props.woa_name = '';
    param_props.edges = 1800:2:2400;
    param_props.units = '(umol kg^{-1})';
    param_props.long_param_name = 'Dissolved Inorganic Carbon';
    param_props.dec_points = 1;
    param_props.cmap = cmocean('matter');
end
