% process input parameter name to GOBAI functions

function param_props = param_config(p)

if strcmp(p,'o2')
    param_props.dir_name = 'O2'; % directory name
    param_props.file_name = 'o2'; % file name
    % param_props.p2 = 'oxygen'; % for saving figures?
    param_props.label = '[O_{2}]'; % axis labels
    param_props.argo_name = 'DOXY'; % Argo name
    % param_props.p5 = 'OXY'; % my name?
    param_props.glodap_name = 'G2oxygen'; % glodap name
    param_props.woa_name = 'OXYGEN'; % woa folder
    param_props.edges = 0:5:500;
    param_props.units = '(umol kg^{-1})';
    param_props.long_param_name = 'Dissolved Oxygen Amount Content';
    param_props.dec_points = 2;
    param_props.cmap = cmocean('ice');
elseif strcmp(p,'no3')
    param_props.p1 = 'NO3'; % directory name
    param_props.p2 = 'nitrate'; % for saving figures?
    param_props.p3 = '[NO_{3}]'; % axis labels
    param_props.p4 = 'NITRATE'; % Argo name
    param_props.p5 = 'NIT'; % my name
    param_props.p6 = 'G2nitrate'; % glodap name
    param_props.p7 = 'NITRATE'; % woa folder
    param_props.edges = 0:0.5:50;
    param_props.units = '(umol kg^{-1})';
    param_props.long_param_name = 'Dissolved Nitrate Amount Content';
    param_props.dec_points = 2;
    param_props.cmap = flipud(cmocean('speed'));
elseif strcmp(p,'ph')
    param_props.p1 = 'pH'; % directory name
    param_props.p2 = 'ph'; % for saving figures?
    param_props.p3 = 'pH_{T}'; % axis labels
    param_props.p4 = 'PH_IN_SITU_TOTAL'; % Argo name
    param_props.p5 = 'PH'; % my name
    param_props.p6 = 'G2phtsinsitutp'; % glodap name
    param_props.p7 = '';
    param_props.edges = 7.2:0.01:8.2;
    param_props.units = '';
    param_props.long_param_name = 'In situ pH on the Total Scale';
    param_props.dec_points = 4;
    param_props.cmap = cmocean('matter');
end
