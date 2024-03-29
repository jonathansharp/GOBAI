% process input parameter name to GOBAI functions

function [p1,p2,p3,p4,p5,p6,p7,edges] = param_name(p)

if strcmp(p,'o2')
    p1 = 'O2'; % directory name
    p2 = 'oxygen'; % for saving figures?
    p3 = '[O_{2}]'; % axis labels
    p4 = 'DOXY'; % Argo name
    p5 = 'OXY'; % my name
    p6 = 'G2oxygen'; % glodap name
    p7 = 'OXYGEN'; % woa folder
    edges = 0:5:500;
elseif strcmp(p,'no3')
    p1 = 'NO3'; % directory name
    p2 = 'nitrate'; % for saving figures?
    p3 = '[NO_{3}]'; % axis labels
    p4 = 'NITRATE'; % Argo name
    p5 = 'NIT'; % my name
    p6 = 'G2nitrate'; % glodap name
    p7 = 'NITRATE'; % woa folder
    edges = 0:0.4:40;
elseif strcmp(p,'ph')
    p1 = 'pH'; % directory name
    p2 = 'ph'; % for saving figures?
    p3 = 'pH_{T}'; % axis labels
    p4 = 'PH_IN_SITU_TOTAL'; % Argo name
    p5 = 'PH'; % my name
    p6 = 'G2phtsinsitutp'; % glodap name
    edges = 7.2:0.01:8.2;
end