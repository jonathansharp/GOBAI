% function to load GOBAI

function gobai = load_gobai(var,ver,path,fields)

% download dimensions
gobai.lon = ncread([path 'GOBAI-' var '-' ver '.nc'],'lon');
gobai.lat = ncread([path 'GOBAI-' var '-' ver '.nc'],'lat');
gobai.pres = ncread([path 'GOBAI-' var '-' ver '.nc'],'pres');
try gobai.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'time'); end
try gobai.time = ncread([path 'GOBAI-' var '-' ver '.nc'],'year'); end

% download variables
for v = 1:length(fields)
    gobai.(fields{v}) = ncread([path 'GOBAI-' var '-' ver '.nc'],fields{v});
end

end