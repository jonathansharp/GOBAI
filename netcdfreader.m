function output = netcdfreader(path)

info = ncinfo(path);
for n = 1:size(info.Variables,2)
    if length(info.Variables(n).Dimensions) > 1
        output.(info.Variables(n).Name) = ...
            ncread(path,info.Variables(n).Name);
    else
        output.(info.Variables(n).Name) = ...
            ncread(path,info.Variables(n).Name);
    end
end

end