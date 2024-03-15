% reformat longitude vector and associated matrix to the imposed structure

function [lonr,zr] = reformat_lon(lon,z,lonmin)

% reformat longitude
idx = find(abs(lon - lonmin) == min(abs(lon - lonmin))); idx = idx(end);
lonr = [lon(idx:end) lon(1:idx-1)+360];

% reformat gridded variable
londim = find(length(lon)==size(z));
if londim == 1
    zr = [z(idx:end,:);z(1:idx-1,:)];
elseif londim == 2
    zr = [z(:,idx:end) z(:,1:idx-1)];
else
    error('Variable matrix dimension does not match longitude dimension.')
end