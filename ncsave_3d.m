function ncsave_3d(fpath,x,y,z,vr)

% delete file if it exists
if isfile(fpath); delete(fpath); end

% write x
nccreate(fpath,x{1},'Dimensions',{x{1},length(x{2})});
ncwrite(fpath,x{1},x{2});
ncwriteatt(fpath,x{1},'long_name',x{3});
ncwriteatt(fpath,x{1},'units',x{4});

% write y
nccreate(fpath,y{1},'Dimensions',{y{1},length(y{2})});
ncwrite(fpath,y{1},y{2});
ncwriteatt(fpath,y{1},'long_name',y{3});
ncwriteatt(fpath,y{1},'units',y{4});

% write z
nccreate(fpath,z{1},'Dimensions',{z{1},length(z{2})});
ncwrite(fpath,z{1},z{2});
ncwriteatt(fpath,z{1},'long_name',z{3});
ncwriteatt(fpath,z{1},'units',z{4});

% write data variable
nccreate(fpath,vr{1},'Dimensions',{x{1},length(x{2}),y{1},length(y{2}),z{1},length(z{2})});
ncwrite(fpath,vr{1},vr{2});
ncwriteatt(fpath,vr{1},'long_name',vr{3});
ncwriteatt(fpath,vr{1},'units',vr{4});
ncwriteatt(fpath,vr{1},'fillvalue','NaN');

end
