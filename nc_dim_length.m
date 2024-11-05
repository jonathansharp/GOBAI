function l = nc_dim_length(fname,dim)
    finfo = ncinfo(fname);
    dimNames = {finfo.Dimensions.Name};
    dimMatch = strncmpi(dimNames,dim,1);
    dimInfo = finfo.Dimensions(dimMatch);
    l = dimInfo.Length;
end