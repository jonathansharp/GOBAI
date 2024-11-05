% This function checks whether a variable exists in a given netcdf file
%
% Input:
%   fname - netcdf file name
%   var - variable name to check
% Output:
%   yn - logical whether the variable exists in the file
%

function yn = nc_var_exist(fname,var)
    finfo = ncinfo(fname);
    varNames = {finfo.Variables.Name};
    varMatch = strncmpi(varNames,var,1);
    yn = any(varMatch);
end