% save variables within parfor loop
% syntax = parsave(filename,var1,'var2name',var2,'var2name',...)

function parsave(fname,varargin)
    if rem(length(varargin),2)
        error('Input must be provided in the form of "variable, variable name"');
    else
        for i = 1:2:length(varargin)-1
            if i == 1
                data.(varargin{i+1}) = varargin{i};
                try
                    save(fname,'-struct','data','-v7');
                catch
                    save(fname,'-struct','data','-v7.3');
                end
            else
                data.(varargin{i+1}) = varargin{i};
                save(fname,'-struct','data','-append');
            end
        end 
    end
end