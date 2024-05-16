% save variables within parfor loop
% syntax = parsave(filename,var1,'var2name',var2,'var2name',...)

function parsave(fname,varargin)

    if rem(length(varargin),2)
        error('Input must be provided in the form of "variable, variable name"');
    else
        % first try saving as normal (v7) file (must be < 2GB)
        for i = 1:2:length(varargin)-1
            data.(varargin{i+1}) = varargin{i};
            if i == 1
                save(fname,'-struct','data','-v7');
            else
                save(fname,'-struct','data','-append');
            end
        end
        % larger files must be saved as v7.3
        % if file doesn't exist
        if ~exist(fname,'file')
            for i = 1:2:length(varargin)-1
                data.(varargin{i+1}) = varargin{i};
                if i == 1
                    save(fname,'-struct','data','-v7.3');
                else
                    save(fname,'-struct','data','-append');
                end
            end
        % or if file is empty
        else
            s = whos('-file',fname);
            if isempty(s)
                for i = 1:2:length(varargin)-1
                    data.(varargin{i+1}) = varargin{i};
                    if i == 1
                        save(fname,'-struct','data','-v7.3');
                    else
                        save(fname,'-struct','data','-append');
                    end
                end
            end
        end
    end
    
end
