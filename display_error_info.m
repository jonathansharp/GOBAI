% function to display relevant information from an error message

function display_error_info(ME)

% display general error message
disp(ME);
% determine levels to error message
levs = length(ME.stack);
% determine functions and line numbers of error
for l = levs:-1:1
    disp(['Function: ' ME.stack(l,1).name]);
    disp(['Line: ' num2str(ME.stack(l,1).line)]);
end

end