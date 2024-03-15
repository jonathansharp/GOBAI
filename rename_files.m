% change .sh file names (o2 to no3)
files = dir;
for f = 1:length(files)
    if contains(files(f).name,'.sh')
        old_file_name = files(f).name;
        new_file_name = replace(files(f).name,'o2','no3');
        movefile(old_file_name,new_file_name);
    end
end