% Runs all scripts to plot GOBAI-O2 RFR gridded fields

try
    dir_base = create_dir_base('RFR',{base_grid;num_clusters;file_date;...
        float_file_ext;numtrees;minLeafSize});
    %plot_gobai_mean;
    plot_gobai_animation('o2',dir_base,base_grid,num_clusters,'RFR')
    disp('success!');
catch ME
    display_error_info(ME);
end
exit