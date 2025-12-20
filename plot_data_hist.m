% plot data histogram
function plot_data_hist(param_props,file_date,float_file_ext,...
    include_float,include_glodap,include_ctd,y1,y2)

% define dataset extensions
if include_float == 1; float_ext = 'f'; else float_ext = ''; end
if include_glodap == 1; glodap_ext = 'g'; else glodap_ext = ''; end
if include_ctd == 1; ctd_ext = 'w'; else ctd_ext = ''; end

% load data
load([param_props.dir_name '/Data/processed_all_' param_props.file_name '_data_' ...
    float_ext glodap_ext ctd_ext '_' file_date float_file_ext '.mat'],'all_data');

% establish figure
figure('visible','on');
set(gcf,'Position',[100 100 1200 400]);
set(gca,'FontSize',14);

% index to unique profiles
[~,prof_idx] = unique(all_data.id);
vars = fieldnames(all_data);
for v = 1:length(vars)
    if ndims(all_data.(vars{v})) < 3
        all_data.(vars{v}) = all_data.(vars{v})(prof_idx);
    end
end

% index to each dataset type
if include_ctd; y_ctd = all_data.year(all_data.type==3); end
if include_glodap; y_osd = all_data.year(all_data.type==2); end
if include_float; y_flt = all_data.year(all_data.type==1); end

% count profiles in each year
if include_ctd; counts_ctd = histc(y_ctd,y1:y2); end
if include_glodap; counts_osd = histc(y_osd,y1:y2); end
if include_float; counts_flt = histc(y_flt,y1:y2); end

% plot histogram
if include_ctd
    bar(start_year:end_year,[counts_osd counts_ctd counts_flt],'stacked');
    legend({'GLODAP' 'CTD' 'Argo Float'},'Location','northwest','FontSize',14);
else
    bar(y1:y2,[counts_osd counts_flt],'stacked');
    legend({'GLODAP' 'Argo Float'},'Location','northwest','FontSize',14);
end
ylabel('Number of Profiles');

% save figure
export_fig(['O2/Figures/dataset_histogram_' num2str(y1) '_' ...
    num2str(y2) '.png'],'-transparent','-silent');
