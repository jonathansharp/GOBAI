% plot Figure 2 for Sharp et al. (in prep) GOBAI-O2 and GOBAI-NO3 Paper

figure('Position',[100 100 800 850],'Visible','off');
dot_size = 2;
clrs = jet(20);

%% load interpolated float and glodap data for o2
load('NO3/Data/GMM/all_data_clusters_20_fgo_Jun-2026_D_A.mat');
clusts = fieldnames(all_data_clusters);
for c = 1:length(clusts)-1
    all_data_clusters = rmfield(all_data_clusters,['c' num2str(c)]);
end
load('NO3/Data/processed_all_no3_data_fgo_Jun-2026_D_A.mat');

%% define pressures
press = [10;200;400;1250];
height = 0.33; width = 0.48;
pos = [0.01 0.68 width height;...
       0.51 0.68 width height;...
       0.01 0.42 width height;...
       0.51 0.42 width height;...
       0.08 0.08 0.9 0.35];

for p = 1:length(press)

    % next tile
    % nexttile(tiles(d),[3 1]); hold on;
    idx = all_data.pressure == press(p);

    %% Map for full period

    % plot map
    ax.(['ax' num2str(p)]) = axes('Position',pos(p,:));
    m_proj('robinson','lon',[20 380]);
    all_data.longitude(all_data.longitude<20) = ...
        all_data.longitude(all_data.longitude<20) + 360;
    m_scatter(all_data.longitude(idx),all_data.latitude(idx),dot_size,...
        all_data_clusters.clusters(idx),'filled');

    % set properties
    title([num2str(press(p)) ' dbar'],'FontSize',12);
    m_coast('patch',rgb('grey'));
    m_grid('linestyle','none','xticklabels',[],...
        'yticklabels',[],'ytick',-90:30:90);
    clim([0.5 20.5]); colormap(clrs);

    % add to main panel
    %set(gca,'Parent',panel_main)

%     if d == 3 || d ==4
%         c = colorbar('location','southoutside');
%         colormap(jet(20)); caxis([0.5 20.5]);
%         c.Label.String = 'Clusters';
%         c.TickLength = 0;
%         c.Label.FontSize = 12;
%     end

end

% define axis
ax.ax5 = axes('Position',pos(5,:));

% define cluster percentages
unique_pressure = unique(all_data.pressure);
clusters = 1:length(unique(all_data_clusters.clusters));
cluster_per = nan(length(unique_pressure),length(clusters));
for d = 1:length(unique_pressure)
    for c = 1:length(clusters)
        idx_pres = all_data.pressure == unique_pressure(d);
        idx_both = idx_pres & all_data_clusters.clusters == clusters(c);
        cluster_per(d,c) = 100*(sum(idx_both)/sum(idx_pres));
    end
end

% plot percentages
set(gca,'TickLength',[0 0]);
b=bar(categorical(unique_pressure),cluster_per,'stacked');
bar_w = 0.9; 
for k=1:numel(b)
    b(k).FaceColor = clrs(k,:);
    b(k).BarWidth = bar_w;
end
ylim([0 100]);
ylabel('Cluster %','FontSize',12);
xlabel('Pressure','FontSize',12);
legend(fliplr(string(clusters)),'Location','northoutside','NumColumns',10,...
    'Direction','normal');

% add text
annotation('textbox', [0 0 1 0.9], ...
    'String', 'Cluster Number', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Find categorical x-indices of the target pressures within unique_pressure
[~, highlight_idx] = intersect(unique_pressure, press);
% outline relevant bars
hold(ax.ax5, 'on');
for i = 1:length(highlight_idx)
    x_center = highlight_idx(i);
    x_left   = x_center - (bar_w / 2);
    % Draw a highlighted border around the entire stacked bar (0 to 100%)
    rectangle(ax.ax5, ...
        'Position', [x_left, 0, bar_w, 100], ...
        'EdgeColor', 'k', ...         % Outline color (Black)
        'LineWidth', 3.0, ...         % Thickness of outline
        'LineStyle', '-');            % Use '--' if you want dashed highlights
end
hold(ax.ax5, 'off');

export_fig(gcf,'Figures/Figure2_no3_cluster_distribution.png','-transparent');
close
