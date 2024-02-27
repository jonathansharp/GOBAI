%% Plot cluster probabilities

% set pressures
pressure = 100;

% set up parallel pool
tic; parpool; fprintf('Pool initiation:'); toc;

parfor cluster_idx = 1:num_clusters

    % load probabilities
    load(['Data/GMM_' base_grid '_' num2str(num_clusters) '/c' ...
        num2str(cluster_idx) '/m1w1],'GMM_probs');

    for depth_idx = 1:length(GMM.pressure)

    % establish figure
    h=figure('visible','on','Position',[100 100 800 400]);
    set(h,'color','white');
    axis tight manual

    filename = ['Figures/Clusters_' num2str(num_clusters) ...
        '/cluster_probability_animation_c' num2str(cluster_idx) ...
        '_' num2str(GMM_probs.pressure(depth_idx)) 'dbar.gif'];
    n=1;
    for m = 1:5%length(GMM_probs.time)
        worldmap([-90 90],[20 380]);
        title(extractAfter(datestr(datenum(2004,0.5+double(GMM_probs.time(m)),1)),'-'));
        pcolorm(double(GMM_probs.latitude),double(GMM_probs.longitude),...
            GMM_cluster_probs(:,:,depth_idx,m)');
        colormap(flipud(gray(100)));
        plot_land('map');
        clim([0 1]);
        c=colorbar;
        c.Label.String = ['Cluster #' num2str(cluster_idx) ' Probability'];
        c.TickLength = 0;
        mlabel off; plabel off;
        % capture frame
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % write to file
        if n == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
        n=n+1;
    end
    close
    
    % clean up
    clear GMM_probs filename h cluster_idx depth_idx n m c frame im imind cm

end