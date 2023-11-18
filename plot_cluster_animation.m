%% Plot clusters
load(['Data/GMM_' num2str(num_clusters) '/GMM'],'GMM');

for depth_idx = 1:length(GMM.pressure)
    % establish figure
    h=figure('visible','on','Position',[100 100 800 400]);
    axis tight manual
    % create folder
    if ~isfolder(['Figures/Clusters_' num2str(num_clusters)])
        mkdir(['Figures/Clusters_' num2str(num_clusters)]);
    end
    % establish file name
    filename = ['Figures/Clusters_' num2str(num_clusters) ...
        '/cluster_animation_' num2str(GMM.pressure(depth_idx)) 'dbar.gif'];
    % set counter ro 1
    n=1;
    % plot clusters each month
    for m = 1:length(GMM.time)
        worldmap([-90 90],[20 380]);
        title(extractAfter(datestr(datenum(2004,0.5+double(GMM.time(m)),1)),'-'));
        pcolorm(double(GMM.latitude),double(GMM.longitude),...
            double(GMM.clusters(:,:,depth_idx,m))');
        colormap([1,1,1;flipud(jet(10))]); % white then jet
        plot_land('map');
        clim([-0.5 10.5]);
        c=colorbar;
        c.Limits = [0.5 10.5];
        c.Label.String = 'Cluster';
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
    clear filename h depth_idx n m c frame im imind cm

end
