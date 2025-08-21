%% Create basin shapes
Poly.Atl = [-60,0;-79,9.4;-81,8.4;-100,22;-100,45;-6,45;-6,35;4,15;25,0;22,-35;-68,-35;-60,0];
Poly.Pac1 = [104,0;104,70;181,70;181,0;181,-35;145,-35;131,-30;131,0;104,0];
Poly.Pac2 = [-180,0;-180,70;-150,70;-150,67;-120,67;-100,22;-81,8.4;-79,9.4;-60,0;-68,-35;-180,-35;-180,0];
Poly.Ind = [22,-35;25,10;38,35;104,35;104,0;131,0;131,-30;116,-35;22,-35];
Poly.Arc = [-180,64;-180,90;181,90;181,67;90,67;0,50;0,40;-6,40;-6,35;-90,35;-120,64;-180,64];
Poly.Med = [-6.5,40;0,40;0,45;20,47;38,35;34,30;-5,30;-6.5,40];
Poly.NSou = [-180,-60;-180,-25;181,-25;181,-60;-180,-60];
Poly.SSou = [-180,-90;-180,-50;181,-50;181,-90;-180,-90];

% Southern Ocean shapes just for plots
Poly.NSou_plot1 = [-180,-60;-180,-25;0,-25;0,-60;-180,-60];
Poly.NSou_plot2 = [0,-60;0,-25;180,-25;180,-60;0,-60];
Poly.SSou_plot1 = [-180,-90;-180,-50;0,-50;0,-90;-180,-90];
Poly.SSou_plot2 = [0,-90;0,-50;180,-50;180,-90;0,-90];

%% Plot shapes
figure; worldmap([-90 90],[0 360]);
hold on;
setm(gca,'ffacecolor','w');
setm(gca,'fontsize',12);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land,'FaceColor',rgb('grey'));
dt = 'polygon'; fa = 0.2;
geoshow(Poly.Atl(:,2),Poly.Atl(:,1),'DisplayType',dt,'FaceColor','r','FaceAlpha',fa)
geoshow(Poly.Pac1(:,2),Poly.Pac1(:,1),'DisplayType',dt,'FaceColor','g','FaceAlpha',fa)
geoshow(Poly.Pac2(:,2),Poly.Pac2(:,1),'DisplayType',dt,'FaceColor','g','FaceAlpha',fa)
geoshow(Poly.Ind(:,2),Poly.Ind(:,1),'DisplayType',dt,'FaceColor','k','FaceAlpha',fa)
geoshow(Poly.Arc(:,2),Poly.Arc(:,1),'DisplayType',dt,'FaceColor','c','FaceAlpha',fa)
geoshow(Poly.Med(:,2),Poly.Med(:,1),'DisplayType',dt,'FaceColor','m','FaceAlpha',fa)
geoshow(Poly.NSou_plot1(:,2),Poly.NSou_plot1(:,1),'DisplayType',dt,'FaceColor','b','FaceAlpha',fa)
geoshow(Poly.NSou_plot2(:,2),Poly.NSou_plot2(:,1),'DisplayType',dt,'FaceColor','b','FaceAlpha',fa)
geoshow(Poly.SSou_plot1(:,2),Poly.SSou_plot1(:,1),'DisplayType',dt,'FaceColor','y','FaceAlpha',fa)
geoshow(Poly.SSou_plot2(:,2),Poly.SSou_plot2(:,1),'DisplayType',dt,'FaceColor','y','FaceAlpha',fa)
mlabel off; plabel off;
