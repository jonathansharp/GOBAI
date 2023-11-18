% plot land areas
function plot_land(type)

if strcmp(type,'xy')

    land = shaperead('landareas', 'UseGeoCoords',false);
    land(2).X(land(2).X>-150) = land(2).X(land(2).X>-150)-360;
    land(163).X(land(163).X>-150) = land(163).X(land(163).X>-150)-360;
    for l = 1:length(land)
        land(l).X = land(l).X+360;
    end
    mapshow(land,'FaceColor',rgb('grey'));

elseif strcmp(type,'map')

    land = shaperead('landareas', 'UseGeoCoords',true);
    geoshow(land,'FaceColor',rgb('grey'));

end

clear land