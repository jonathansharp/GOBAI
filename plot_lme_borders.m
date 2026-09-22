function plot_lme_borders(region,lme_shape,lme_idx)

for n = 1:length(region)
    tmp_lon = convert_lon(lme_shape(lme_idx.(region{n})).X');
    tmp_lat = lme_shape(lme_idx.(region{n})).Y';
    if n == 3 % remove line at 180 for AI
      idx=tmp_lon<180.01&tmp_lon>179.99;
      plotm(tmp_lat(~idx),tmp_lon(~idx),'k','linewidth',1);
    elseif n == 4 % remove line at 180 for EBS
      tmp_lon(1) = NaN; tmp_lat(1) = NaN;
      tmp_lon(3315) = NaN; tmp_lat(3315) = NaN;
      plotm(tmp_lat,tmp_lon,'k','linewidth',1);
    elseif n == 11 % remove line at 180 for EBS
      tmp_lon(11418:11429) = NaN; tmp_lat(11418:11429) = NaN;
      tmp_lon(13358:13370) = NaN; tmp_lat(13358:13370) = NaN;
      plotm(tmp_lat,tmp_lon,'k','linewidth',1);
    else
        m_plot(tmp_lat,tmp_lon,'k','linewidth',1);
    end
end