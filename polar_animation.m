function tree_movie_temp_100(TreeSetUp)



% min_layer=0;
% max_layer=2000;
% syear=2021;
% fyear=syear-1;
% 
% 
% path_figs='C:\data\OHCA\figs\tree_paper\'
% year_of_oco_pub=2022;
% slope_min_year=1993;
% file_name='argo_2021_01_01_QC'
% % this is the time range of the maps that are to be saved and outputted
% max_year_maps_out=2021;
% min_year_maps_out=1993;
cold_to_hot_colormap=diverging_map([0:1/200:1],[20 43 140]/255,[204 0 51]/255);
% min_depth and max_depth must be a layer bounds

nbasins_use=TreeSetUp.nbasins_use;
file_name=TreeSetUp.file_name;
var_type=TreeSetUp.var_type;

file_name_season=TreeSetUp.file_name_season;
file_name_season_anom=TreeSetUp.file_name_season_anom;
file_WOD_suf=TreeSetUp.file_WOD_suf;
file_path_hdata=TreeSetUp.file_path_hdata;
fname_nc_season=TreeSetUp.fname_nc_season;
fname_nc=TreeSetUp.fname_nc;

path_ERDDAP=TreeSetUp.path_ERDDAP;

% path_Fig_data=TreeSetUp.path_Fig_data;

tree_prefix=TreeSetUp.tree_prefix;
tree_model_file_name_season=TreeSetUp.tree_model_file_name_season;
tree_model_file_name_yearly=TreeSetUp.tree_model_file_name_yearly;
tree_model_file_name_all_year=TreeSetUp.tree_model_file_name_all_year;
tree_model_file_name_combined=TreeSetUp.tree_model_file_name_combined;
tree_model_file_name_combined_withcycle=TreeSetUp.tree_model_file_name_combined_withcycle;

% 
% file_name_season=[file_name,'_seasonal'];
% file_name_season_anom=[file_name_season,'_anom'];

path_oisst=TreeSetUp.path_oisst;
path_OHCA_data_out=TreeSetUp.path_OHCA_data_out;
path_OHCA_data_in=TreeSetUp.path_OHCA_data_in;
path_ssh=TreeSetUp.path_ssh;

path_tree=TreeSetUp.path_tree;
path_new_tree_season=TreeSetUp.path_new_tree_season;
path_new_tree_yearly=TreeSetUp.path_new_tree_yearly;
path_new_tree_all_year=TreeSetUp.path_new_tree_all_year;
path_new_tree_combined=TreeSetUp.path_new_tree_combined;
path_new_tree_combined_withcycle=TreeSetUp.path_new_tree_combined_withcycle;

path_Figs=TreeSetUp.path_Figs;
path_tree_junk=TreeSetUp.path_tree_junk;
path_curve=TreeSetUp.path_curve;

layer_bounds=TreeSetUp.layer_bounds;

start_year=TreeSetUp.start_year;
end_year=TreeSetUp.end_year;

start_year_mean=TreeSetUp.start_year_mean;
end_year_mean=TreeSetUp.end_year_mean;
max_year_fit=TreeSetUp.max_year_fit;
min_year_fit=TreeSetUp.min_year_fit;
center_year=TreeSetUp.center_year;

start_yearly_maps=TreeSetUp.start_yearly_maps;
end_yearly_maps=TreeSetUp.end_yearly_maps;

start_all_year=TreeSetUp.start_all_year;
end_all_year=TreeSetUp.end_all_year;

start_year_trans=TreeSetUp.start_year_trans;
end_year_trans=TreeSetUp.end_year_trans;

% OUTOUT_type=TreeSetUp.OUTOUT_type;%%


function polar_animation(param,vrs,year)
    movie_name=['temp_100m_' vrs '_' year '_SO.gif'];
    filename=[path_Figs,movie_name];

    
end





movie_name=['temp_100m_',num2str(floor(end_year)),'_2015_newg.gif'];
filename=[path_Figs,movie_name];
% v=VideoWriter([path_Figs,movie_name],'MPEG-4');
% v.FrameRate=4;
% open(v)
% nframe=length(time_aviso);
% 
% 
% year_aviso=floor(time_aviso);
% aviso_day=1+(time_aviso-year_aviso).*yeardays(year_aviso);
% sday=round(aviso_day+datenum(year_aviso,1,1)-1);
% 
% days_since_1950=sday-datenum(1950,1,1);%NEED TO CHECK THIS

subdir='yearly_withcycle';
start_year_file=start_year;
end_year_file=end_year;
path_new_tree=path_new_tree_combined_withcycle;
tree_model=tree_model_file_name_combined_withcycle;
% tree_file_name=tree_file_name_in;
%%

path_nc_erddap=[path_ERDDAP,'netcdf\',tree_prefix,'\',subdir,'\'];
if var_type=='s'
     file_prefix='RFROMV2_SAL_';
elseif var_type=='t'
     file_prefix='RFROMV2_TEMP_';
elseif var_type=='h'
     file_prefix='RFROMV2_OHC_';
end
%%


start_year_movie=2015;
start_month_movie=1;
end_year_movie=floor(end_year);
end_year_movie=2015;
end_month_movie=12;

nyears_movie=end_year_movie-start_year_movie+1;

years_movie=repelem(start_year_movie:end_year_movie,12);
months_movie=repmat(1:12,[1 nyears_movie]);

months_movie=months_movie(start_month_movie:end-12+end_month_movie);
years_movie=years_movie(start_month_movie:end-12+end_month_movie);

nfiles=length(months_movie);





for ifile=1:nfiles

    %load in data
    year_load=years_movie(ifile);
    month_load=months_movie(ifile);

     if month_load>=10
          file_name_nc= [path_nc_erddap,file_prefix,num2str(year_load),'_',num2str(month_load),'.nc'];
       else
          file_name_nc= [path_nc_erddap,file_prefix,num2str(year_load),'_0',num2str(month_load),'.nc'];
     end

     lon_tpx=double(ncread(file_name_nc,'longitude'));
     lat_tpx=double(ncread(file_name_nc,'latitude'));

     time=double(ncread(file_name_nc,'time'));

     temp_total=ncread(file_name_nc,'ocean_temperature');

     ntime=length(time);





    for itime=1:ntime


        figure(1)
        clf
            set(gcf,'color','white');
            m_proj('Equidistant cylindrical','long',[30 390],'lat',[-90 90]);
            
            
            
            
            corrhc=double(squeeze((temp_total(:,:,13,itime))));
            date_name=datestr(time(itime)+datenum(1950,1,1));

           
            
            
            lon=lon_tpx';
            lat=lat_tpx';
            scale_fac=1;
            
            % put into the proper coordinates




%             del_val=2/scale_fac;
%             del_cont=.5/scale_fac;
%             min_val=-3./scale_fac;
%             max_val=3./scale_fac;
%             del_val=.5./scale_fac;
%             del_cont=.25./scale_fac;

             min_val=-2;
             max_val=28.;
            cont_level=[min_val:.2:5 6:max_val];
            val_lev=[min_val:1:5 7:2:max_val];
            del_cont=abs(cont_level(end)-cont_level(end-1));
            ii=find(lon<30);
            jj=find(lon>=30);
            lon=[lon(jj),lon(ii)+360];
            corrhc=[corrhc(jj,:);corrhc(ii,:)];
            %colormap jet(256)
            
            colormap(cold_to_hot_colormap) 
            
            
            m_proj('Equidistant cylindrical','long',[30 390],'lat',[-90 90]);
            %colormap(fresh_to_salty_colormap) 
            [cs1,h1]=m_contourf(lon,lat,corrhc',[-1000,min_val:del_cont:max_val]);
            % % 'cat'
        %     save 'OHCA_2019_tpx.mat' lon lat corrhc
            hold on
            set(h1,'linecolor','none')
            hold on
            m_grid('tickdir','out','xtick',[30:60:390],'ytick',[-90:30:90],'linestyle','none');
            
            
            %[cs12,h12]=m_contour(lon_cont,lat_cont,sal_cont',[32 32],'k');
            a=gca;
            hold on
            m_coast('patch',[0 0 0]);
            t1=m_text(30,100,[' ', date_name],'fontsize',12);
            %t3=m_text(50,45,[num2str(depth_top_plot),'-',num2str(depth_bot_plot)],'fontsize',12);
             t2=m_text(170,-150,'         (C^{\circ})','fontsize',10);
           
            caxis([min_val max_val-del_cont])
            
            hold off
            
            ja=axes('pos',[.262 .90-.0475 .51 .01]);
            %colormap jet(256)
            
            colormap(cold_to_hot_colormap) 
            
            [cs,h]=contourf([min_val:.01:max_val],[0 1],[1 1]'*[min_val:.01:max_val],[min_val:del_cont:max_val]);
            set(h,'edgecolor','none')
            set(ja,'tickdir','out','xaxisl','top','xtick',val_lev,'ytick',[])
            caxis([min_val max_val-del_cont])
            japos=get(ja,'pos');
            set(ja,'XAxisLocation','bottom','pos',japos-[.075 .7 -.15 -.015])
        
            frame=getframe(gcf);
            im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if itime== 1 && ifile==1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.25);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',.25,'WriteMode','append');
        end


%             writeVideo(v,frame)
        end


end
% close(v)

